version 1.0

workflow ValidateVariants {
    input {
        String sample_id

        File cram
        File crai

        File ref_fa_uncompressed
        File ref_fai_uncompressed

        File counts_jf
        File bed

        File? sample_map

        File locityper_db_tar_gz

        Int locityper_n_cpu = 32
        Int locityper_max_gts = 100000
    }

    # call GunzipReference { input: ref_gz = ref_fa_with_alt }

    if (defined(sample_map)) {
        call SubsetBed {
            input:
                bed = bed,
                sample_map = select_first([sample_map]),
                sample_id = sample_id
        }
    }

    call SplitBed {
        input:
            bed = select_first([SubsetBed.subset_bed, bed]),
            N = 2000
    }

    scatter (split_bed in SplitBed.split_beds) {
        call LocityperPreprocessAndGenotype {
            input:
                sample_id = sample_id,
                cram = cram,
                crai = crai,
                reference = ref_fa_uncompressed,
                reference_index = ref_fai_uncompressed,
                db_targz = locityper_db_tar_gz,
                counts_file = counts_jf,
                bed = split_bed,
                locityper_n_cpu = locityper_n_cpu,
                locityper_max_gts = locityper_max_gts
        }
    }

    call CombineTarFiles { input: tar_files = LocityperPreprocessAndGenotype.genotype_tar, sample_id = sample_id }

    call Summarize { input: sample_id = sample_id, genotype_tar = CombineTarFiles.combined_tar_gz }

    output {
        File summary_csv = Summarize.summary_csv
        File results_tar_gz = CombineTarFiles.combined_tar_gz
    }
}

task GenerateDBFromVCF {
    input {
        File reference
        File reference_index
        File counts_jf
        File vcf
        File vcf_tbi
        File bed
    }

    Int disk_size = 1 + 10*ceil(size([reference, vcf, counts_jf, bed], "GiB"))
    String output_tar = "vcf_db.tar.gz"

    command <<<
        set -euxo pipefail

        nthreads=$(nproc)
        echo "using ${nthreads} threads"

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        locityper add \
            -@ ${nthreads} \
            -d vcf_db \
            -v ~{vcf} \
            -r reference.fa \
            -j ~{counts_jf} \
            -L ~{bed}

        echo "compressing DB"
        tar -czf ~{output_tar} vcf_db
        echo "done compressing DB"
    >>>

    output {
        File db_tar = output_tar
    }

    runtime {
        memory: "32 GB"
        cpu: "32"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
        docker: "eichlerlab/locityper:0.19.1"
    }
}

task GunzipReference {
    input {
        File ref_gz
    }

    Int disk_size_gb = 1 + 4*ceil(size(ref_gz, "GiB"))

    command <<<
        set -euxo pipefail

        gunzip -c ~{ref_gz} > reference.fa
        samtools faidx reference.fa
    >>>

    output {
        File ref_fa = "reference.fa"
        File ref_fai = "reference.fa.fai"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-locityper:1.0.0"
        memory: "2 GB"
        cpu: 2
        preemptible: 1
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task LocityperPreprocessAndGenotype {
    input {
        File cram
        File crai
        File reference
        File reference_index
        File counts_file
        File db_targz
        String sample_id
        # Array[String] locus_names

        File bed

        Int locityper_n_cpu
        Int locityper_max_gts
    }

    parameter_meta {
        cram: { localization_optional: true }
    }

    Int disk_size = 1 + 4*ceil(size([cram, crai, counts_file, reference, reference_index, db_targz, bed], "GiB"))
    Int locityper_mem_gb = ceil(4.0 * locityper_n_cpu)

    String output_tar = sample_id + ".locityper.tar.gz"

    command <<<
        set -euxo pipefail

        mv ~{reference} reference.fa
        mv ~{reference_index} reference.fa.fai
        mv ~{crai} full.cram.crai

        mkdir -p locityper_preproc

        cp ~{bed} new.bed
        echo -e "chr17\t72062001\t76562000" >> new.bed
        awk '{print $1 "\t" $2+5000 "\t" $3-5000 "\t" $4}' new.bed > extended.bed
        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
        samtools view -h -T reference.fa --regions-file extended.bed -b -o subset.bam -X ~{cram} full.cram.crai 
        samtools index subset.bam

        locityper preproc -a subset.bam \
            -r reference.fa \
            -j ~{counts_file} \
            -@ ~{locityper_n_cpu} \
            --technology illumina \
            -o locityper_preproc

        tar -xzf ~{db_targz}

        # Create out_dir and ensure it exists before parallel processing
        mkdir -p out_dir
        mkdir -p out_dir/loci

        process_single_locus() {
            line="$1"
            locus_name=$(echo "$line" | cut -f4)
            echo "Processing locus: ${locus_name}"
            
            # Ensure the locus-specific directory exists
            mkdir -p "out_dir/loci/${locus_name}"
            
            locityper genotype -a subset.bam \
                -r reference.fa \
                -d vcf_db \
                -p locityper_preproc \
                --subset-loci "${locus_name}" \
                -o out_dir
        }
        export -f process_single_locus

        cat ~{bed} | /usr/bin/parallel --line-buffer -j ~{locityper_n_cpu} process_single_locus {}

        find out_dir -type f -name "*.bam" -exec rm -f {} \;
        
        rm -f subset.bam
        rm -f new.bed

        tar -czf ~{output_tar} out_dir

    
    >>>

    output {
        File genotype_tar = output_tar
    }

    runtime {
        memory: "~{locityper_mem_gb} GB"
        cpu: locityper_n_cpu
        disks: "local-disk ~{disk_size} HDD"
        preemptible: 3
        docker: "docker://eichlerlab/locityper:1.4.5.0"
    }
}

task SplitBedNames {
    input {
        File bed
        Int N = 20
    }

    Int disk_size = 1 + ceil(size(bed, "GiB"))

    command <<<
        set -euxo pipefail

        cut -f4 ~{bed} | split -l ~{N} - split_part_ && wc -l split_part_*
    >>>

    output {
        Array[File] name_parts = glob("split_part_*")
    }

    runtime {
        memory: "1 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: "staphb/bcftools:1.22"
    }
}

task SubsetBed {
    input {
        File bed
        File sample_map
        String sample_id
    }

    Int disk_size = 1 + ceil(size([bed, sample_map], "GiB"))

    command <<<
        set -euxo pipefail

        grep '^~{sample_id}' ~{sample_map} | cut -f2 | sed 's/,/\n/g' > loci.txt
        cat loci.txt
        grep -f loci.txt ~{bed} > subset.bed
    >>>

    output {
        File subset_bed = "subset.bed"
    }

    runtime {
        memory: "1 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: "staphb/bcftools:1.22"
    }
}

task SplitBed {
    input {
        File bed
        Int N = 20
    }

    Int disk_size = 1 + ceil(size(bed, "GiB"))

    command <<<
        set -euxo pipefail

        cat ~{bed} | split -l ~{N} - split_part_ && wc -l split_part_*
    >>>

    output {
        Array[File] split_beds = glob("split_part_*")
    }

    runtime {
        memory: "1 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: "staphb/bcftools:1.22"
    }
}

task FilterNames {
    input {
        File bed
        Array[String] names_to_remove
    }

    Int disk_size = 1 + 2*ceil(size([bed], "GiB"))

    command <<<
        set -euxo pipefail

        grep -v -f ~{write_lines(names_to_remove)} ~{bed} | cut -f4 > filtered.txt
    >>>

    output {
        Array[String] filtered_names = read_lines("filtered.txt")
    }

    runtime {
        memory: "1 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: "staphb/bcftools:1.22"
    }
}

task FilterBed {
    input {
        File locus_names
        File bed
    }

    Int disk_size = 1 + 2*ceil(size([locus_names], "GiB"))

    command <<<
        set -euxo pipefail

        grep -f ~{locus_names} ~{bed} > filtered.bed
    >>>

    output {
        File filtered_bed = "filtered.bed"
    }

    runtime {
        memory: "1 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: "staphb/bcftools:1.22"
    }
}

task CombineTarFiles {
    input {
        Array[File] tar_files
        String sample_id
    }

    Int disk_size = 1 + 10*ceil(size(tar_files, "GiB"))
    String combined_tar = sample_id + ".combined.tar.gz"

    command <<<
        set -euxo pipefail

        # Create a temporary directory to extract all tar files
        mkdir -p combined_temp
        
        # Extract all tar files into the temporary directory
        for tar_file in ~{sep=" " tar_files}; do
            echo "Extracting $tar_file"
            tar -xzf "$tar_file" -C combined_temp
        done
        
        # Create a new combined tar file from the extracted contents
        echo "Creating combined tar file: ~{combined_tar}"
        tar -czf ~{combined_tar} -C combined_temp .
        
        echo "Combined tar file created successfully"
        ls -lh ~{combined_tar}
        
        # Clean up
        rm -rf combined_temp
    >>>

    output {
        File combined_tar_gz = combined_tar
    }

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.122"
    }
}

task Summarize {
    input {
        String sample_id
        File genotype_tar
    }

    Int disk_size = 1 + 10*ceil(size(genotype_tar, "GiB"))

    command <<<
        set -euxo pipefail

        tar -xzvf ~{genotype_tar}

        mv out_dir ~{sample_id}

        # Remove subdirectories that don't have res.json.gz
        for dir in ~{sample_id}/loci/*/; do
            if [ ! -f "${dir}/res.json.gz" ]; then
                echo "Removing directory ${dir} - no res.json.gz found"
                rm -rf "${dir}"
            fi
        done

        python3 /locityper/extra/into_csv.py -i ./~{sample_id} -o gts.csv

        grep -v '^#' gts.csv > gts.filtered.csv
    >>>

    output {
        File summary_csv = "gts.filtered.csv"
    }

    runtime {
        memory: "4 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.122"
    }
}
