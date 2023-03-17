version 1.1


workflow intersect_pav {
    input {
      File ordered_vcf_shards_list
      File ordered_sample_names
      File background_bed
      File intersect_script
      File svpoplib
    }

    meta {
        workflow_description: "Intersects with HGSVC+HPRC"
    }
    parameter_meta {
        ordered_vcf_shards_list: "List of VCFs to process"
        background_bed: "Bed file to intersect the background with"
        intersect_script: "Script to run the intersects"
        ordered_sample_names: "Sample names associated with VCFs"
        svpoplib: "SVPOP Library"
    }

    Array[File] input_vcfs = read_lines(ordered_vcf_shards_list)
    Array[File] input_sample_names = read_lines(ordered_sample_names)
    
    scatter (vcf_shard_pair in zip(input_vcfs, input_sample_names)) {
        call ProcessPavVcf {
            input:
                vcfIn = vcf_shard_pair.left,
                sample = vcf_shard_pair.right,
        }
        call IntersectWithBackground {
            input:
                bed = ProcessPavVcf.bed,
                sample = vcf_shard_pair.right,
                threads = "4",
                intersect_py = intersect_script,
                background = background_bed,
                svpoplib_tgz = svpoplib,
        }

    }
    output {
        Array[File] sample_intersect_beds = IntersectWithBackground.int_bed
        Array[File] sample_pav_beds = ProcessPavVcf.bed
    }
}

 
task ProcessPavVcf {
  input {
    File vcfIn
    String sample 
  }

  command <<<
    set -e
    echo -e "#CHROM\tPOS\tEND\tREF\tALT\tGT\tID\tSVTYPE" > ~{sample + "_PAV.bed"}
    bcftools query -f "%CHROM\t%POS\t%END\t%REF\t%ALT\t[%GT]\t%INFO/ID\t%INFO/SVTYPE\n" ~{vcfIn} | grep -v SNV >> ~{sample + "_PAV.bed"}
  >>>

  runtime {
    memory: "8 GiB"
    cpu: "1"
    disks: "local-disk 100 HDD"
    docker: "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.1.17"
  }

  output {
    File bed = "~{sample}_PAV.bed"
  }
}
 


task IntersectWithBackground {

    input {
        File bed
        File background
        File intersect_py
        File svpoplib_tgz
        String threads
        String sample
    }
    
    command <<<
        set -e
        tar zxvf ~{svpoplib_tgz}
        python ~{intersect_py} -a ~{bed} -b ~{background} -o ~{sample + "_bg_int.bed"} -t ~{threads}
    >>>

    runtime {
        memory: "8 GiB"
        cpu: threads
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.1.17"
    }

    output {
        File int_bed = "~{sample}_bg_int.bed"
    }    

}
