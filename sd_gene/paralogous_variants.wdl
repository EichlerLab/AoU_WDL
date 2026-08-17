version 1.0

workflow ParalogousVariants {
  input {
    String sample_name
    File haplotig_fasta_hap1
    File haplotig_fasta_hap2
    File gene_seq
    File fixed_variants
  }

  meta {
    author: "Luyao Ren"
    workflow_description: "Detect paralog variants from assemblies"
  }

  Array[String] haps = ["hap1", "hap2"]

  scatter (hap in haps) {
    String sample_w_hap = sample_name + "_" + hap
    File haplotig_fasta = if (hap == "hap1") then haplotig_fasta_hap1 else haplotig_fasta_hap2
    
    call extract_sample_sequence {
      input:
        asm = haplotig_fasta,
        gene_seq = gene_seq,
        sample_id = sample_w_hap
    }
    call best_match {
      input:
        paf = extract_sample_sequence.paf,
        asm = haplotig_fasta,
        fixed_variants = fixed_variants,
        sample_id = sample_w_hap
    }
  }
  call call_variants {
    input:
      gene_seq = gene_seq,
      best_match_bed = best_match.best_match_bed,
      best_match_fa = best_match.best_match_fa,
      sample_id = sample_name
  }

  call normalize_merge_vcfs {
    input:
      gene_seq = gene_seq,
      raw_vcfs = call_variants.raw_vcfs,
      sample_id = sample_name
  }

  output {
    File final_bed = call_variants.bed
    File final_vcf = normalize_merge_vcfs.vcf
    File final_vcf_idx = normalize_merge_vcfs.vcf_idx
  }
}


task extract_sample_sequence {
  input {
    File asm
    File gene_seq
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    minimap2 -x asm20 -c --secondary=yes -p 0.3 -N 20 --eqx -t 4 -r 500 -K 200M ~{asm} ~{gene_seq} > ~{sample_id}.paf
  >>>

  output {
    File paf = "~{sample_id}.paf"
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          4,
    mem_gb:             32,
    disk_gb:            32,
    boot_disk_gb:       32,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/align-basics:0.2"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 select_first([runtime_attr.docker,            default_attr.docker])
  }
}

task best_match {
  input {
    File paf
    File asm
    File fixed_variants
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    rb stat --paf ~{paf} > ~{sample_id}.bed
    bedtools getfasta -fi ~{asm} -bed ~{sample_id}.bed -fo ~{sample_id}.fa
    seqkit locate -f ~{fixed_variants} ~{sample_id}.fa |  awk '
      NR==1 { next }
      {
        split($1, a, /[:\-]/)
        key = a[1] "\t" a[2] "\t" a[3]
        count[key]++
      }
      END {
        for (k in count) {
          if (count[k] >= 5)
            print k
        }
      }
      ' | awk -v sid="~{sample_id}" 'BEGIN{OFS="\t"} {print $0, sid}' > ~{sample_id}.best_match.bed
    bedtools getfasta -fi ~{asm} -bed ~{sample_id}.best_match.bed -fo ~{sample_id}.best_match.fa
  >>>

  output {
    File bed = "~{sample_id}.bed"
    File sample_fa = "~{sample_id}.fa"
    File best_match_bed = "~{sample_id}.best_match.bed"
    File best_match_fa = "~{sample_id}.best_match.fa"
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          1,
    mem_gb:             8,
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/assembly_eval:0.5"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 select_first([runtime_attr.docker,            default_attr.docker])
  }
}

task call_variants {
  input {
    Array[File] best_match_fa
    Array[File] best_match_bed
    File gene_seq
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    cat ~{sep=" " best_match_fa} > merged.fa
    cat ~{sep=" " best_match_bed} > ~{sample_id}.bed
    mkdir -p split_fa raw_vcfs
    seqkit split -i -O split_fa merged.fa
    for fa in split_fa/*.fa; do
      fasta_id=$(grep -m1 '^>' "$fa" | sed 's/^>//' | cut -d' ' -f1)
      clean_fasta_id=$(echo "$fasta_id" | sed 's/[^A-Za-z0-9._-]/_/g')
      clean_id="~{sample_id}_${clean_fasta_id}"
      minimap2 -c --cs -t 4 ~{gene_seq} "$fa" | sort -k6,6 -k8,8n | paftools.js call -f ~{gene_seq} -L50 -l1 - > raw_vcfs/${clean_id}.raw.vcf
    done
  >>>

  output {
    File bed = "~{sample_id}.bed"
    Array[File] raw_vcfs = glob("raw_vcfs/*.raw.vcf")
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          4,
    mem_gb:             16,
    disk_gb:            16,
    boot_disk_gb:       16,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/assembly_eval:0.5"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 select_first([runtime_attr.docker,            default_attr.docker])
  }
}

task normalize_merge_vcfs {
  input {
    Array[File] raw_vcfs
    File gene_seq
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    mkdir -p normed_vcfs
    for vcf in ~{sep=" " raw_vcfs}; do
      base=$(basename "$vcf" .raw.vcf)
      printf "sample\t%s\n" "$base" > sample_map.txt
      bcftools reheader -s sample_map.txt "$vcf" | bcftools norm -f ~{gene_seq} -m -any -Oz -o normed_vcfs/${base}.norm.vcf.gz
      bcftools index -t normed_vcfs/${base}.norm.vcf.gz
    done
    shopt -s nullglob
    vcfs=(normed_vcfs/*.norm.vcf.gz)
    if [ ${#vcfs[@]} -eq 0 ]; then
      echo "ERROR: no normalized VCFs generated" >&2
      exit 1
    elif [ ${#vcfs[@]} -eq 1 ]; then
      echo "INFO: only one VCF found, skipping merge"
      cp "${vcfs[0]}" ~{sample_id}.merged.vcf.gz
      cp "${vcfs[0]}.tbi" ~{sample_id}.merged.vcf.gz.tbi
    else
      bcftools merge "${vcfs[@]}" -Oz -o ~{sample_id}.merged.vcf.gz
      bcftools index -t ~{sample_id}.merged.vcf.gz
    fi
  >>>

  output {
    File vcf = "~{sample_id}.merged.vcf.gz"
    File vcf_idx = "~{sample_id}.merged.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          1,
    mem_gb:             8,
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/binf-basics:0.2"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 select_first([runtime_attr.docker,            default_attr.docker])
  }
}

struct RuntimeAttr {
  Float? mem_gb
  Int? cpu_cores
  Int? disk_gb
  Int? boot_disk_gb
  Int? preemptible_tries
  Int? max_retries
  String? docker
}
