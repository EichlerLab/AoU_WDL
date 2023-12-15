version 1.0

workflow align_reads_to_contig_hifi {
    input {
      File input_ordered_sample_names_list
      File input_ordered_asm_h1_list
      File input_ordered_asm_h2_list
      File input_ordered_reads_list
    }
    meta {
        workflow_description: "Align reads fastq for each sample to hap1 and hap2 haplotigs"
    }
    Array[String] input_sample_names = read_lines(input_ordered_sample_names_list)
    Array[String] input_asm_h1 = read_lines(input_ordered_asm_h1_list)
    Array[String] input_asm_h2 = read_lines(input_ordered_asm_h2_list)
    Array[String] input_reads = read_lines(input_ordered_reads_list)

    scatter (i in range(length(input_sample_names))) {
        String sample = input_sample_names[i]
        File asm_h1 = input_asm_h1[i]
        File asm_h2 = input_asm_h2[i]
        File reads = input_reads[i]
        call alignToAsm {
            input:
                sample = sample,
                asmH1In = asm_h1,
                asmH2In = asm_h2,
                readsIn = reads
        }
    }
    output {
        Array[File] sample_crams = alignToAsm.cramOut
        Array[File] sample_idx = alignToAsm.craiOut
    }
}

 
task alignToAsm {
  input {
    String sample
    File asmH1In
    File asmH2In
    File readsIn

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    zcat ~{asmH1In} ~{asmH2In} | bgzip -c > ~{sample}_both_haps.fa.gz
    samtools faidx ~{sample}_both_haps.fa.gz
    minimap2 -a -t 4 -I 10G -Y -x map-pb --eqx -L --cs ~{sample}_both_haps.fa.gz ~{readsIn} | \
    samtools sort -o ~{sample}-asm.minimap2.cram --reference ~{sample}_both_haps.fa.gz
    samtools index ~{sample}-asm.minimap2.cram
  >>>

  output {
    File cramOut = "~{sample}-asm.minimap2.cram"
    File craiOut = "~{sample}-asm.minimap2.cram.crai"
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          4,
      mem_gb:             70,
      disk_gb:            120,
      boot_disk_gb:       120,
      preemptible_tries:  1,
      max_retries:        0,
      docker:             "us.gcr.io/broad-dsp-lrma/lr-asm:0.1.14"
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