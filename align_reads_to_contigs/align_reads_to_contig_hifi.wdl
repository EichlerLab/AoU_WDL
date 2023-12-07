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
        call alignToAsm as alignToAsmH1 {
            input:
                sample = sample,
                hap = "hap1",
                asmIn = asm_h1,
                readsIn = reads
        }
        call alignToAsm as alignToAsmH2 {
            input:
                sample = sample,
                hap = "hap2",
                asmIn = asm_h2,
                readsIn = reads
        }
    }
    output {
        Array[File] sample_h1_bams = alignToAsmH1.bamOut
        Array[File] sample_h1_idx = alignToAsmH1.baiOut
        Array[File] sample_h2_bams = alignToAsmH2.bamOut
        Array[File] sample_h2_idx = alignToAsmH2.baiOut
    }
}

 
task alignToAsm {
  input {
    String sample
    String hap
    File asmIn
    File readsIn

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    minimap2 -a -t 4 -I 10G -Y -x map-pb --eqx -L --cs ~{asmIn} ~{readsIn} | samtools sort -o ~{sample}-asm_~{hap}.minimap2.bam -;samtools index ~{sample}-asm_~{hap}.minimap2.bam
  >>>


  output {
    File bamOut = "~{sample}-asm_~{hap}.minimap2.bam"
    File baiOut = "~{sample}-asm_~{hap}.minimap2.bam.bai"
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          4,
      mem_gb:             32,
      disk_gb:            30,
      boot_disk_gb:       30,
      preemptible_tries:  2,
      max_retries:        1,
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
