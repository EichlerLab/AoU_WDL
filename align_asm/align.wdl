version 1.0


workflow align_asm {
    input {
      File input_sample_names_list
      File input_asm_h1_list
      File input_asm_h2_list
      File ref
    }
    meta {
        workflow_description: "Creates callerset for a single sample"
    }
    Array[String] input_sample_names = read_lines(input_sample_names_list)
    Array[String] input_asm_h1 = read_lines(input_asm_h1_list)
    Array[String] input_asm_h2 = read_lines(input_asm_h2_list)

    scatter (asm_sample_pair in zip(input_asm_h1, input_sample_names)) {
        call alignToRef as alignToRefH1 {
            input:
                asmIn = asm_sample_pair.left,
                sample = asm_sample_pair.right,
                hap = "h1",
                ref = ref,
                threads = 4,
        }
        call compressAndIndex as compressAndIndexH1 {
            input:
                samIn = alignToRefH1.samOut,
                sample = asm_sample_pair.right,
                hap = "h1",
        }
    }
    scatter (asm_sample_pair in zip(input_asm_h2, input_sample_names)) {
        call alignToRef as alignToRefH2 {
            input:
                asmIn = asm_sample_pair.left,
                sample = asm_sample_pair.right,
                hap = "h2",
                ref = ref,
                threads = 4,
        }
        call compressAndIndex as compressAndIndexH2 {
            input:
                samIn = alignToRefH2.samOut,
                sample = asm_sample_pair.right,
                hap = "h2",
        }
    }
    output {
        Array[File] sample_h1_bams = compressAndIndexH1.bamOut
        Array[File] sample_h1_pafs = compressAndIndexH1.pafOut
        Array[File] sample_h1_idx = compressAndIndexH1.indexOut
        Array[File] sample_h2_bams = compressAndIndexH2.bamOut
        Array[File] sample_h2_pafs = compressAndIndexH2.pafOut
        Array[File] sample_h2_idx = compressAndIndexH2.indexOut
    }
}

 
task alignToRef {
  input {
    File asmIn
    String sample
    File ref
    String hap
    Int threads

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    minimap2 -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 --secondary=no -a -t ~{threads} --eqx -Y ~{ref} ~{asmIn} > ~{sample + "-asm_" + hap + ".minimap2.sam"}
  >>>


  output {
    File samOut = "~{sample}-asm_~{hap}.minimap2.sam"
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          threads,
      mem_gb:             32,
      disk_gb:            11,
      boot_disk_gb:       11,
      preemptible_tries:  1,
      max_retries:        1,
      docker:             "eichlerlab/assembly_eval:0.2"
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


task compressAndIndex {
  input {
    File samIn
    String sample
    String hap

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    samtools view -b ~{samIn} | samtools sort -O bam -o ~{sample + "-asm_" + hap + ".minimap2.bam"} -
    samtools index ~{sample + "-asm_" + hap + ".minimap2.bam"}

    k8 $(which paftools.js) sam2paf -L ~{samIn} > ~{sample + "-asm_" + hap + ".minimap2.paf"}
  >>>


  output {
    File bamOut = "~{sample}-asm_~{hap}.minimap2.bam"
    File pafOut = "~{sample}-asm_~{hap}.minimap2.paf"
    File indexOut = "~{sample}-asm_~{hap}.minimap2.bam.bai"
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             4,
      disk_gb:            11,
      boot_disk_gb:       11,
      preemptible_tries:  1,
      max_retries:        1,
      docker:             "eichlerlab/assembly_eval:0.2"
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
