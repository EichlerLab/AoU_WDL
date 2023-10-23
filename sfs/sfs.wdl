version 1.0


workflow sfs_sv {
    input {
      File input_sample_names_list
      File input_paf_h1_list
      File input_paf_h2_list
      File all_calls_bed
      File input_callerset_list
      File aouWDL
    }
    meta {
        workflow_description: "Creates callerset for a single sample"
    }
    Array[String] input_sample_names = read_lines(input_sample_names_list)
    Array[String] input_paf_h1 = read_lines(input_paf_h1_list)
    Array[String] input_paf_h2 = read_lines(input_paf_h2_list)
    Array[String] input_callerset = read_lines(input_callerset_list)

    scatter (paf_sample_pair in zip(input_paf_h1, input_sample_names)) {
        call mergeBed as mergeBedH1 {
            input:
                pafIn = paf_sample_pair.left,
                sample = paf_sample_pair.right,
                hap = "h1",
        }
    }
    scatter (paf_sample_pair in zip(input_paf_h2, input_sample_names)) {
        call mergeBed as mergeBedH2 {
            input:
                pafIn = paf_sample_pair.left,
                sample = paf_sample_pair.right,
                hap = "h2",
        }
    }
    scatter (callerset_sample_pair in zip(input_callerset, input_sample_names)) {
        call callersetUnique {
            input:
                callersetIn = callerset_sample_pair.left,
                sample = callerset_sample_pair.right,
                allCallsBed = all_calls_bed,
                aouWDL = aouWDL,
        }
    }
    output {
        Array[File] sample_h1_bed = mergeBedH1.bedOut
        Array[File] sample_h2_bed = mergeBedH2.bedOut
        Array[File] sample_callerset_bed = callersetUnique.unfoundOut
        Array[File] sample_callerset_gt = callersetUnique.gtOut
    }
}

 
task mergeBed {
  input {
    File pafIn
    String sample
    String hap

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    cut -f 6,8,9 | sort -k1,1 -k2,2n | bedtools merge -i - > ~{sample + "-asm_" + hap + ".callable.bed"} 
  >>>


  output {
    File bedOut = "~{sample}-asm_~{hap}.callable.bed"
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  2,
      max_retries:        1,
      docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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


task callersetUnique {
  input {
    File callersetIn
    File allCallsBed
    File aouWDL
    String sample
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    tar zxvf ~{aouWDL}
    python AoU_WDL/sfs/callserset_unfound_gt.py ~{callersetIn} ~{allCallsBed} ~{sample + "-all_calls_unfound.bed"} ~{sample + ".gt.tsv"} ~{sample} 
  >>>


  output {
    File unfoundOut = "~{sample}-all_calls_unfound.bed"
    File gtOut = "~{sample}.gt.tsv"
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  2,
      max_retries:        1,
      docker:             "us.gcr.io/broad-dsp-lrma/lr-talon:5.0"
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
