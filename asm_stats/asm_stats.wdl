version 1.0


workflow asm_stats {
    input {
      File input_sample_names_list
      File input_paf_h1_list
      File input_paf_h2_list
    }
    meta {
        workflow_description: "Assembly stats"
    }
    Array[String] input_sample_names = read_lines(input_sample_names_list)
    Array[String] input_paf_h1 = read_lines(input_paf_h1_list)
    Array[String] input_paf_h2 = read_lines(input_paf_h2_list)

    scatter (paf_sample_pair in zip(input_paf_h1, input_sample_names)) {
        call asmAlign as asmAlignH1 {
            input:
                pafIn = paf_sample_pair.left,
                sample = paf_sample_pair.right,
                hap = "h1",
        }
    }
    scatter (paf_sample_pair in zip(input_paf_h1, input_sample_names)) {
        call asmAlign as asmAlignH2 {
            input:
                pafIn = paf_sample_pair.left,
                sample = paf_sample_pair.right,
                hap = "h2",
        }
    }
    call combineAsm {
        input:
            hap1_list = asmAlignH1.tabOut,
            hap2_list = asmAlignH2.tabOut,
    }
    output {
        File all_lens = combineAsm.tabOut 
    }
}


task asmAlign {
  input {
    File pafIn
    String sample
    String hap

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    cut -f 1,2 ~{pafIn} | sort -u | awk -vOFS="\t" '{print $1,$2,"~{sample}"}' >  ~{sample + "-asm_" + hap + ".ctg.tab"}  
  >>>


  output {
    File tabOut = "~{sample}-asm_~{hap}.ctg.tab"
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



task combineAsm {
  input {
    Array[File] hap1_list
    Array[File] hap2_list

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    paste <( cat ~{sep=" " hap1_list} ) <( cat ~{sep=" " hap2_list} ) > all_lens.ctg.tab
  >>>


  output {
    File tabOut = "all_lens.ctg.tab"
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

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

