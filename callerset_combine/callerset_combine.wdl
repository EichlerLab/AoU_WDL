version 1.0


workflow callerset_combine {
    input {
      String sample
      File pav_tab
      File pbsv_tab
      File sniffles_tab
      File pav_bg
      File pav_pbsv
      File pav_sniffles
      File pbsv_bg
      File pbsv_sniffles
      File sniffles_bg
      File callerset_py
    }
    meta {
        workflow_description: "Combines callers Sniffles, PAV, PBSV "
    }

    call callersetCombine {
      input:
        sample = sample,
        pav_tab = pav_tab,
        pbsv_tab = pbsv_tab,
        sniffles_tab = sniffles_tab,
        pav_bg = pav_bg,
        pav_pbsv = pav_pbsv,
        pav_sniffles = pav_sniffles,
        pbsv_bg = pbsv_bg,
        pbsv_sniffles = pbsv_sniffles,
        sniffles_bg = sniffles_bg,
        callerset_py = callerset_py
    }
    output {
        File callerset = callersetCombine.bed
    }
}

 
task callersetCombine {
  input {
    String sample
    File pav_tab
    File pbsv_tab
    File sniffles_tab
    File pav_bg
    File pav_pbsv
    File pav_sniffles
    File pbsv_bg
    File pbsv_sniffles
    File sniffles_bg
    File callerset_py
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -e
    python ~{callerset_py} --pav_tab ~{pav_tab} \
    --pbsv_tab ~{pbsv_tab} --sniffles_tab ~{sniffles_tab} \
    --pav_bg ~{pav_bg} --pav_pbsv ~{pav_pbsv} \
    --pav_sniffles ~{pav_sniffles} --pbsv_bg ~{pbsv_bg} \
    --pbsv_sniffles ~{pbsv_sniffles} --pbsv_sniffles ~{pbsv_sniffles} \
    --sniffles_bg ~{sniffles_bg} --output ~{sample}_callerset.tab
  >>>


  output {
    File bed = "~{sample}_callerset.tab"
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             24,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  1,
      max_retries:        0,
      docker:             "us.gcr.io/sound-memory-266616/auo-utils:0.1"
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