version 1.0


workflow sfs_plot {
    input {
      File input_all_callerset
      File input_matrix
      File aouWDL
    }
    meta {
        workflow_description: "Plots SFS SV results"
    }

    call allMatrix as allMatrix_one {
        input:
            callerset = input_all_callerset,
            matrix = input_matrix,
            support = "3",
            aouWDL = aouWDL,
    }
    call allMatrix as allMatrix_two {
        input:
            callerset = input_all_callerset,
            matrix = input_matrix,
            support = "2",
            aouWDL = aouWDL,
    }
    call allMatrix as allMatrix_three {
        input:
            callerset = input_all_callerset,
            matrix = input_matrix,
            support = "1",
            aouWDL = aouWDL,
    }


    output {
        File all_matrix_three_plot = allMatrix_three.plot 
        File all_matrix_two_plot = allMatrix_two.plot 
        File all_matrix_one_plot = allMatrix_one.plot 
    }
}

task allMatrix {
  input {
    File matrix
    File callerset
    File aouWDL
    String support
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    tar zxvf ~{aouWDL}
    python AoU_WDL/sfs/sfs_plot.py -m ~{matrix} -c ~{callerset} -s ~{support} -o ~{support + "_class.tab"}
  >>>

  output {
    File plot = "~{support}_class.tab"
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             512,
      disk_gb:            64,
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
