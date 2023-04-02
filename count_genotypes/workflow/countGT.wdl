version 1.0


workflow count_gt {
    input {
      File snv_vcf
      File windows_file
      File process_script
      String sample
    }

    meta {
        workflow_description: "Counts Haplotype Presence over 1kb windows"
    }
#    parameter_meta {
#        ordered_vcf_shards_list: "List of VCFs to process"
#        background_bed: "Bed file to intersect the background with"
#        intersect_script: "Script to run the intersects"
#        ordered_sample_names: "Sample names associated with VCFs"
#    }

    call intersectWindows {
        input:
            vcfIn = snv_vcf,
            sample = sample,
            windows_bed = windows_file,
        }
    call countGenotypes {
        input:
            bed = intersectWindows.bed,
            sample = sample,
            process_script = process_script,
    }

    }
    output {
        Array[File] sample_intersect_beds = IntersectWithBackground.int_bed
        Array[File] sample_pav_beds = ProcessPavVcf.bed
    }
}

 
task bedtoolsIntVCF {
  input {
    File vcfIn
    String sample
    File windows_bed

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -e
    bedtools intersect -a ~{windows_bed} -b ~{vcfIn} -wa -wb > ~{sample}_GT.bed
  >>>


  output {
    File bed = "~{sample}_GT.bed"
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  2,
      max_retries:        1,
      docker:             "us.gcr.io/broad-dsp-lrma/lr-base:0.1.1"
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
 


task IntersectWithBackground {

    input {
        File bed
        File background
        File intersect_py
        String sample

        RuntimeAttr? runtime_attr_override 
    }
    
    command <<<
        set -e
        python ~{intersect_py} -a ~{bed} -b ~{background} -o ~{sample + "_bg_int.bed"} -t ~{threads}
    >>>

    output {
        File int_bed = "~{sample}_bg_int.bed"
    } 



  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          4,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  2,
      max_retries:        1,
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