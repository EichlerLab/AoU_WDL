version 1.0


workflow intersect_pav {
    input {
      File ordered_vcf_shards_list
      File ordered_sample_names
      File background_bed
      File intersect_script
      File windows_file
    }

    meta {
        workflow_description: "Counts Haplotype Presence over "
    }
    parameter_meta {
        ordered_vcf_shards_list: "List of VCFs to process"
        background_bed: "Bed file to intersect the background with"
        intersect_script: "Script to run the intersects"
        ordered_sample_names: "Sample names associated with VCFs"
    }

    Array[File] input_vcfs = read_lines(ordered_vcf_shards_list)
    Array[File] input_sample_names = read_lines(ordered_sample_names)
    
    scatter (vcf_shard_pair in zip(input_vcfs, input_sample_names)) {
        call intersectWindows {
            input:
                vcfIn = vcf_shard_pair.left,
                sample = vcf_shard_pair.right,
                windows_bed = windows_file,
        }
        call countGenotypes {
            input:
                bed = intersectWindows.bed,
                sample = vcf_shard_pair.right,
                process_script = intersect_script,
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
    File windows_bed

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -e
    bedtools intersect -a ~{windows_bed} -b ~{vcfIn} -wa -wb > ~{bed}
  >>>


  output {
    File bed = "~{sample}_PAV.bed"
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