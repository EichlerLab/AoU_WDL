version 1.0


workflow intersect_pav {
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
    }
    meta {
        workflow_description: "Creates callerset for a single sample"
    }

    scatter (vcf_shard_pair in zip(input_vcfs, input_sample_names)) {
        call ProcessPavVcf {
            input:
                vcfIn = vcf_shard_pair.left,
                sample = vcf_shard_pair.right,
        }
        call IntersectWithBackground {
            input:
                bed = ProcessPavVcf.bed,
                sample = vcf_shard_pair.right,
                intersect_py = intersect_script,
                background = background_bed,
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

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -e
    echo -e "#CHROM\tPOS\tEND\tREF\tALT\tGT\tID\tSVTYPE" > ~{sample + "_PAV.bed"}
    bcftools query -f "%CHROM\t%POS\t%END\t%REF\t%ALT\t[%GT]\t%INFO/ID\t%INFO/SVTYPE\n" ~{vcfIn} | grep -v SNV >> ~{sample + "_PAV.bed"}
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