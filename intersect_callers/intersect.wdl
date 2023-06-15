version 1.0

workflow intersect {
    input {
        File intersect_script
        File ordered_bed_shards_list_caller_1
        File ordered_bed_shards_list_caller_2
        File ordered_sample_names
        String caller_string_1
        String caller_string_2
    }
    meta {
        workflow_description: "Performs svpop interesects"
    }
    
    Array[File] input_beds_1 = read_lines(ordered_bed_shards_list_caller_1)
    Array[File] input_beds_2 = read_lines(ordered_bed_shards_list_caller_2)
    Array[File] input_sample_names = read_lines(ordered_sample_names)
    
    scatter (shard_pair in zip(zip (input_beds_1, input_beds_2),input_sample_names) ) {
        call IntersectWithBackground {
            input:
                bedIn_1 = shard_pair.left.left,
                sample = shard_pair.right,
                intersect_py = intersect_script,
                bedIn_2 = shard_pair.left.right,  
                caller_1 = caller_string_1,
                caller_2 = caller_string_2
        
        }
    }
        
    output {
        Array[File] sample_intersect_beds = IntersectWithBackground.int_bed
    }
    
}

task IntersectWithBackground {
    input {
        File bedIn_1
        File bedIn_2
        String sample
        File intersect_py
        String caller_1
        String caller_2
        
        RuntimeAttr? runtime_attr_override
    }
    
    command <<<
        set -e
        python ~{intersect_py} -a ~{bedIn_1} -b ~{bedIn_2} -o ~{sample +"_"+ caller_1 + "_"+ caller_2 + "_int.bed"} -t 4
    >>>
    
    output {
        File int_bed="~{sample}_~{caller_1}_~{caller_2}_int.bed"
    }
    
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/genial-venture-385019/auo-utils:0.1"
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