version 1.0


workflow svpop_inter {
    input {
      File sample_list
      File svpop_tsv
      File input_ref
      File svpop_vcf_header
      Array[File] svpop_vcfs
      Array[File] svpop_indexes
    }
    meta {
        workflow_description: "Creates sampleset from a set of svpop_vcfs"
    }
    call svpop_make_config {
        input:
            sample_list = sample_list,
    }
    call svpop_make_sampleset {
        input:
            ref = input_ref,
            svpop_config = svpop_make_config.config_file,
            svpop_tsv = svpop_tsv,
            svpop_vcfs = svpop_vcfs,
            svpop_vcfs = svpop_indexes,
            threads = 16,
    }
    call svpop_vcf_final {
        input:
            svpop_inter = svpop_make_sampleset.sampleset,
            svpop_vcf_header = svpop_vcf_header,
            sample_list = sample_list,
    }
    output {
        File svpop_vcf = svpop_vcf_final.vcf
        File svpop_vcf_index = svpop_vcf_final.index
    }
}


task svpop_make_config {
    input {
        File sample_list
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        git clone https://github.com/EichlerLab/AoU_WDL.git
        python AoU_WDL/svpop/scripts/svpop_make_config.py -i ~{sample_list}
    >>>

    output {
        File config_file = "config.json"
    }

    #########################
    RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  2,
      max_retries:        1,
      docker:             "us.gcr.io/broad-dsp-lrma/aou-lr/svpop:0.0.17"
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

task svpop_make_sampleset {
    input {
        Array[File] svpop_vcfs
        Array[File] svpop_indexes
        File svpop_config
        File ref
        File svpop_tsv
        Int threads
        RuntimeAttr? runtime_attr_override
    }
    command <<<
        mkdir -p config
        cp -l ~{svpop_config} config/config.json
        cp -l ~{svpop_tsv} config/samples.tsv
        mkdir -p link_data/ref/
        cp -l ~{sep=" " svpop_vcfs} link_data/
        cp -l ~{sep=" " svpop_indexes} link_data/
        cp -l ~{ref} link_data/ref/ref.fa
        snakemake -s svpop/Snakefile -j ~{threads} -k --restart-times 1 results/variant/sampleset/aou/aou/all/all/bed/sv_insdel.bed.gz
    >>>
    output {
        File sampleset = "results/variant/sampleset/aou/aou/all/all/bed/sv_insdel.bed.gz"
    }
    RuntimeAttr default_attr = object {
      cpu_cores:          threads,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  2,
      max_retries:        1,
      docker:             "us.gcr.io/broad-dsp-lrma/aou-lr/svpop:0.0.17"
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
task svpop_vcf_final {
    input {
        File sample_list        
        File svpop_vcf_header
        File svpop_inter
        RuntimeAttr? runtime_attr_override
    }
    command <<<
        git clone https://github.com/EichlerLab/AoU_WDL.git
        python AoU_WDL/svpop/scripts/svpop_inter_vcf.py -s ~{sample_list} -i ~{svpop_inter} -o svpop.body 
        cat ~{svpop_vcf_header} svpop.body | bcftools sort - | bgzip -c > AoU_svpop-inter.vcf.gz
        tabix -p vcf AoU_svpop-inter.vcf.gz
    >>>
    output {
        File vcf = "AoU_svpop-inter.vcf.gz"
        File index = "AoU_svpop-inter.vcf.gz.tbi"
    }
    RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  2,
      max_retries:        1,
      docker:             "us.gcr.io/broad-dsp-lrma/aou-lr/svpop:0.0.17"
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
