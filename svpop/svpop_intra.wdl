version 1.0


workflow sfs_sv {
    input {
      String sample
      File sniffles_vcf
      File pbsv_vcf
      File pav_vcf
      File input_ref
      File svpop_config
      File svpop_tsv
      File svpop_vcf_header
    }
    meta {
        workflow_description: "Creates callerset for a single sample"
    }

    call bcftools_norm_initial {
        input:
            ref = input_ref,
            sniffles_vcf = sniffles_vcf, 
            sample = sample,
    }
    call retrieve_seq {
        input:
            ref = input_ref,
            sniffles_nowm = bcftools_norm_initial.norm_vcf,
            sample = sample,
    }
    call concat_norm_sniffles {
        input:
            ref = input_ref,
            sniffles_norm = bcftools_norm_initial.norm_vcf,
            new_seq = retrieve_seq.body,
            sample = sample,
    }
    call run_svpop {
        input:
            sniffles_vcf = concat_norm_sniffles.vcf,
            sample = sample,
            pbsv_vcf = pbsv_vcf,
            ref = input_ref,
            pav_vcf = pav_vcf,
            svpop_config = svpop_config,
            svpop_tsv = svpop_tsv,
            threads = 16,
    }
    call svpop_vcf_body {
        input:
            sample = sample,
            callerset = run_svpop.callerset,
            pbsv_vcf = pbsv_vcf,
    }
    call svpop_vcf_final {
        input:
            sample = sample,
            svpop_body = svpop_vcf_body.body,
            svpop_vcf_header = svpop_vcf_header,
            ref = input_ref,
    }
    output {
        File svpop_vcf = svpop_vcf_final.vcf 
    }
}


task bcftools_norm_initial {
    input {
        File ref
        File sniffles_vcf
        String sample
    }

    command <<<
        less ~{sniffles_vcf} | bcftools norm -f ~{ref} -c s - > ~{sample + ".sniffles_norm.vcf"}
    >>>

    output {
        File norm_vcf = "~{sample}.sniffles_norm.vcf"
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



task retrieve_seq {
    input { 
        File ref
        File sniffles_norm
        String sample
    }
    command <<<
    git clone https://github.com/EichlerLab/AoU_WDL.git
    python AoU_WDL/svpop/scripts/extract_seq.py -i ~{sniffles_norm} -s ~{sample} -r ~{input_ref} -o ~{sample + ".sniffles_norm.vcf"}
    >>>
    output {
        File body = "~{sample}.vcf.body"
    }
    #########################
    RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  2,
      max_retries:        1,
      docker:             "us.gcr.io/broad-dsp-lrma/aou-lr/svpop:0.0.14"
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

task concat_norm_sniffles {
    input {
        File ref
        File sniffles_norm
        File new_seq
        File sample
    }
    command <<<
        cat <( grep "\#" ~{sniffles_norm} ) ~{new_seq} | bcftools norm - -f ~{input_ref} -c s | awk '$5 !~ /<|>/' | bgzip -c > ~{sample + ".sniffles-filt.vcf.gz"}
    >>>
    }
    output {
        File vcf = "~{sample}.sniffles-filt.vcf.gz"
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

call run_svpop {
    input {
        File sniffles_vcf
        String sample
        File pbsv_vcf
        File pav_vcf
        File svpop_config
        File ref
        File svpop_tsv
        Int threads
    }
    command <<<
        mkdir -p config
        cp -l ~{svpop_config} config/config.json
        cp -l ~{svpop_tsv} config/samples.tsv
        mkdir -p link_data/ref/
        cp -l ~{sniffles_vcf} ~{"link_data/" + sample + ".sniffles.vcf.gz"}
        cp -l ~{pav_vcf} ~{"link_data/" + sample + ".pav.vcf.gz"}
        cp -l ~{pbsv_vcf} ~{"link_data/" + sample + ".pbsv.vcf.gz"}
        cp -l ~{ref} link_data/ref/ref.fa
        find link_data -type f | xargs -i tabix -p vcf {}
        snakemake -s svpop/Snakefile -j ~{threads} -k --restart-times 1 ~{"results/variant/callerset/aou/" + sample + "/all/all/bed/sv_insdel.bed.gz"}
    >>>
    output {
        File callerset = "results/variant/callerset/aou/~{sample}/all/all/bed/sv_insdel.bed.gz"
    }
    RuntimeAttr default_attr = object {
      cpu_cores:          threads,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  2,
      max_retries:        1,
      docker:             "us.gcr.io/broad-dsp-lrma/aou-lr/svpop:0.0.14"
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
task svpop_vcf_body {
    input {
        String sample
        File callerset
        File pbsv_vcf
    }
    command <<<
        git clone https://github.com/EichlerLab/AoU_WDL.git
        python AoU_WDL/svpop/scripts/svpop_vcf.py -i ~{callerset} -p ~{pbsv_vcf} -s ~{sample} -o ~{sample + "_svpop.body"}
    >>>
    output {
        File body = "~{sample}_svpop.body" 
    }
    RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  2,
      max_retries:        1,
      docker:             "us.gcr.io/broad-dsp-lrma/aou-lr/svpop:0.0.14"
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
        String sample
        File svpop_body
        File svpop_vcf_header
        File ref
    }
    command <<<
        cat ~{svpop_vcf_header} ~{svpop_body} | bcftools sort - | bcftools norm -f ~{ref} -c s - | bgzip -c > ~{sample + "_svpop.vcf.gz"}
        tabix -p vcf ~{sample + "_svpop.vcf.gz"}
    >>>
    output {
        File vcf = "~{sample}_svpop.vcf.gz"
        File index = "~{sample}_svpop.vcf.gz.tbi"
    }
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
