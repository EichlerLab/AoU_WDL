version 1.0

workflow vcf_to_bed {

    input {
        File ordered_sample_names
        File ordered_sniffles_vcf_shards_list
        File sniffles_format_script
        File ordered_pav_vcf_shards_list
        File ordered_pbsv_vcf_shards_list
        File filter_indel_script
        File filter_inv_dup_script
        File combine_pbsv_script
    }
    meta {
        workflow_description: "Converts Sniffles Vcfs to Beds for intersect"
    }
    parameter_meta {
        ordered_sample_names: "Sample names associated with VCFs"
        ordered_sniffles_vcf_shards_list: "List of sniffles VCFs to process"
        sniffles_format_script: "Script to format sniffles vcf for intersect"
    }
    
    Array[File] input_sniffles_vcfs = read_lines(ordered_sniffles_vcf_shards_list)
    Array[File] input_sample_names = read_lines(ordered_sample_names)
    Array[File] input_pbsv_vcfs = read_lines(ordered_pbsv_vcf_shards_list)
    Array[File] input_pav_vcfs = read_lines(ordered_pav_vcf_shards_list)

    scatter (vcf_shard_pair_sniff in zip(input_sniffles_vcfs, input_sample_names)) {
        call ProcessSnifflesVcf {
            input:
                snifflesvcfIn = vcf_shard_pair_sniff.left,
                sample = vcf_shard_pair_sniff.right
        }
        call FiltSnifflesSV {
            input:
                snifflesbedIn = ProcessSnifflesVcf.bed,
                filter_sniffles_sv_py = sniffles_format_script,
                svpoplib_tgz = svpoplib,
                sample = vcf_shard_pair_sniff.right,
        }
        
    }
    scatter (vcf_shard_pair_pbsv in zip(input_pbsv_vcfs, input_sample_names)) {
        call ProcessPBSVVcf {
            input:
                pbsvvcfIn = vcf_shard_pair_pbsv.left,
                sample = vcf_shard_pair_pbsv.right,
        }
        call FiltSVIndels {
            input:
                pbsvbedIn = ProcessPBSVVcf.bed,
                sample = vcf_shard_pair_pbsv.right,
                filter_indel_py = filter_indel_script,
                svpoplib_tgz = svpoplib,
        }
        
        call FiltSVDupInv {
            input:
                pbsvbedIn = ProcessPBSVVcf.bed,
                filter_inv_dup_py = filter_inv_dup_script,
                svpoplib_tgz = svpoplib,
                sample = vcf_shard_pair_pbsv.right,
        }
        call CombineSVs {
            input:
                sv_inv_dup_bedIn = FiltSVDupInv.sv_inv_dup_bed,
                sv_indel_bedIn = FiltSVIndels.sv_indel_bed,
                combine_pbsv_py = combine_pbsv_script,
                sample = vcf_shard_pair_pbsv.right,
        }

        
    }

    scatter (vcf_shard_pair_pav in zip(input_pav_vcfs, input_sample_names)) {
            call ProcessPAVVcf {
                input:
                    pavvcfIn = vcf_shard_pair_pav.left,
                    sample = vcf_shard_pair_pav.right
            }
    }
    
    output {
        Array[File] sample_sniffles_beds = ProcessSnifflesVcf.bed
        Array[File] sample_sniffles_sv = FiltSnifflesSV.sniffles_sv_bed
        Array[File] sample_pbsv_beds = ProcessPBSVVcf.bed
        Array[File] sample_pbsv_indels = FiltSVIndels.sv_indel_bed
        Array[File] sample_pbsv_inv_dup = FiltSVDupInv.sv_inv_dup_bed
        Array[File] sample_pbsv_combine = CombineSVs.pbsv_bed
        Array[File] sample_pav_beds = ProcessPAVVcf.bed
    }
    
}

task ProcessSnifflesVcf {
    input {
        File snifflesvcfIn
        String sample
        
        RuntimeAttr? runtime_attr_override
    }
    
    command <<<
        set -e
        bcftools query -H -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%FILTER\\t%INFO/SVTYPE\\t%INFO/END\\t%INFO/SVLEN[\\t%GT][\\t%DR][\\t%DV][]\\n" ~{snifflesvcfIn} > ~{sample + "_Sniffles.bed"} 
    >>>
    
    
    output {
        File bed = "~{sample}_Sniffles.bed"
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

task FiltSnifflesSV {
    input {
        File snifflesbedIn
        File filter_sniffles_sv_py
        File svpoplib_tgz
        String sample
        
        RuntimeAttr? runtime_attr_override
    }
    command <<<
        set -e
        python ~{filter_sniffles_sv_py} -b ~{snifflesbedIn} -s ~{sample} -o ~{sample + "_sniffles_sv.bed"}
    >>>
    
    
    output {
        File sniffles_sv_bed = "~{sample}_sniffles_sv.bed"
    }
    
    #####################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
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

task ProcessPBSVVcf {
    input {
        File pbsvvcfIn
        String sample
        
        RuntimeAttr? runtime_attr_override
    }
    
    command <<<
        set -e
        bcftools query -H -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%FILTER\\t%INFO/SVTYPE\\t%INFO/END\\t%INFO/SVLEN\\t%INFO/CIPOS\\t%INFO/MATEID\\t%INFO/MATEDIST\\t%INFO/IMPRECISE[\\t%GT][\\t%AD][\\t%DP][\\t%SAC][]\\n" ~{pbsvvcfIn} >> ~{sample + "_PBSV.bed"}
    >>>
    

    
    output {
        File bed = "~{sample}_PBSV.bed"
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

task FiltSVIndels {
    input {
        File pbsvbedIn
        File filter_indel_py
        File svpoplib_tgz
        String sample
        
        RuntimeAttr? runtime_attr_override
    }
    command <<<
        set -e
        python ~{filter_indel_py} -b ~{pbsvbedIn} -s ~{sample} -o ~{sample + "_sv_indels.bed"}
    >>>
    
    output {
        File sv_indel_bed = "~{sample}_sv_indels.bed"
    }
    
    #####################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
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

task FiltSVDupInv {
    input {
        File pbsvbedIn
        File filter_inv_dup_py
        File svpoplib_tgz
        String sample
        
        RuntimeAttr? runtime_attr_override
    }
    command <<<
        set -e
        python ~{filter_inv_dup_py} -b ~{pbsvbedIn} -s ~{sample} -o ~{sample + "_sv_inv_dup.bed"}
    >>>
    
    output {
        File sv_inv_dup_bed = "~{sample}_sv_inv_dup.bed"
    }
    #####################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
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

task CombineSVs {
    input {
        File sv_inv_dup_bedIn
        File sv_indel_bedIn
        File combine_pbsv_py
        String sample
        
        RuntimeAttr? runtime_attr_override
    }
    command <<<
        set -e
        python ~{combine_pbsv_py} -a ~{sv_inv_dup_bedIn} -b ~{sv_indel_bedIn} -o ~{sample + "_PBSV_SV.bed"}
    >>>

    output {
        File pbsv_bed = "~{sample}_PBSV_SV.bed"
    }
    #####################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
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

task ProcessPAVVcf {
    input {
        File pavvcfIn
        File sample

        RuntimeAttr? runtime_attr_override
    }
    command <<<
        set -e
        echo -e "#CHROM\tPOS\tEND\tREF\tALT\tGT\tID\tSVTYPE" > ~{sample + "_PAV.bed"}
        bcftools query -f "%CHROM\t%POS\t%END\t%REF\t%ALT\t[%GT]\t%INFO/ID\t%INFO/SVTYPE\n" ~{pavvcfIn} | grep -v SNV >> ~{sample + "_PAV.bed"}
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

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
