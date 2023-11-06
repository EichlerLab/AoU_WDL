version 1.0

workflow pav_to_contig {
  input {
    File pav_to_contig_bash_script
    File make_pav_contig_bed_script
    File ordered_pav_vcfs
    String ordered_sample_names
    File regions_bed
    Int gene_len
    File make_faidx_inputs_script
    File run_faidx_script
    File ordered_haplotig_fastas
  }
  meta{
      workflow_description: "Converts PAV VCF to original haplotigs in reference genome region, mapped from contigs via TIG_REGION PAV INFO field."
  }

  Array[String] haps = ["hap1", "hap2"]
  Array[File] input_sample_names = read_lines(ordered_sample_names)
  Array[File] input_vcfs = read_lines(ordered_pav_vcfs)
  Map[String, File] input_haplotig_fastas = read_map(ordered_haplotig_fastas)
    
  scatter (sample_vcf_pair in zip(input_sample_names, input_vcfs)) {
    String sample = sample_vcf_pair.left
    call RunPavToContig {
      input:
        pav_to_contig_bash_script = pav_to_contig_bash_script,
        make_pav_contig_bed_script = make_pav_contig_bed_script,
        sample = sample,
        vcfIn = sample_vcf_pair.right,
        regions_bed = regions_bed,
        gene_len = gene_len,
    }
    call MakeFaidxInputs {
      input:
        make_faidx_inputs_script = make_faidx_inputs_script,
        sample = sample,
        in_file = RunPavToContig.out
    }

    scatter(hap in haps) {
      String curr_sample_w_hap = sample + "_" + hap
      File curr_sample_haplotig_fasta = input_haplotig_fastas[curr_sample_w_hap]

      call RunFaidx {
          input:
            run_faidx_script = run_faidx_script,
            contig_regions_file = MakeFaidxInputs.out_contig_regions[hap],
            sample_w_hap = curr_sample_w_hap,
            haplotig_fasta = curr_sample_haplotig_fasta
      }
    }
  }
  call ConcatContigs {
      input:
        contig_files = select_all(flatten(RunFaidx.out_fasta))
  }

  output{
      File final_out = ConcatContigs.out
  }
}

task RunPavToContig {
    input {
      File pav_to_contig_bash_script
      File make_pav_contig_bed_script
      String sample
      File vcfIn
      File regions_bed
      Int gene_len
      
      RuntimeAttr? runtime_attr_override
    }
    command <<<
      sh ~{pav_to_contig_bash_script} ~{vcfIn} ~{regions_bed} ~{make_pav_contig_bed_script} ~{gene_len} > ~{sample}_PAV_contigs.tsv
    >>>

    output {
      File out = "~{sample}_PAV_contigs.tsv"
    }
      
    #########################
    RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  1,
      max_retries:        1,
      docker: "us.gcr.io/broad-dsp-lrma/lr-talon:5.0"
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

task MakeFaidxInputs {
    input {
      File make_faidx_inputs_script
      String sample
      File in_file
      
      RuntimeAttr? runtime_attr_override
    }
    command <<<
      python ~{make_faidx_inputs_script} ~{in_file} ~{sample}_hap1.txt ~{sample}_hap2.txt
    >>>

    output {
      Map[String, File] out_contig_regions = {"hap1": "~{sample}_hap1.txt", "hap2": "~{sample}_hap2.txt"}
    }
      
    #########################
    RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  1,
      max_retries:        1,
      docker: "us.gcr.io/broad-dsp-lrma/lr-talon:5.0"
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

task RunFaidx {
    input {
      File run_faidx_script
      File contig_regions_file
      String sample_w_hap
      File haplotig_fasta
      
      RuntimeAttr? runtime_attr_override
    }
    command <<<
      sh ~{run_faidx_script} ~{contig_regions_file} ~{haplotig_fasta} ~{sample_w_hap}.fa
    >>>

    output {
      File out_fasta = "~{sample_w_hap}.fa"
    }
      
    #########################
    RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  1,
      max_retries:        1,
      docker: "us.gcr.io/broad-dsp-lrma/lr-talon:5.0"
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

task ConcatContigs {
    input {
      Array[File] contig_files
      
      RuntimeAttr? runtime_attr_override
    }
    command <<<
      cat ~{sep=" " contig_files} > all_contigs_concat.fa
    >>>

    output {
      File out = "all_contigs_concat.fa"
    }
      
    #########################
    RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  1,
      max_retries:        1,
      docker: "us.gcr.io/broad-dsp-lrma/lr-talon:5.0"
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
