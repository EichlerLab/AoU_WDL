version 1.0

workflow bam_to_contig {
  input {
    File bam_to_contig_bash_script
    File make_contig_bed_script
    File ordered_sample_names_w_hap
    File ordered_bams
    File ordered_bais
    File ordered_haplotig_fastas
    File regions_bed
    Int flank_bp
    File run_faidx_script
    File concat_contigs_script
  }
  meta{
      workflow_description: "Converts bam to original haplotigs in reference genome region."
  }

  Array[String] haps = ["hap1", "hap2"]
  Array[File] input_sample_names_w_hap = read_lines(ordered_sample_names_w_hap)
  Array[File] input_bams = read_lines(ordered_bams)
  Array[File] input_bais = read_lines(ordered_bais)
  Array[File] input_haplotig_fastas = read_lines(ordered_haplotig_fastas)

  scatter (i in range(length(input_sample_names_w_hap))) {
      String sample_w_hap = input_sample_names_w_hap[i]
      File bam_in = input_bams[i]
      File bai_in = input_bais[i]
      File haplotig_fasta = input_haplotig_fastas[i]
      
      call RunBamToContig {
          input:
            bam_to_contig_bash_script = bam_to_contig_bash_script,
            make_contig_bed_script = make_contig_bed_script,
            sample_w_hap = sample_w_hap,
            bam_in = bam_in,
            bai_in = bai_in,
            regions_bed = regions_bed,
            flank_bp = flank_bp
      }
      call RunFaidx {
          input:
            run_faidx_script = run_faidx_script,
            contig_regions_file = RunBamToContig.out,
            sample_w_hap = sample_w_hap,
            haplotig_fasta = haplotig_fasta
      }
  }

  call ConcatContigs {
      input:
        concat_contigs_script = concat_contigs_script,
        contig_files = select_all(RunFaidx.out_fasta)
  }

  output{
      File final_out = ConcatContigs.out
  }
}

task RunBamToContig {
    input {
      File bam_to_contig_bash_script
      File make_contig_bed_script
      String sample_w_hap = sample_w_hap
      File bam_in
      File bai_in
      File regions_bed
      Int flank_bp
      String sample_w_hap

      RuntimeAttr? runtime_attr_override
    }
    command <<<
      sh ~{bam_to_contig_bash_script} ~{make_contig_bed_script} ~{bam_in} ~{regions_bed} ~{flank_bp} > ~{sample_w_hap}.txt
    >>>

    output {
      File out = "~{sample_w_hap}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             8,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  1,
      max_retries:        0,
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
      max_retries:        0,
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
      File concat_contigs_script
      Array[File] contig_files

      RuntimeAttr? runtime_attr_override
    }
    command <<<
        python ~{concat_contigs_script} all_contigs_concat.fa ~{sep=" " contig_files}
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
      max_retries:        0,
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
