version 1.0

workflow ParalogousVariants {
  input {
    # TSV with header: sample_name  site  hap1  hap2
    # hap1/hap2 columns hold file paths to each haplotype's assembly FASTA.
    File sample_tsv

    # Which "site" values (2nd column of sample_tsv) to run.
    # ["ALL"] (the default) runs every sample regardless of site.
    # ["UW"] runs only site == "UW". ["BI","UW"] runs BI and UW samples.
    Array[String] sites_to_run = ["ALL"]

    File gene_seq
    File fixed_variants

    # Standalone conversion script (kept alongside this WDL); reused as-is
    # by convert_coordinates instead of being generated on every run.
    # It derives the chromosome name/offset directly from gene_seq's own
    # FASTA header (">chrN:start-end"), so no chrom.sizes table is needed.
    File conversion_script = "gene_to_chrom_vcf.py"

    # VEP annotation inputs (see VEP.wdl)
    File gff
    File gff_idx
    File fasta
    File AlphaMissense
    File AlphaMissense_idx
  }

  meta {
    author: "Luyao Ren"
    workflow_description: "Detect paralog variants from assemblies across all samples and haplotypes (optionally restricted to a set of sites), convert gene-based coordinates to chromosome-based coordinates, merge all samples into one cohort VCF, and annotate with VEP."
  }

  # ----------------------------------------------------------------------
  # Step -1 (NEW): filter the sample sheet down to the requested site(s)
  # and pull out the sample_name / hap1 / hap2 columns as parallel arrays.
  # ----------------------------------------------------------------------
  call filter_sample_sheet {
    input:
      sample_tsv = sample_tsv,
      sites_to_run = sites_to_run
  }

  Array[String] sample_names          = filter_sample_sheet.sample_names
  Array[File]   haplotig_fasta_hap1   = filter_sample_sheet.hap1_files
  Array[File]   haplotig_fasta_hap2   = filter_sample_sheet.hap2_files

  Int n_samples = length(sample_names)

  # ----------------------------------------------------------------------
  # Step -0.5 (NEW): parse gene_seq's FASTA header ONCE (not once per
  # sample/hap). Every sample shares the same gene_seq, so the resulting
  # chrom/start/contig_length are identical for all of them -- computing
  # them here and passing the plain values into convert_coordinates avoids
  # re-parsing gene_seq and re-doing the chrom.sizes lookup on every one of
  # the (potentially thousands of) per-sample convert_coordinates calls.
  # ----------------------------------------------------------------------
  call parse_gene_seq {
    input:
      gene_seq = gene_seq
  }

  # ----------------------------------------------------------------------
  # Step 0: build a single FLAT array of (sample, hap, fasta) combinations.
  # WDL does not allow a scatter nested inside another scatter, so instead
  # of scatter(sample){ scatter(hap){...} } we build, for each sample, a
  # 2-element array [hap1_entry, hap2_entry] in one (non-nested) scatter,
  # then flatten() the resulting Array[Array[SampleHap]] into a flat
  # Array[SampleHap] of length 2*n_samples ([s0_hap1, s0_hap2, s1_hap1, ...]).
  # ----------------------------------------------------------------------
  scatter (i in range(n_samples)) {
    SampleHap sh_hap1 = object {
      sample_id: sample_names[i],
      hap: "hap1",
      haplotig_fasta: haplotig_fasta_hap1[i]
    }
    SampleHap sh_hap2 = object {
      sample_id: sample_names[i],
      hap: "hap2",
      haplotig_fasta: haplotig_fasta_hap2[i]
    }
    Array[SampleHap] sample_hap_pair = [sh_hap1, sh_hap2]
  }
  Array[SampleHap] all_sample_haps = flatten(sample_hap_pair)

  # ----------------------------------------------------------------------
  # Step 1: single flat scatter across ALL samples AND haplotypes together
  # (2 * n_samples iterations), instead of one scatter per sample.
  # ----------------------------------------------------------------------
  scatter (sh in all_sample_haps) {
    String sample_w_hap = sh.sample_id + "_" + sh.hap

    call extract_sample_sequence {
      input:
        asm = sh.haplotig_fasta,
        gene_seq = gene_seq,
        sample_id = sample_w_hap
    }
    call best_match {
      input:
        paf = extract_sample_sequence.paf,
        asm = sh.haplotig_fasta,
        fixed_variants = fixed_variants,
        sample_id = sample_w_hap
    }
  }

  # ----------------------------------------------------------------------
  # Step 2: re-group the flat hap1/hap2 results back per sample. Since
  # all_sample_haps was built as [s0_hap1, s0_hap2, s1_hap1, s1_hap2, ...],
  # sample i's two entries sit at indices 2*i and 2*i+1. This is a second,
  # sibling (not nested) scatter over samples.
  # ----------------------------------------------------------------------
  scatter (i in range(n_samples)) {
    Int hap1_idx = 2 * i
    Int hap2_idx = 2 * i + 1
    Array[File] sample_best_match_fa  = [best_match.best_match_fa[hap1_idx], best_match.best_match_fa[hap2_idx]]
    Array[File] sample_best_match_bed = [best_match.best_match_bed[hap1_idx], best_match.best_match_bed[hap2_idx]]

    call call_variants {
      input:
        gene_seq = gene_seq,
        best_match_bed = sample_best_match_bed,
        best_match_fa = sample_best_match_fa,
        sample_id = sample_names[i]
    }

    call normalize_merge_vcfs {
      input:
        gene_seq = gene_seq,
        raw_vcfs = call_variants.raw_vcfs,
        sample_id = sample_names[i]
    }

    # NEW: convert this sample's gene-based (region) coordinates to
    # real chromosome-based coordinates, reusing the chrom/start/length
    # already parsed once by parse_gene_seq above.
    call convert_coordinates {
      input:
        vcf = normalize_merge_vcfs.vcf,
        chrom = parse_gene_seq.chrom,
        start = parse_gene_seq.start,
        contig_length = parse_gene_seq.contig_length,
        sample_id = sample_names[i],
        conversion_script = conversion_script
    }
  }

  # ----------------------------------------------------------------------
  # Step 3 (NEW): merge every sample's chromosome-based VCF into one
  # cohort-level VCF.
  # ----------------------------------------------------------------------
  call merge_all_vcfs {
    input:
      vcfs = convert_coordinates.converted_vcf,
      vcf_idxs = convert_coordinates.converted_vcf_idx,
      cohort_id = "cohort"
  }

  # ----------------------------------------------------------------------
  # Step 4 (NEW): annotate the merged cohort VCF with VEP.
  # ----------------------------------------------------------------------
  call vep {
    input:
      vcf = merge_all_vcfs.merged_vcf,
      gff = gff,
      gff_idx = gff_idx,
      fasta = fasta,
      AlphaMissense = AlphaMissense,
      AlphaMissense_idx = AlphaMissense_idx,
      sample_id = "cohort"
  }

  output {
    Array[File] final_bed             = call_variants.bed
    Array[File] per_sample_vcf        = normalize_merge_vcfs.vcf
    Array[File] per_sample_vcf_idx    = normalize_merge_vcfs.vcf_idx
    Array[File] per_sample_chrom_vcf  = convert_coordinates.converted_vcf
    Array[File] per_sample_chrom_vcf_idx = convert_coordinates.converted_vcf_idx
    File final_merged_vcf     = merge_all_vcfs.merged_vcf
    File final_merged_vcf_idx = merge_all_vcfs.merged_vcf_idx
    File vep_output            = vep.vep_output
  }
}


# ---------------------------------------------------------------------------
# NEW TASK: filter the sample_tsv (sample_name, site, hap1, hap2) down to
# the site(s) requested in sites_to_run, and emit parallel manifests of
# sample_name / hap1 path / hap2 path (one entry per line, matching order)
# for the workflow to read back in as Array[String] / Array[File].
#
# sites_to_run semantics:
#   ["ALL"]        -> keep every sample regardless of its site
#   ["UW"]         -> keep only samples where site == "UW"
#   ["BI","UW"]    -> keep samples where site == "BI" OR site == "UW"
# ---------------------------------------------------------------------------
task filter_sample_sheet {
  input {
    File sample_tsv
    Array[String] sites_to_run

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    sites_csv="~{sep=',' sites_to_run}"

    run_all=0
    IFS=',' read -ra sites_arr <<< "$sites_csv"
    for s in "${sites_arr[@]}"; do
      if [ "$s" == "ALL" ]; then
        run_all=1
      fi
    done

    : > sample_names.txt
    : > hap1_paths.txt
    : > hap2_paths.txt

    awk -F'\t' -v OFS='\t' -v run_all="$run_all" -v sites_csv="$sites_csv" '
      BEGIN { n = split(sites_csv, wanted, ",") }
      NR==1 {
        for (i=1;i<=NF;i++) col[$i]=i
        print $0 > "filtered_samples.tsv"
        next
      }
      {
        site = $(col["site"])
        keep = 0
        if (run_all == 1) {
          keep = 1
        } else {
          for (j=1;j<=n;j++) {
            if (wanted[j] == site) { keep = 1; break }
          }
        }
        if (keep) {
          print $0 >> "filtered_samples.tsv"
          print $(col["sample_name"]) >> "sample_names.txt"
          print $(col["hap1"]) >> "hap1_paths.txt"
          print $(col["hap2"]) >> "hap2_paths.txt"
        }
      }
    ' ~{sample_tsv}

    n_kept=$(wc -l < sample_names.txt)
    echo "INFO: ${n_kept} sample(s) matched sites_to_run=[${sites_csv}]" >&2
    if [ "$n_kept" -eq 0 ]; then
      echo "ERROR: no samples matched sites_to_run=[${sites_csv}]. Check the 'site' column values in the sample sheet." >&2
      exit 1
    fi
  >>>

  output {
    File filtered_tsv = "filtered_samples.tsv"
    Array[String] sample_names = read_lines("sample_names.txt")
    Array[File] hap1_files = read_lines("hap1_paths.txt")
    Array[File] hap2_files = read_lines("hap2_paths.txt")
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          1,
    mem_gb:             2,
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/binf-basics:0.2"
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

task extract_sample_sequence {
  input {
    File asm
    File gene_seq
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    minimap2 -x asm20 -c --secondary=yes -p 0.3 -N 20 --eqx -t 4 -r 500 -K 200M ~{asm} ~{gene_seq} > ~{sample_id}.paf
  >>>

  output {
    File paf = "~{sample_id}.paf"
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          4,
    mem_gb:             32,
    disk_gb:            32,
    boot_disk_gb:       32,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/align-basics:0.2"
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

task best_match {
  input {
    File paf
    File asm
    File fixed_variants
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    rb stat --paf ~{paf} > ~{sample_id}.bed
    bedtools getfasta -fi ~{asm} -bed ~{sample_id}.bed -fo ~{sample_id}.fa
    seqkit locate -f ~{fixed_variants} ~{sample_id}.fa |  awk '
      NR==1 { next }
      {
        split($1, a, /[:\-]/)
        key = a[1] "\t" a[2] "\t" a[3]
        count[key]++
      }
      END {
        for (k in count) {
          if (count[k] >= 5)
            print k
        }
      }
      ' | awk -v sid="~{sample_id}" 'BEGIN{OFS="\t"} {print $0, sid}' > ~{sample_id}.best_match.bed
    bedtools getfasta -fi ~{asm} -bed ~{sample_id}.best_match.bed -fo ~{sample_id}.best_match.fa
  >>>

  output {
    File bed = "~{sample_id}.bed"
    File sample_fa = "~{sample_id}.fa"
    File best_match_bed = "~{sample_id}.best_match.bed"
    File best_match_fa = "~{sample_id}.best_match.fa"
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          1,
    mem_gb:             8,
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/assembly_eval:0.5"
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

task call_variants {
  input {
    Array[File] best_match_fa
    Array[File] best_match_bed
    File gene_seq
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    cat ~{sep=" " best_match_fa} > merged.fa
    cat ~{sep=" " best_match_bed} > ~{sample_id}.bed
    mkdir -p split_fa raw_vcfs
    seqkit split -i -O split_fa merged.fa
    for fa in split_fa/*.fa; do
      fasta_id=$(grep -m1 '^>' "$fa" | sed 's/^>//' | cut -d' ' -f1)
      clean_fasta_id=$(echo "$fasta_id" | sed 's/[^A-Za-z0-9._-]/_/g')
      clean_id="~{sample_id}_${clean_fasta_id}"
      minimap2 -c --cs -t 4 ~{gene_seq} "$fa" | sort -k6,6 -k8,8n | paftools.js call -f ~{gene_seq} -L50 -l1 - > raw_vcfs/${clean_id}.raw.vcf
    done
  >>>

  output {
    File bed = "~{sample_id}.bed"
    Array[File] raw_vcfs = glob("raw_vcfs/*.raw.vcf")
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          4,
    mem_gb:             16,
    disk_gb:            16,
    boot_disk_gb:       16,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/assembly_eval:0.5"
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

task normalize_merge_vcfs {
  input {
    Array[File] raw_vcfs
    File gene_seq
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    mkdir -p normed_vcfs
    for vcf in ~{sep=" " raw_vcfs}; do
      base=$(basename "$vcf" .raw.vcf)
      printf "sample\t%s\n" "$base" > sample_map.txt
      bcftools reheader -s sample_map.txt "$vcf" | bcftools norm -f ~{gene_seq} -m -any -Oz -o normed_vcfs/${base}.norm.vcf.gz
      bcftools index -t normed_vcfs/${base}.norm.vcf.gz
    done
    shopt -s nullglob
    vcfs=(normed_vcfs/*.norm.vcf.gz)
    if [ ${#vcfs[@]} -eq 0 ]; then
      echo "ERROR: no normalized VCFs generated" >&2
      exit 1
    elif [ ${#vcfs[@]} -eq 1 ]; then
      echo "INFO: only one VCF found, skipping merge"
      cp "${vcfs[0]}" ~{sample_id}.merged.vcf.gz
      cp "${vcfs[0]}.tbi" ~{sample_id}.merged.vcf.gz.tbi
    else
      bcftools merge "${vcfs[@]}" -Oz -o ~{sample_id}.merged.vcf.gz
      bcftools index -t ~{sample_id}.merged.vcf.gz
    fi
  >>>

  output {
    File vcf = "~{sample_id}.merged.vcf.gz"
    File vcf_idx = "~{sample_id}.merged.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          1,
    mem_gb:             8,
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/binf-basics:0.2"
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

# ---------------------------------------------------------------------------
# NEW TASK: parse gene_seq's FASTA header ONCE for the whole workflow.
# Header format: ">chrN:start-end[ optional description]" e.g.
#   >chr15:28419321-28443019
# Emits the real chromosome name, the region's start offset (for POS
# conversion), and the TRUE full chromosome length (looked up from the
# built-in T2T-CHM13v2 table below, keyed by the parsed chromosome name --
# no separate chrom.sizes file input needed). These three plain values are
# then reused by every convert_coordinates call instead of each of the
# (potentially thousands of) per-sample calls re-parsing gene_seq itself.
# ---------------------------------------------------------------------------
task parse_gene_seq {
  input {
    File gene_seq

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    # ">chr15:28419321-28443019 optional description" -> "chr15:28419321-28443019"
    seq_id=$(head -n1 ~{gene_seq} | sed 's/^>//' | awk '{print $1}')
    chrom=$(echo "$seq_id" | cut -d':' -f1)
    region=$(echo "$seq_id" | cut -d':' -f2)
    start=$(echo "$region" | cut -d'-' -f1)
    end=$(echo "$region" | cut -d'-' -f2)

    declare -A CHROM_SIZES=(
      [chr1]=248387328   [chr2]=242696752   [chr3]=201105948   [chr4]=193574945
      [chr5]=182045439   [chr6]=172126628   [chr7]=160567428   [chr8]=146259331
      [chr9]=150617247   [chr10]=134758134  [chr11]=135127769  [chr12]=133324548
      [chr13]=113566686  [chr14]=101161492  [chr15]=99753195   [chr16]=96330374
      [chr17]=84276897   [chr18]=80542538   [chr19]=61707364   [chr20]=66210255
      [chr21]=45090682   [chr22]=51324926   [chrX]=154259566   [chrY]=62460029
    )

    if [ -n "${CHROM_SIZES[$chrom]+x}" ]; then
      contig_length="${CHROM_SIZES[$chrom]}"
    else
      echo "Warning: '$chrom' not found in the built-in chrom.sizes table; using the region end ($end) as the contig length instead of the true chromosome length." >&2
      contig_length="$end"
    fi

    echo "$chrom" > chrom.txt
    echo "$start" > start.txt
    echo "$contig_length" > contig_length.txt
  >>>

  output {
    String chrom = read_string("chrom.txt")
    Int start = read_int("start.txt")
    Int contig_length = read_int("contig_length.txt")
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          1,
    mem_gb:             1,
    disk_gb:            5,
    boot_disk_gb:       5,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/binf-basics:0.2"
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

# ---------------------------------------------------------------------------
# NEW TASK: convert this sample's merged, gene/region-based VCF
# (CHROM = "chrN:start-end", POS = 1-based offset within that region) into
# real chromosome-based coordinates (CHROM = "chrN", POS = genomic position).
# Also rewrites the ##contig header line(s) to the true chromosome length.
#
# chrom / start / contig_length are computed ONCE by parse_gene_seq (above)
# and simply passed in here, so this task does no FASTA parsing or
# chrom.sizes lookup itself -- important since it runs once per sample
# and may be scattered across thousands of samples.
#
# The conversion logic lives in the standalone gene_to_chrom_vcf.py file
# (passed in as `conversion_script`) instead of being written out fresh
# inside the command block on every run.
# ---------------------------------------------------------------------------
task convert_coordinates {
  input {
    File vcf
    String chrom
    Int start
    Int contig_length
    String sample_id
    File conversion_script

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    bcftools view ~{vcf} > input.vcf

    python3 ~{conversion_script} input.vcf \
      --chrom ~{chrom} \
      --start ~{start} \
      --length ~{contig_length} \
      -o ~{sample_id}.chrom.vcf

    bgzip -c ~{sample_id}.chrom.vcf > ~{sample_id}.chrom.unsorted.vcf.gz
    bcftools sort -Oz -o ~{sample_id}.chrom.vcf.gz ~{sample_id}.chrom.unsorted.vcf.gz
    bcftools index -t ~{sample_id}.chrom.vcf.gz
  >>>

  output {
    File converted_vcf = "~{sample_id}.chrom.vcf.gz"
    File converted_vcf_idx = "~{sample_id}.chrom.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          1,
    mem_gb:             4,
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/binf-basics:0.2"
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

# ---------------------------------------------------------------------------
# NEW TASK: merge every sample's chromosome-based VCF into one cohort VCF.
# ---------------------------------------------------------------------------
task merge_all_vcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String cohort_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    vcf_arr=(~{sep=" " vcfs})
    idx_arr=(~{sep=" " vcf_idxs})
    # Cromwell may localize each File input into its own directory, so make
    # sure each index sits right next to its VCF before calling bcftools.
    for i in "${!idx_arr[@]}"; do
      ln -sf "${idx_arr[$i]}" "${vcf_arr[$i]}.tbi"
    done

    if [ "${#vcf_arr[@]}" -eq 1 ]; then
      echo "INFO: only one VCF found, skipping merge"
      cp "${vcf_arr[0]}" ~{cohort_id}.merged.vcf.gz
      cp "${vcf_arr[0]}.tbi" ~{cohort_id}.merged.vcf.gz.tbi
    else
      bcftools merge "${vcf_arr[@]}" -Oz -o ~{cohort_id}.merged.vcf.gz
      bcftools index -t ~{cohort_id}.merged.vcf.gz
    fi
  >>>

  output {
    File merged_vcf = "~{cohort_id}.merged.vcf.gz"
    File merged_vcf_idx = "~{cohort_id}.merged.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          2,
    mem_gb:             8,
    disk_gb:            20,
    boot_disk_gb:       10,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/binf-basics:0.2"
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

# ---------------------------------------------------------------------------
# VEP annotation (inlined from VEP.wdl)
# ---------------------------------------------------------------------------
task vep {
  input {
    File vcf
    File gff
    File gff_idx
    File fasta
    File AlphaMissense
    File AlphaMissense_idx
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    cp ~{fasta} ref.fasta
    cp ~{fasta}.fai ref.fasta.fai

    vep -i ~{vcf} -o ~{sample_id}.vep_output.txt \
            --gff ~{gff} \
            --fasta ref.fasta \
            --plugin AlphaMissense,file=~{AlphaMissense}
  >>>

  output {
    File vep_output = "~{sample_id}.vep_output.txt"
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          1,
    mem_gb:             8,
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "ensemblorg/ensembl-vep:release_115.2"
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

struct SampleHap {
  String sample_id
  String hap
  File haplotig_fasta
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
