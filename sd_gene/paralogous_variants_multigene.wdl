version 1.0

workflow ParalogousVariantsMultiGene {
  input {
    # TSV with header: sample_name  site  hap1  hap2
    # hap1/hap2 columns hold file paths to each haplotype's assembly FASTA.
    File sample_tsv

    # Which "site" values (2nd column of sample_tsv) to run.
    # ["ALL"] (the default) runs every sample regardless of site.
    # ["UW"] runs only site == "UW". ["BI","UW"] runs BI and UW samples.
    Array[String] sites_to_run

    # NEW: TSV with header: gene  gene_seq  kmer_fa
    #   gene_seq = <GENE>.chm13.fa  (one record, header ">chrN:start-end", 1-based)
    #   kmer_fa  = <GENE>.k31.fa    (paralog-specific 31-mers, header carries
    #                                gene= pos= own= other=; SNP at the centre)
    # All genes are merged into one panel and run together.
    File gene_tsv

    # Standalone multi-region conversion script (gene_to_chrom_vcf_multi.py,
    # kept alongside this WDL). It converts every region contig in a VCF to
    # chromosome coordinates using the gene_map.tsv built by prepare_gene_panel.
    File conversion_script

    # NEW: standalone script (vep_lollipop.py, kept alongside this WDL) that
    # combines per-variant carrier counts with the VEP output into per-gene
    # .lollipop.txt / .detailed.lollipop.txt (as rule lollipop_df).
    File lollipop_script

    # NEW: k-mer assignment thresholds (see best_match).
    # A candidate assembly region is assigned to a gene when
    #   informative >= required  AND  own / informative >= min_own_frac
    # where informative = diagnostic sites showing either the own or the other
    # allele, and
    #   required = min(total_sites, max(min_sites_floor, ceil(min_site_cov * total_sites)))
    # so the requirement scales with how many diagnostic k-mers the gene has.
    Float min_own_frac = 0.9
    Float min_site_cov = 0.5
    Int   min_sites_floor = 3
    # Genes with fewer diagnostic sites than this are labelled confidence=low.
    Int   low_conf_sites = 10

    # VEP annotation inputs (see VEP.wdl)
    File gff
    File gff_idx
    File fasta
  }

  meta {
    author: "Luyao Ren"
    workflow_description: "Multi-gene version: merge several SD genes' CHM13 sequences into one panel, extract candidate copies from every sample/haplotype assembly with a single minimap2 run, assign each copy to its gene with paralog-specific own/other k-mers, call variants against the matching gene, convert to chromosome coordinates, merge all samples into one cohort VCF, and annotate with VEP."
  }

  # ----------------------------------------------------------------------
  # Step -1: filter the sample sheet down to the requested site(s)
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
  # Step -0.5 (NEW): read the gene manifest and build the gene panel ONCE:
  #   panel.fa         all genes' CHM13 sequences (headers kept as-is)
  #   panel.kmers.fa   own + other 31-mers for every diagnostic site
  #   gene_map.tsv     gene -> region_id / chrom / start / contig_length / n_sites
  # gene_map.tsv replaces the old parse_gene_seq task.
  # ----------------------------------------------------------------------
  call read_gene_manifest {
    input:
      gene_tsv = gene_tsv
  }

  call prepare_gene_panel {
    input:
      gene_names = read_gene_manifest.gene_names,
      gene_seqs = read_gene_manifest.gene_seqs,
      kmer_fas = read_gene_manifest.kmer_fas
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
  # (2 * n_samples iterations). Each assembly is aligned against the whole
  # gene panel in ONE minimap2 run, then every candidate region is assigned
  # to a gene (or dropped) using the own/other k-mers.
  # ----------------------------------------------------------------------
  scatter (sh in all_sample_haps) {
    String sample_w_hap = sh.sample_id + "_" + sh.hap
    String cn_sample = sh.sample_id
    String cn_hap = sh.hap

    call extract_sample_sequence {
      input:
        asm = sh.haplotig_fasta,
        gene_seq = prepare_gene_panel.panel_fa,
        sample_id = sample_w_hap
    }
    call best_match {
      input:
        paf = extract_sample_sequence.paf,
        asm = sh.haplotig_fasta,
        panel_kmers = prepare_gene_panel.panel_kmers,
        gene_map = prepare_gene_panel.gene_map,
        min_own_frac = min_own_frac,
        min_site_cov = min_site_cov,
        min_sites_floor = min_sites_floor,
        low_conf_sites = low_conf_sites,
        sample_id = sample_w_hap
    }
  }

  # ----------------------------------------------------------------------
  # Step 1.5 (NEW): per-haplotype copy number table. CN of a gene on a
  # haplotype = number of regions assigned to it (status == PASS in that
  # haplotype's kmer_stats.tsv); genes with no assigned region get 0.
  # One row per sample/haplotype, one column per gene (gene_tsv order).
  # ----------------------------------------------------------------------
  call summarize_cn {
    input:
      kmer_stats = best_match.kmer_stats,
      sample_ids = cn_sample,
      haps = cn_hap,
      gene_map = prepare_gene_panel.gene_map,
      cohort_id = "cohort"
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
        panel_fa = prepare_gene_panel.panel_fa,
        best_match_bed = sample_best_match_bed,
        best_match_fa = sample_best_match_fa,
        sample_id = sample_names[i]
    }

    call normalize_merge_vcfs {
      input:
        panel_fa = prepare_gene_panel.panel_fa,
        raw_vcfs = call_variants.raw_vcfs,
        sample_id = sample_names[i]
    }

    # Convert every gene region contig in this sample's VCF to real
    # chromosome-based coordinates, using gene_map.tsv built once above.
    call convert_coordinates {
      input:
        vcf = normalize_merge_vcfs.vcf,
        gene_map = prepare_gene_panel.gene_map,
        sample_id = sample_names[i],
        conversion_script = conversion_script
    }
  }

  # ----------------------------------------------------------------------
  # Step 3: merge every sample's chromosome-based VCF into one
  # cohort-level VCF.
  # ----------------------------------------------------------------------
  call merge_all_vcfs {
    input:
      vcfs = convert_coordinates.converted_vcf,
      vcf_idxs = convert_coordinates.converted_vcf_idx,
      cohort_id = "cohort"
  }

  # ----------------------------------------------------------------------
  # Step 3.5 (NEW): make the VEP input joinable back to the genotypes:
  # split multiallelic records (bcftools merge re-joins them) and give every
  # record a unique ID (CHROM_POS_REF_ALT), which VEP reports as
  # #Uploaded_variation. final_merged_vcf itself is left as it was.
  # ----------------------------------------------------------------------
  call prepare_vep_input {
    input:
      vcf = merge_all_vcfs.merged_vcf,
      cohort_id = "cohort"
  }

  # ----------------------------------------------------------------------
  # Step 4: annotate the merged cohort VCF with VEP.
  # ----------------------------------------------------------------------
  call vep {
    input:
      vcf = prepare_vep_input.vep_input_vcf,
      gff = gff,
      gff_idx = gff_idx,
      fasta = fasta,
      sample_id = "cohort"
  }

  # ----------------------------------------------------------------------
  # Step 5 (NEW): combine variant carrier counts with the VEP output into
  # one lollipop table + one detailed lollipop table for all genes (as rule
  # lollipop_df in SD_genes/paralog_variants_merged, plus a gene column).
  # ----------------------------------------------------------------------
  call lollipop_df {
    input:
      vcf = prepare_vep_input.vep_input_vcf,
      vep_output = vep.vep_output,
      gene_map = prepare_gene_panel.gene_map,
      lollipop_script = lollipop_script,
      cohort_id = "cohort"
  }

  output {
    File gene_map                     = prepare_gene_panel.gene_map
    File panel_fa                     = prepare_gene_panel.panel_fa
    File panel_kmers                  = prepare_gene_panel.panel_kmers
    Array[File] kmer_stats            = best_match.kmer_stats
    File gene_cn_table                = summarize_cn.cn_table
    Array[File] final_bed             = call_variants.bed
    Array[File] per_sample_vcf        = normalize_merge_vcfs.vcf
    Array[File] per_sample_vcf_idx    = normalize_merge_vcfs.vcf_idx
    Array[File] per_sample_chrom_vcf  = convert_coordinates.converted_vcf
    Array[File] per_sample_chrom_vcf_idx = convert_coordinates.converted_vcf_idx
    File final_merged_vcf     = merge_all_vcfs.merged_vcf
    File final_merged_vcf_idx = merge_all_vcfs.merged_vcf_idx
    File vep_input_vcf         = prepare_vep_input.vep_input_vcf
    File vep_output            = vep.vep_output
    File lollipop              = lollipop_df.lollipop
    File detailed_lollipop     = lollipop_df.detailed_lollipop
  }
}


# ---------------------------------------------------------------------------
# Filter the sample_tsv (sample_name, site, hap1, hap2) down to the site(s)
# requested in sites_to_run, and emit parallel manifests of
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

# ---------------------------------------------------------------------------
# NEW TASK: read the gene manifest (gene, gene_seq, kmer_fa) into parallel
# arrays. Done in its own task (like filter_sample_sheet) so the file paths
# are handed back to the workflow as Array[File] and get localized properly
# by the next task.
# ---------------------------------------------------------------------------
task read_gene_manifest {
  input {
    File gene_tsv

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    : > gene_names.txt
    : > gene_seq_paths.txt
    : > kmer_fa_paths.txt

    awk -F'\t' '
      { sub(/\r$/, "") }
      NR==1 {
        for (i=1;i<=NF;i++) col[$i]=i
        if (!("gene" in col) || !("gene_seq" in col) || !("kmer_fa" in col)) {
          print "ERROR: gene_tsv header must contain the columns gene, gene_seq, kmer_fa" > "/dev/stderr"
          bad = 1
          exit 1
        }
        next
      }
      /^[[:space:]]*$/ { next }
      {
        print $(col["gene"])     >> "gene_names.txt"
        print $(col["gene_seq"]) >> "gene_seq_paths.txt"
        print $(col["kmer_fa"])  >> "kmer_fa_paths.txt"
      }
      END { if (bad) exit 1 }
    ' ~{gene_tsv}

    n_genes=$(wc -l < gene_names.txt)
    echo "INFO: ${n_genes} gene(s) in gene_tsv" >&2
    if [ "$n_genes" -eq 0 ]; then
      echo "ERROR: gene_tsv has no gene rows" >&2
      exit 1
    fi
  >>>

  output {
    Array[String] gene_names = read_lines("gene_names.txt")
    Array[File] gene_seqs = read_lines("gene_seq_paths.txt")
    Array[File] kmer_fas = read_lines("kmer_fa_paths.txt")
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
# NEW TASK: build the gene panel ONCE for the whole workflow.
#
#   panel.fa        every gene_seq concatenated; headers stay ">chrN:start-end"
#                   (the region_id), so minimap2/paftools/bcftools see one
#                   contig per gene.
#   panel.kmers.fa  for every diagnostic site in every <GENE>.k31.fa, two
#                   31-mers named "<gene>|<chrom>:<pos>|own" and
#                   "<gene>|<chrom>:<pos>|other". The other-kmer is the
#                   own-kmer with the SNP base (at offset pos - kmer_start - 1)
#                   swapped to the header's other= allele.
#   gene_map.tsv    gene  region_id  chrom  start  end  contig_length  n_sites
#                   (replaces parse_gene_seq; contig_length comes from the
#                   built-in T2T-CHM13v2 table below).
#
# Sanity checks: one record per gene_seq, unique gene names / region_ids,
# kmer gene= must equal the manifest gene, each site must lie inside the
# gene region, and the base at the SNP offset must equal own=.
# ---------------------------------------------------------------------------
task prepare_gene_panel {
  input {
    Array[String] gene_names
    Array[File] gene_seqs
    Array[File] kmer_fas

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    genes=(~{sep=" " gene_names})
    seqs=(~{sep=" " gene_seqs})
    kmers=(~{sep=" " kmer_fas})

    if [ "${#genes[@]}" -ne "${#seqs[@]}" ] || [ "${#genes[@]}" -ne "${#kmers[@]}" ]; then
      echo "ERROR: gene / gene_seq / kmer_fa columns have different lengths" >&2
      exit 1
    fi
    dup_genes=$(printf "%s\n" "${genes[@]}" | sort | uniq -d)
    if [ -n "$dup_genes" ]; then
      echo "ERROR: duplicated gene name(s) in gene_tsv: ${dup_genes}" >&2
      exit 1
    fi

    # T2T-CHM13v2 chromosome lengths
    chrom_sizes="chr1:248387328 chr2:242696752 chr3:201105948 chr4:193574945 chr5:182045439 chr6:172126628 chr7:160567428 chr8:146259331 chr9:150617247 chr10:134758134 chr11:135127769 chr12:133324548 chr13:113566686 chr14:101161492 chr15:99753195 chr16:96330374 chr17:84276897 chr18:80542538 chr19:61707364 chr20:66210255 chr21:45090682 chr22:51324926 chrX:154259566 chrY:62460029"

    : > panel.fa
    : > panel.kmers.fa
    printf "gene\tregion_id\tchrom\tstart\tend\tcontig_length\tn_sites\n" > gene_map.tsv

    for i in "${!genes[@]}"; do
      gene="${genes[$i]}"
      seq_fa="${seqs[$i]}"
      kmer_fa="${kmers[$i]}"

      if ! [[ "$gene" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "ERROR: gene name '${gene}' may only contain letters, digits, '.', '_' or '-'" >&2
        exit 1
      fi

      n_rec=$(grep -c '^>' "$seq_fa" || true)
      if [ "$n_rec" -ne 1 ]; then
        echo "ERROR: ${gene}: gene_seq must contain exactly one FASTA record (found ${n_rec})" >&2
        exit 1
      fi

      # ">chr15:28419321-28443019 optional description" -> "chr15:28419321-28443019"
      region_id=$(grep -m1 '^>' "$seq_fa" | sed 's/^>//' | tr -d '\r' | awk '{print $1}')
      chrom="${region_id%%:*}"
      range="${region_id#*:}"
      start="${range%%-*}"
      end="${range##*-}"
      if ! [[ "$start" =~ ^[0-9]+$ && "$end" =~ ^[0-9]+$ ]]; then
        echo "ERROR: ${gene}: gene_seq header '${region_id}' is not of the form chrN:start-end" >&2
        exit 1
      fi

      contig_length=$(awk -v c="$chrom" -v s="$chrom_sizes" 'BEGIN { n = split(s, a, " "); for (i=1;i<=n;i++) { split(a[i], b, ":"); if (b[1] == c) print b[2] } }')
      if [ -z "$contig_length" ]; then
        echo "Warning: ${gene}: '${chrom}' not found in the built-in chrom.sizes table; using the region end (${end}) as the contig length instead of the true chromosome length." >&2
        contig_length="$end"
      fi

      if [ "$seq_fa" == "$kmer_fa" ] || [ "$(basename "$seq_fa")" == "$(basename "$kmer_fa")" ]; then
        echo "ERROR: ${gene}: gene_seq and kmer_fa point to the same file (${kmer_fa}); check the kmer_fa column of gene_tsv" >&2
        exit 1
      fi

      awk '{ sub(/\r$/, ""); print }' "$seq_fa" >> panel.fa

      # own + other k-mers for every diagnostic site of this gene
      awk -v gene="$gene" -v chrom="$chrom" -v rs="$start" -v re="$end" -v nsfile="n_sites.txt" -v ref_fa="$seq_fa" '
        BEGIN {
          # this gene'"'"'s own CHM13 sequence, to drop sites whose other-kmer
          # also occurs in it (internal repeats): such a site would read as
          # "other" in every true copy and is not diagnostic
          while ((getline line < ref_fa) > 0) {
            sub(/\r$/, "", line)
            if (line !~ /^>/) ref = ref toupper(line)
          }
          close(ref_fa)
          comp["A"] = "T"; comp["C"] = "G"; comp["G"] = "C"; comp["T"] = "A"
        }
        function revcomp(s,   i, r, b) {
          r = ""
          for (i = length(s); i >= 1; i--) {
            b = substr(s, i, 1)
            r = r ((b in comp) ? comp[b] : "N")
          }
          return r
        }
        function fail(msg) {
          print "ERROR: " gene ": " msg " (" hdr ")" > "/dev/stderr"
          failed = 1
          exit 1
        }
        function emit(   f, n, i, kv, id, kc, ks, pos, own, oth, kg, kmer, off, other_kmer, site) {
          if (hdr == "") return
          n = split(hdr, f, /[ \t]+/)
          id = substr(f[1], 2)
          kc = id; sub(/:.*/, "", kc)
          ks = id; sub(/^[^:]*:/, "", ks); sub(/-.*/, "", ks)
          pos = ""; own = ""; oth = ""; kg = ""
          for (i = 2; i <= n; i++) {
            split(f[i], kv, "=")
            if (kv[1] == "gene") kg = kv[2]
            else if (kv[1] == "pos") pos = kv[2]
            else if (kv[1] == "own") own = toupper(kv[2])
            else if (kv[1] == "other") oth = toupper(kv[2])
          }
          if (kg == "") fail("k-mer header has no gene= field; " FILENAME " does not look like a <GENE>.k31.fa, check the kmer_fa column of gene_tsv")
          if (kg != gene) fail("k-mer gene=" kg " does not match the manifest gene; check the kmer_fa column of gene_tsv (" FILENAME ")")
          if (pos !~ /^[0-9]+$/ || ks !~ /^[0-9]+$/) fail("cannot parse k-mer start / pos=")
          if (kc != chrom || pos + 0 < rs + 0 || pos + 0 > re + 0) fail("site lies outside the gene region " chrom ":" rs "-" re)
          kmer = toupper(seq)
          # k-mer header start is 0-based, pos is 1-based -> 0-based offset of the SNP
          off = pos - ks - 1
          if (off < 0 || off >= length(kmer) || substr(kmer, off + 1, 1) != own) fail("base at the SNP offset does not equal own=" own)
          if (oth !~ /^[ACGT]$/ || oth == own) {
            print "Warning: " gene ": skipping site with other=" oth " (" hdr ")" > "/dev/stderr"
            return
          }
          site = chrom ":" pos
          if (site in seen) {
            print "Warning: " gene ": duplicated site " site ", keeping the first k-mer" > "/dev/stderr"
            return
          }
          seen[site] = 1
          other_kmer = substr(kmer, 1, off) oth substr(kmer, off + 2)
          if (index(ref, other_kmer) > 0 || index(ref, revcomp(other_kmer)) > 0) {
            print "Warning: " gene ": dropping site " site ", its other-kmer also occurs in the gene'"'"'s own CHM13 sequence" > "/dev/stderr"
            return
          }
          print ">" gene "|" site "|own"
          print kmer
          print ">" gene "|" site "|other"
          print other_kmer
          nsites++
        }
        { sub(/\r$/, "") }
        /^>/ { emit(); hdr = $0; seq = ""; next }
        { seq = seq $0 }
        END {
          if (failed) exit 1
          emit()
          if (failed) exit 1
          print nsites + 0 > nsfile
        }
      ' "$kmer_fa" >> panel.kmers.fa

      n_sites=$(cat n_sites.txt)
      if [ "$n_sites" -eq 0 ]; then
        echo "Warning: ${gene}: no usable diagnostic k-mers; this gene can never be assigned and will be skipped" >&2
      fi
      echo "INFO: ${gene}: region ${region_id}, ${n_sites} diagnostic site(s)" >&2

      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$gene" "$region_id" "$chrom" "$start" "$end" "$contig_length" "$n_sites" >> gene_map.tsv
    done

    dup_regions=$(tail -n +2 gene_map.tsv | cut -f2 | sort | uniq -d)
    if [ -n "$dup_regions" ]; then
      echo "ERROR: the same gene_seq region is used by more than one gene: ${dup_regions}" >&2
      exit 1
    fi
  >>>

  output {
    File panel_fa = "panel.fa"
    File panel_kmers = "panel.kmers.fa"
    File gene_map = "gene_map.tsv"
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

# ---------------------------------------------------------------------------
# CHANGED: assign every candidate assembly region to a gene with own/other
# k-mers instead of a fixed ">= 5 hits" cut-off.
#
#   1. rb stat -> candidate regions (first 3 columns = assembly coordinates),
#      de-duplicated by coordinates (several gene queries / secondary
#      alignments can hit the same assembly interval).
#   2. seqkit locate (both strands, case-insensitive: assemblies and CHM13
#      may be soft-masked) with panel.kmers.fa.
#   3. per (region, gene), each diagnostic site is counted ONCE as
#        own    only the own-kmer was found
#        other  only the other-kmer was found
#        mixed  both were found (e.g. an internal duplication) -> counts as
#               informative but not as own, i.e. conservatively
#      informative = own + other + mixed, own_frac = own / informative
#   4. required = min(total_sites, max(min_sites_floor, ceil(min_site_cov * total_sites)))
#      PASS when informative >= required and own_frac >= min_own_frac.
#      If a region passes for more than one gene, keep the gene with the
#      highest own_frac (then most own sites) and warn.
#   5. confidence = low when the gene has < low_conf_sites diagnostic sites.
#
# A sample/hap with no region for some gene simply has no row for it
# (skipped); kmer_stats.tsv lists every (region, gene) pair that had any hit.
#
# best_match.bed columns:
#   contig start end sample_hap gene region_id own other mixed informative own_frac confidence
# ---------------------------------------------------------------------------
task best_match {
  input {
    File paf
    File asm
    File panel_kmers
    File gene_map
    Float min_own_frac
    Float min_site_cov
    Int min_sites_floor
    Int low_conf_sites
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    rb stat --paf ~{paf} > ~{sample_id}.bed

    awk -v OFS='\t' '!/^#/ && NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ { print $1, $2, $3 }' ~{sample_id}.bed \
      | sort -k1,1 -k2,2n -k3,3n -u > ~{sample_id}.candidates.bed

    : > ~{sample_id}.fa
    : > ~{sample_id}.kmer_hits.tsv
    if [ -s ~{sample_id}.candidates.bed ]; then
      bedtools getfasta -fi ~{asm} -bed ~{sample_id}.candidates.bed -fo ~{sample_id}.fa
      seqkit locate -i -f ~{panel_kmers} ~{sample_id}.fa > ~{sample_id}.kmer_hits.tsv
    else
      echo "Warning: ~{sample_id}: no candidate regions in the PAF" >&2
    fi

    awk -F'\t' -v OFS='\t' -v sid="~{sample_id}" \
        -v min_frac="~{min_own_frac}" -v min_cov="~{min_site_cov}" \
        -v floor_sites="~{min_sites_floor}" -v low_sites="~{low_conf_sites}" \
        -v gmap="~{gene_map}" -v cands="~{sample_id}.candidates.bed" \
        -v stats_out="~{sample_id}.kmer_stats.tsv" -v bed_out="~{sample_id}.best_match.unsorted.bed" '
      FILENAME == gmap {
        if (FNR > 1) { total[$1] = $7; region_of[$1] = $2 }
        next
      }
      FILENAME == cands {
        key = $1 ":" $2 "-" $3
        c_chr[key] = $1; c_s[key] = $2; c_e[key] = $3
        next
      }
      # seqkit locate output: seqID patternName pattern strand start end matched
      FNR == 1 { next }
      {
        np = split($2, p, "|")
        if (np != 3) next
        k = $1 SUBSEP p[1] SUBSEP p[2]
        if (!(k in state)) state[k] = p[3]
        else if (state[k] != p[3]) state[k] = "mixed"
        pair[$1 SUBSEP p[1]] = 1
      }
      END {
        eps = 1e-9
        for (k in state) {
          split(k, a, SUBSEP)
          rg = a[1] SUBSEP a[2]
          if (state[k] == "own") n_own[rg]++
          else if (state[k] == "other") n_oth[rg]++
          else n_mix[rg]++
        }

        n = 0
        for (rg in pair) {
          split(rg, a, SUBSEP)
          seqid = a[1]; g = a[2]
          if (!(seqid in c_chr)) {
            print "ERROR: k-mer hit on unknown sequence " seqid > "/dev/stderr"
            exit 1
          }
          tot = total[g] + 0
          o = n_own[rg] + 0; t = n_oth[rg] + 0; m = n_mix[rg] + 0
          inf = o + t + m

          x = min_cov * tot
          req = int(x - eps); if (req < x - eps) req++
          if (req < floor_sites) req = floor_sites
          if (req > tot) req = tot

          frac = (inf > 0) ? o / inf : 0
          conf = (tot < low_sites) ? "low" : "high"

          if (inf == 0 || inf < req) status = "FAIL_sites"
          else if (frac < min_frac - eps) status = "FAIL_own_frac"
          else status = "PASS"

          n++
          r_seq[n] = seqid; r_gene[n] = g; r_tot[n] = tot; r_own[n] = o; r_oth[n] = t
          r_mix[n] = m; r_inf[n] = inf; r_req[n] = req; r_frac[n] = frac; r_conf[n] = conf; r_stat[n] = status

          if (status == "PASS") {
            npass[seqid]++
            if (!(seqid in best) || frac > best_frac[seqid] + eps || (frac >= best_frac[seqid] - eps && o > best_own[seqid])) {
              best[seqid] = n; best_frac[seqid] = frac; best_own[seqid] = o
            }
          }
        }

        print "contig", "start", "end", "sample_hap", "gene", "total_sites", "own", "other", "mixed", "informative", "required", "own_frac", "confidence", "status" > stats_out
        for (i = 1; i <= n; i++) {
          s = r_seq[i]
          status = r_stat[i]
          if (status == "PASS" && best[s] != i) status = "PASS_not_best"
          print c_chr[s], c_s[s], c_e[s], sid, r_gene[i], r_tot[i], r_own[i], r_oth[i], r_mix[i], r_inf[i], r_req[i], sprintf("%.4f", r_frac[i]), r_conf[i], status > stats_out
          if (status == "PASS") {
            print c_chr[s], c_s[s], c_e[s], sid, r_gene[i], region_of[r_gene[i]], r_own[i], r_oth[i], r_mix[i], r_inf[i], sprintf("%.4f", r_frac[i]), r_conf[i] > bed_out
          }
        }
        for (s in npass) {
          if (npass[s] > 1) print "Warning: " sid ": region " s " passes for " npass[s] " genes; kept " r_gene[best[s]] > "/dev/stderr"
        }
        close(stats_out); close(bed_out)
      }
    ' ~{gene_map} ~{sample_id}.candidates.bed ~{sample_id}.kmer_hits.tsv

    touch ~{sample_id}.best_match.unsorted.bed
    sort -k1,1 -k2,2n -k3,3n ~{sample_id}.best_match.unsorted.bed > ~{sample_id}.best_match.bed

    if [ -s ~{sample_id}.best_match.bed ]; then
      cut -f1-3 ~{sample_id}.best_match.bed > ~{sample_id}.best_match.bed3
      bedtools getfasta -fi ~{asm} -bed ~{sample_id}.best_match.bed3 -fo ~{sample_id}.best_match.fa
    else
      echo "Warning: ~{sample_id}: no region assigned to any gene" >&2
      : > ~{sample_id}.best_match.fa
    fi

    echo "INFO: ~{sample_id}: assigned regions per gene:" >&2
    cut -f5 ~{sample_id}.best_match.bed | sort | uniq -c >&2 || true
  >>>

  output {
    File bed = "~{sample_id}.bed"
    File sample_fa = "~{sample_id}.fa"
    File kmer_stats = "~{sample_id}.kmer_stats.tsv"
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

# ---------------------------------------------------------------------------
# NEW TASK: per-haplotype gene copy number table.
#   sample  hap  <gene1>  <gene2> ...
#   S1      h1   1        0
#   S1      h2   2        1
# CN = number of rows with status == PASS for that gene in the haplotype's
# kmer_stats.tsv (= rows of its best_match.bed). Every haplotype gets a row
# and every gene in gene_map.tsv a column, so absent genes show up as 0.
# Rows follow the workflow's sample/hap order; hap1/hap2 are written h1/h2.
# ---------------------------------------------------------------------------
task summarize_cn {
  input {
    Array[File] kmer_stats
    Array[String] sample_ids
    Array[String] haps
    File gene_map
    String cohort_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    stats_arr=(~{sep=" " kmer_stats})
    sample_arr=(~{sep=" " sample_ids})
    hap_arr=(~{sep=" " haps})

    if [ "${#stats_arr[@]}" -ne "${#sample_arr[@]}" ] || [ "${#stats_arr[@]}" -ne "${#hap_arr[@]}" ]; then
      echo "ERROR: kmer_stats / sample_ids / haps have different lengths" >&2
      exit 1
    fi

    : > manifest.tsv
    for i in "${!stats_arr[@]}"; do
      printf "%s\t%s\t%s\n" "${stats_arr[$i]}" "${sample_arr[$i]}" "${hap_arr[$i]}" >> manifest.tsv
    done

    awk -F'\t' -v OFS='\t' '
      FILENAME == ARGV[1] { if (FNR > 1) gene[++ng] = $1; next }
      { nf++; path[nf] = $1; smp[nf] = $2; hp[nf] = $3 }
      END {
        line = "sample" OFS "hap"
        for (g = 1; g <= ng; g++) line = line OFS gene[g]
        print line
        for (i = 1; i <= nf; i++) {
          split("", cn)
          split("", col)
          first = 1
          while ((getline l < path[i]) > 0) {
            n = split(l, f, "\t")
            if (first) {
              for (j = 1; j <= n; j++) col[f[j]] = j
              if (!("gene" in col) || !("status" in col)) {
                print "ERROR: " path[i] " has no gene/status columns" > "/dev/stderr"
                exit 1
              }
              first = 0
              continue
            }
            if (f[col["status"]] == "PASS") cn[f[col["gene"]]]++
          }
          close(path[i])
          h = hp[i]; sub(/^hap/, "h", h)
          line = smp[i] OFS h
          for (g = 1; g <= ng; g++) line = line OFS (cn[gene[g]] + 0)
          print line
        }
      }
    ' ~{gene_map} manifest.tsv > ~{cohort_id}.gene_cn.tsv
  >>>

  output {
    File cn_table = "~{cohort_id}.gene_cn.tsv"
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

# ---------------------------------------------------------------------------
# CHANGED: each assigned sequence is aligned ONLY to the reference of the
# gene it was assigned to (pulled out of panel.fa by region_id), driven by
# best_match.bed rather than by splitting a merged FASTA. hap1 and hap2 are
# processed from their own FASTA files, so identical contig names in the two
# haplotype assemblies cannot collide.
# Raw VCF / sample column name: <sample>_<hap>_<gene>_<contig_start-end>
# ---------------------------------------------------------------------------
task call_variants {
  input {
    Array[File] best_match_fa
    Array[File] best_match_bed
    File panel_fa
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    fa_arr=(~{sep=" " best_match_fa})
    bed_arr=(~{sep=" " best_match_bed})
    cat "${bed_arr[@]}" > ~{sample_id}.bed
    mkdir -p raw_vcfs refs

    for i in "${!fa_arr[@]}"; do
      fa="${fa_arr[$i]}"
      bed="${bed_arr[$i]}"
      while IFS=$'\t' read -r -u 3 contig start end sample_hap gene region_id rest; do
        [ -z "$contig" ] && continue
        ref="refs/${gene}.fa"
        if [ ! -s "$ref" ]; then
          seqkit grep -p "$region_id" ~{panel_fa} > "$ref"
        fi
        seq_key="${contig}:${start}-${end}"
        seqkit grep -p "$seq_key" "$fa" > query.fa
        if [ ! -s query.fa ] || [ ! -s "$ref" ]; then
          echo "ERROR: could not extract ${seq_key} or reference ${region_id} for ${gene}" >&2
          exit 1
        fi
        clean_key=$(echo "$seq_key" | sed 's/[^A-Za-z0-9._-]/_/g')
        clean_id="${sample_hap}_${gene}_${clean_key}"
        minimap2 -c --cs -t 4 "$ref" query.fa | sort -k6,6 -k8,8n | paftools.js call -f "$ref" -L50 -l1 - > raw_vcfs/${clean_id}.raw.vcf
      done 3< "$bed"
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

# Same as before, normalizing against panel.fa (it holds every gene region
# contig) instead of a single gene_seq.
task normalize_merge_vcfs {
  input {
    Array[File] raw_vcfs
    File panel_fa
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    mkdir -p normed_vcfs
    for vcf in ~{sep=" " raw_vcfs}; do
      base=$(basename "$vcf" .raw.vcf)
      printf "sample\t%s\n" "$base" > sample_map.txt
      bcftools reheader -s sample_map.txt "$vcf" | bcftools norm -f ~{panel_fa} -m -any -Oz -o normed_vcfs/${base}.norm.vcf.gz
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
# CHANGED: convert this sample's merged, region-based VCF (one CHROM per gene
# region, "chrN:start-end", POS = 1-based offset in that region) into real
# chromosome coordinates. Every region is looked up in gene_map.tsv (built
# once by prepare_gene_panel), so this task does no FASTA parsing itself.
# The logic lives in gene_to_chrom_vcf_multi.py (passed in as
# `conversion_script`).
# ---------------------------------------------------------------------------
task convert_coordinates {
  input {
    File vcf
    File gene_map
    String sample_id
    File conversion_script

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    bcftools view ~{vcf} > input.vcf

    python ~{conversion_script} input.vcf \
      --map ~{gene_map} \
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
    docker: "eichlerlab/plink-phewas:0.2"
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
# Merge every sample's chromosome-based VCF into one cohort VCF.
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
# NEW TASK: VEP input that can be joined back to the genotypes.
#   1. bcftools norm -m -any: split multiallelic records (bcftools merge
#      joins per-sample biallelic records back into multiallelic ones), so
#      each record = one ALT allele and GT 1 = carrier of that allele.
#   2. ID = CHROM_POS_REF_ALT for every record. VEP reports a VCF's ID as
#      #Uploaded_variation, so lollipop_df can join on it exactly instead of
#      re-deriving VEP's own (indel-trimmed) variant names.
# ---------------------------------------------------------------------------
task prepare_vep_input {
  input {
    File vcf
    String cohort_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    bcftools norm -m -any ~{vcf} -Oz -o split.vcf.gz
    bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' split.vcf.gz -Oz -o ~{cohort_id}.vep_input.vcf.gz
    bcftools index -t ~{cohort_id}.vep_input.vcf.gz

    n_dup=$(bcftools query -f '%ID\n' ~{cohort_id}.vep_input.vcf.gz | sort | uniq -d | wc -l)
    if [ "$n_dup" -gt 0 ]; then
      echo "Warning: ${n_dup} duplicated variant ID(s) in ~{cohort_id}.vep_input.vcf.gz; lollipop_df merges their carriers" >&2
    fi
  >>>

  output {
    File vep_input_vcf = "~{cohort_id}.vep_input.vcf.gz"
    File vep_input_vcf_idx = "~{cohort_id}.vep_input.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          1,
    mem_gb:             4,
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
    String sample_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    cp ~{fasta} ref.fasta
    cp ~{fasta}.fai ref.fasta.fai

    vep -i ~{vcf} -o ~{sample_id}.vep_output.txt \
            --gff ~{gff} \
            --fasta ref.fasta
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

# ---------------------------------------------------------------------------
# NEW TASK: lollipop tables, as rule lollipop_df / lollipop.py in
# SD_genes/paralog_variants_merged, with ALL genes in one file each
# (gene_map.tsv order):
#   <cohort>.lollipop.txt           gene CHROM POS Consequence tag
#                                   Protein_position Amino_acids num hap
#                                   (hap = carriers, ", ")
#   <cohort>.detailed.lollipop.txt  same, one row per carrier copy
# = lollipop.py's columns with a leading gene column.
# num = number of assembled copies (VCF sample columns,
# <sample>_<hap>_<gene>_<contig_start-end>) carrying the ALT allele.
# Logic lives in vep_lollipop.py (passed in as `lollipop_script`).
# ---------------------------------------------------------------------------
task lollipop_df {
  input {
    File vcf
    File vep_output
    File gene_map
    File lollipop_script
    String cohort_id

    RuntimeAttr? runtime_attr_override
  }
  command <<<
    set -euo pipefail

    python ~{lollipop_script} \
      -f ~{vcf} \
      -v ~{vep_output} \
      -g ~{gene_map} \
      -l ~{cohort_id}.lollipop.txt \
      -d ~{cohort_id}.detailed.lollipop.txt
  >>>

  output {
    File lollipop = "~{cohort_id}.lollipop.txt"
    File detailed_lollipop = "~{cohort_id}.detailed.lollipop.txt"
  }

  RuntimeAttr default_attr = object {
    cpu_cores:          1,
    mem_gb:             4,
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  0,
    max_retries:        0,
    docker: "eichlerlab/plink-phewas:0.2"
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
