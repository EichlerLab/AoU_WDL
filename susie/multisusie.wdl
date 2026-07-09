version 1.0

## Multi-population fine-mapping with MultiSuSiE, scattered over (SV, phecode) pairs.
##
## Workflow:
##   1. SplitGenoBySv — single task; pre-splits geno_tsv into one small per-SV TSV
##      (only for SVs actually needed), so the potentially-huge full table is read
##      and localized ONCE instead of once per (SV, phecode) shard.
##   2. RunMultiSuSiE — scattered; wraps multisusie_finemap.py (pull SV + flanking-SNP
##      genotypes, run a per-ancestry logistic-regression scan, feed the summary stats
##      to multisusie_rss).
## See the script's own docstring for full methodology / caveats.
##
## Input TSV (no header): sv_id <TAB> phenotype
##   e.g. chr1-1000-DEL-10000<TAB>phecode_100.0

workflow multisusie_finemap {
    input {
        ## One (sv_id, phenotype) pair per line, tab-separated, no header.
        File sv_phecode_tsv

        ## The multisusie_finemap.py script to run.
        File multisusie_finemap_script

        ## Long-format SV genotype TSV (columns: SVID, Sample, GT)
        File geno_tsv
        ## Covariate file (IID, age, SEX, ancestry column, principal_component_*)
        File covariate_file
        ## Phenotype file (IID + one column per phecode, plink2 1=control/2=case coding)
        File pheno_tsv
        ## e.g. "gs://path/to/CHROM_file.vcf.gz" -- CHROM is replaced with the chromosome
        ## parsed from each sv_id
        String vcf_pattern
        ## Optional file of sample IDs (one per line) to subset the VCF to before pulling
        ## genotypes -- strongly recommended at biobank scale
        File? sample_list

        ## GCP project for requester-pays GCS access to vcf_pattern
        String google_project

        String ancestry_col = "ancestry"
        Array[String] exclude_ancestry = []
        Int window = 10000

        ## Overrides for parse_sv_id -- only needed if sv_id doesn't follow the
        ## chrom-pos-type-length convention
        String? sv_chrom
        String? sv_pos
        String? sv_end

        Float case_val = 2
        Float control_val = 1
        Int min_n = 100
        Int min_cases = 20
        Float single_population_mac_thresh = 20
        Float min_case_mac = 20
        Float multi_population_maf_thresh = 0.01
        Float max_missingness = 0.05
        Int L = 10
        Float rho = 0.75
        Float coverage = 0.95

        String bcftools_bin = "bcftools"

        ## Runtime -- default docker must contain: python3 + pandas/numpy/scipy/
        ## statsmodels/MultiSuSiE, bcftools, and the gcloud CLI (for GCS_OAUTH_TOKEN)
        String docker = "eichlerlab/multisusie:0.1"

        ## SplitGenoBySv resources -- runs once over the full geno_tsv, so disk scales
        ## with its actual size (input + external-sort spill + output splits) rather
        ## than a fixed guess. awk/sort are single-pass/streaming, so this stays cheap
        ## even for a very large geno_tsv.
        Int    prep_cpu          = 2
        Int    prep_memory_gb    = 4
        Int    prep_disk_gb      = ceil(size(geno_tsv, "GB")) * 3 + 20
        Int    prep_preemptible  = 2
        Int    prep_sort_memory_mb = 3072

        ## RunMultiSuSiE resources -- one shard per (SV, phecode) pair, each now only
        ## handling one SV's already-filtered genotypes plus the covariate/phenotype
        ## tables, so this is deliberately lightweight. preemptible=3 is safe because
        ## each shard is idempotent (writes only its own out_prefix files) and short.
        ## disk uses HDD, not SSD -- the task is network/CPU-bound (GCS streaming +
        ## in-memory pandas), not local-disk-IO-bound, and HDD is far cheaper per GB.
        Int    cpu         = 2
        Int    memory_gb   = 8
        Int    disk_gb     = 20
        Int    preemptible = 3
    }

    Array[Array[String]] sv_phecode_pairs = read_tsv(sv_phecode_tsv)

    ## ── Step 1: split geno_tsv into one small per-SV TSV (once, not per shard) ──
    call SplitGenoBySv {
        input:
            geno_tsv        = geno_tsv,
            sv_phecode_tsv  = sv_phecode_tsv,
            sort_memory_mb  = prep_sort_memory_mb,
            cpu             = prep_cpu,
            memory_gb       = prep_memory_gb,
            docker          = docker,
            disk_gb         = prep_disk_gb,
            preemptible     = prep_preemptible
    }

    ## ── Step 2: per-(SV, phecode) scatter ────────────────────────────────────
    ## Indexed rather than a Map[String, File] lookup: SplitGenoBySv.per_pair_geno_tsv
    ## is already aligned 1:1 with sv_phecode_pairs by position (same order, one entry
    ## per input row, duplicates preserved), so no runtime lookup is needed.
    scatter (i in range(length(sv_phecode_pairs))) {
        String shard_sv_id     = sv_phecode_pairs[i][0]
        String shard_phenotype = sv_phecode_pairs[i][1]
        File   shard_geno_tsv  = SplitGenoBySv.per_pair_geno_tsv[i]

        call RunMultiSuSiE {
            input:
                sv_id                        = shard_sv_id,
                phenotype                    = shard_phenotype,
                script                       = multisusie_finemap_script,
                geno_tsv                     = shard_geno_tsv,
                covariate_file               = covariate_file,
                pheno_tsv                    = pheno_tsv,
                vcf_pattern                  = vcf_pattern,
                sample_list                  = sample_list,
                google_project               = google_project,
                ancestry_col                 = ancestry_col,
                exclude_ancestry             = exclude_ancestry,
                window                       = window,
                sv_chrom                     = sv_chrom,
                sv_pos                       = sv_pos,
                sv_end                       = sv_end,
                case_val                     = case_val,
                control_val                  = control_val,
                min_n                        = min_n,
                min_cases                    = min_cases,
                single_population_mac_thresh = single_population_mac_thresh,
                min_case_mac                 = min_case_mac,
                multi_population_maf_thresh  = multi_population_maf_thresh,
                max_missingness              = max_missingness,
                L                            = L,
                rho                          = rho,
                coverage                     = coverage,
                bcftools_bin                 = bcftools_bin,
                docker                       = docker,
                cpu                          = cpu,
                memory_gb                    = memory_gb,
                disk_gb                      = disk_gb,
                preemptible                  = preemptible
        }
    }

    output {
        Array[File] variant_results  = RunMultiSuSiE.variant_results_tsv
        Array[File] multisusie_pkl   = RunMultiSuSiE.multisusie_pkl
    }
}

# ---------------------------------------------------------------------------
task RunMultiSuSiE {
    input {
        String sv_id
        String phenotype
        File   script
        File   geno_tsv
        File   covariate_file
        File   pheno_tsv
        String vcf_pattern
        File?  sample_list
        String google_project
        String ancestry_col
        Array[String] exclude_ancestry
        Int    window
        String? sv_chrom
        String? sv_pos
        String? sv_end
        Float  case_val
        Float  control_val
        Int    min_n
        Int    min_cases
        Float  single_population_mac_thresh
        Float  min_case_mac
        Float  multi_population_maf_thresh
        Float  max_missingness
        Int    L
        Float  rho
        Float  coverage
        String bcftools_bin
        String docker
        Int    cpu
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    ## Sanitize sv_id/phenotype for use as a filename prefix (both can contain
    ## characters like "." that are fine on disk, but keep this defensive against
    ## anything stranger showing up in the input TSV).
    String out_prefix = sub(sv_id, "[^A-Za-z0-9_.-]", "_") + "__" + sub(phenotype, "[^A-Za-z0-9_.-]", "_")

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
        export GCS_REQUESTER_PAYS_PROJECT="~{google_project}"

        python3 "~{script}" \
            --sv-id "~{sv_id}" \
            --phenotype "~{phenotype}" \
            --geno-tsv "~{geno_tsv}" \
            --covariate-file "~{covariate_file}" \
            --pheno-tsv "~{pheno_tsv}" \
            --vcf-pattern "~{vcf_pattern}" \
            --bcftools-bin "~{bcftools_bin}" \
            ~{"--sample-list " + sample_list} \
            --ancestry-col "~{ancestry_col}" \
            ~{sep=" " prefix("--exclude-ancestry ", exclude_ancestry)} \
            --window ~{window} \
            ~{"--sv-chrom " + sv_chrom} \
            ~{"--sv-pos " + sv_pos} \
            ~{"--sv-end " + sv_end} \
            --case-val ~{case_val} \
            --control-val ~{control_val} \
            --min-n ~{min_n} \
            --min-cases ~{min_cases} \
            --single-population-mac-thresh ~{single_population_mac_thresh} \
            --min-case-mac ~{min_case_mac} \
            --multi-population-maf-thresh ~{multi_population_maf_thresh} \
            --max-missingness ~{max_missingness} \
            --L ~{L} \
            --rho ~{rho} \
            --coverage ~{coverage} \
            --out-prefix "~{out_prefix}"
    >>>

    output {
        File variant_results_tsv = "~{out_prefix}_variant_results.tsv"
        File multisusie_pkl      = "~{out_prefix}_multisusie_result.pkl"
    }

    runtime {
        docker:      docker
        cpu:         cpu
        memory:      "~{memory_gb} GB"
        disks:       "local-disk ~{disk_gb} HDD"
        preemptible: preemptible
    }
}

# ---------------------------------------------------------------------------
task SplitGenoBySv {
    input {
        File geno_tsv
        File sv_phecode_tsv
        Int  sort_memory_mb
        Int    cpu
        Int    memory_gb
        String docker
        Int    disk_gb
        Int    preemptible
    }

    command <<<
        set -euo pipefail
        mkdir -p sv_geno

        header=$(head -n1 "~{geno_tsv}")

        ## Unique SV ids actually needed (col 1 of sv_phecode_tsv, deduplicated) --
        ## SVs tested against multiple phecodes are only extracted once.
        cut -f1 "~{sv_phecode_tsv}" | sort -u > unique_sv_ids.txt

        ## Single streaming pass: filter geno_tsv to wanted SVs, sort by SVID (external
        ## sort spills to disk so this scales to arbitrarily large geno_tsv without
        ## loading it into memory), then write one row-group per file, never holding
        ## more than one output file handle open at a time.
        awk -F'\t' 'NR==FNR { want[$1]=1; next } ($1 in want)' \
                unique_sv_ids.txt <(tail -n +2 "~{geno_tsv}") \
            | sort -t $'\t' -k1,1 -S "~{sort_memory_mb}M" \
            | awk -F'\t' -v header="$header" '
                BEGIN { cur = "" }
                {
                    if ($1 != cur) {
                        if (cur != "") close(outfile)
                        cur = $1
                        safe = cur
                        gsub(/[^A-Za-z0-9_.-]/, "_", safe)
                        outfile = "sv_geno/" safe ".tsv"
                        print header > outfile
                    }
                    print $0 > outfile
                }
            '

        missing=""
        while IFS= read -r sv; do
            [ -z "$sv" ] && continue
            safe=$(printf '%s' "$sv" | sed 's/[^A-Za-z0-9_.-]/_/g')
            [[ -s "sv_geno/${safe}.tsv" ]] || missing="${missing}${sv}\n"
        done < unique_sv_ids.txt
        if [[ -n "$missing" ]]; then
            echo "ERROR: the following sv_id(s) from sv_phecode_tsv have no rows in geno_tsv:" >&2
            printf "%b" "$missing" >&2
            exit 1
        fi

        : > ordered_paths.txt
        while IFS=$'\t' read -r sv _rest; do
            safe=$(printf '%s' "$sv" | sed 's/[^A-Za-z0-9_.-]/_/g')
            printf '%s\n' "$PWD/sv_geno/${safe}.tsv" >> ordered_paths.txt
        done < "~{sv_phecode_tsv}"
    >>>

    output {
        Array[File] per_pair_geno_tsv = read_lines("ordered_paths.txt")
    }

    runtime {
        docker:      docker
        cpu:         cpu
        memory:      "~{memory_gb} GB"
        disks:       "local-disk ~{disk_gb} HDD"
        preemptible: preemptible
    }
}
