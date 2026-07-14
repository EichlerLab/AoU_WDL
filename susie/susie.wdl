version 1.0

## Fine-mapping pipeline: per-SV scatter for susieR.
##
## Input TSV (no header): sv_id <TAB> phecode_1 phecode_2 ... (space-separated phecodes).
## Workflow:
##   1. PrepVcfs          — single task; splits geno TSV into one VCF per unique SV
##   2. PrepMergedPgen    — scattered per unique SV (outer scatter); fetches flanking SNPs,
##                          merges with SV GT, builds pgen, and computes pairwise unphased
##                          LD (r^2) across the merged set — runs once per SV regardless of
##                          how many phecodes it has
##   3. RunSusie          — scattered per phecode within each SV (inner scatter), reusing
##                          that SV's pgen; susie(X, y) → PIPs + credible sets
##   4. MergeSusieResults — single task; gathers all per-(sv, phecode) susie outputs
##                          into one directory

workflow susieR_finemap_prep {
    input {
        ## Pre-computed per-SV VCFs — if provided, PrepVcfs is skipped entirely.
        ## A text file with one VCF path per line, parallel (same order) to the phecode
        ## groupings in precomputed_phenotypes_grouped_file, a JSON array of arrays (one
        ## inner array of phecodes per VCF), e.g. [["T2D","Obesity"],["Asthma"]].
        File? precomputed_vcf_files_list
        File? precomputed_phenotypes_grouped_file

        ## Pre-built merged pgen/pvar/psam per SV — if provided, PrepMergedPgen (including its
        ## GCS SNP fetch) is skipped entirely for every SV. Each is a text file with one path
        ## per line, parallel to each other and to vcf_files/phenotypes_grouped (i.e. sorted
        ## by sv_id). Useful to re-run RunSusie after a PrepMergedPgen shard was lost to
        ## preemption elsewhere in the scatter, without re-paying for the GCS fetch on shards
        ## that already succeeded.
        File? precomputed_merged_pgens_list
        File? precomputed_merged_pvars_list
        File? precomputed_merged_psams_list

        ## PrepVcfs inputs — bonferroni_hits_tsv always required; others only needed
        ## when precomputed VCFs are not supplied
        File  bonferroni_hits_tsv
        File? sv_geno_tsv
        File? prep_vcfs_script
        File? relatedness_flagged_samples

        ## Phenotype + covariate tables
        File phenotype_file
        File covariate_file

        ## Covariates passed to --covar-name in plink2
        String covar_names = "age,SEX,principal_component_1,principal_component_2,principal_component_3,principal_component_4,principal_component_5,principal_component_6,principal_component_7,principal_component_8,principal_component_9,principal_component_10"

        ## Flanking window around SV start (bp)
        Int flank_bp = 100000

        ## Minimum AC for SNP inclusion
        Int min_ac = 100

        ## Minimum number of cases carrying a variant for it to be included in susie —
        ## drops variants whose signal hinges on a handful of case carriers
        Int min_case_carriers = 5

        ## susieR: maximum number of causal signals per locus (L=10 is standard for GWAS)
        Int   susie_L        = 10
        ## susieR: credible set coverage
        Float susie_coverage = 0.8

        ## GCP project for requester-pays GCS access (needed by PrepMergedPgen)
        String google_project

        ## Resources — PrepVcfs needs more memory (full geno TSV in RAM)
        Int    prep_cpu         = 1
        Int    prep_memory_gb   = 32
        Int    prep_disk_gb     = 50
        Int    prep_preemptible = 1

        Int    cpu         = 8
        Int    memory_gb   = 85
        String docker      = "eichlerlab/plink-phewas:0.2"
        Int    disk_gb     = 50
        Int    preemptible = 1
    }

    if (defined(precomputed_vcf_files_list)) {
        Array[File] precomputed_vcf_files = read_lines(select_first([precomputed_vcf_files_list]))
    }

    ## ── Step 1: split geno TSV into per-SV VCFs (skipped if precomputed supplied) ──
    if (!defined(precomputed_vcf_files_list)) {
        call PrepVcfs {
            input:
                bonferroni_hits_tsv         = select_first([bonferroni_hits_tsv]),
                sv_geno_tsv                 = select_first([sv_geno_tsv]),
                prep_vcfs_script            = select_first([prep_vcfs_script]),
                relatedness_flagged_samples = select_first([relatedness_flagged_samples]),
                phenotype_file              = phenotype_file,
                cpu                         = prep_cpu,
                memory_gb                   = prep_memory_gb,
                docker                      = docker,
                disk_gb                     = prep_disk_gb,
                preemptible                 = prep_preemptible
        }
    }

    if (defined(precomputed_phenotypes_grouped_file)) {
        Array[Array[String]] precomputed_phenotypes_grouped = read_json(select_first([precomputed_phenotypes_grouped_file]))
    }

    Array[File]           vcf_files          = select_first([precomputed_vcf_files,          PrepVcfs.vcf_files])
    Array[Array[String]]  phenotypes_grouped = select_first([precomputed_phenotypes_grouped, PrepVcfs.phenotypes_grouped])

    if (defined(precomputed_merged_pgens_list)) {
        Array[File] precomputed_merged_pgens = read_lines(select_first([precomputed_merged_pgens_list]))
        Array[File] precomputed_merged_pvars = read_lines(select_first([precomputed_merged_pvars_list]))
        Array[File] precomputed_merged_psams = read_lines(select_first([precomputed_merged_psams_list]))
    }

    Boolean skip_prep_merged_pgen = defined(precomputed_merged_pgens_list)

    ## ── Step 2–4: outer scatter per unique SV, inner scatter per phecode ──────
    scatter (i in range(length(vcf_files))) {
        File   sv_gt_vcf = vcf_files[i]
        String sv_id     = basename(sv_gt_vcf, "_GT.vcf")

        if (!skip_prep_merged_pgen) {
            call PrepMergedPgen {
                input:
                    sv_id          = sv_id,
                    sv_gt_vcf      = sv_gt_vcf,
                    google_project = google_project,
                    flank_bp       = flank_bp,
                    min_ac         = min_ac,
                    cpu            = cpu,
                    memory_gb      = memory_gb,
                    docker         = docker,
                    disk_gb        = disk_gb,
                    preemptible    = preemptible
            }
        }

        File this_merged_pgen = if skip_prep_merged_pgen then select_first([precomputed_merged_pgens])[i] else select_first([PrepMergedPgen.merged_pgen])
        File this_merged_pvar = if skip_prep_merged_pgen then select_first([precomputed_merged_pvars])[i] else select_first([PrepMergedPgen.merged_pvar])
        File this_merged_psam = if skip_prep_merged_pgen then select_first([precomputed_merged_psams])[i] else select_first([PrepMergedPgen.merged_psam])

        scatter (phenotype in phenotypes_grouped[i]) {
            call RunSusie {
                input:
                    sv_id          = sv_id,
                    phenotype      = phenotype,
                    merged_pgen    = this_merged_pgen,
                    merged_pvar    = this_merged_pvar,
                    merged_psam    = this_merged_psam,
                    phenotype_file = phenotype_file,
                    covariate_file = covariate_file,
                    covar_names        = covar_names,
                    min_case_carriers  = min_case_carriers,
                    susie_L            = susie_L,
                    susie_coverage     = susie_coverage,
                    cpu                = cpu,
                    memory_gb          = memory_gb,
                    docker             = docker,
                    disk_gb            = disk_gb,
                    preemptible       = preemptible
            }
        }
    }

    Array[File] susie_out_flat = flatten(RunSusie.susie_out)

    ## ── Step 5: gather every susie result into one output directory ──────────
    call MergeSusieResults {
        input:
            susie_files = susie_out_flat,
            docker      = docker
    }

    output {
        ## Empty when precomputed_merged_pgens is supplied (PrepMergedPgen skipped, so no new
        ## snp_tsv/LD files are produced — the caller already has them from the prior run).
        Array[File] snp_tsv          = select_all(PrepMergedPgen.snp_tsv)
        Array[File] ld_files         = select_all(PrepMergedPgen.merged_ld)
        Array[File] susie_out        = susie_out_flat
        File        susie_merged_dir = MergeSusieResults.merged_dir
    }
}

# ---------------------------------------------------------------------------
task PrepVcfs {
    input {
        File   bonferroni_hits_tsv
        File   sv_geno_tsv
        File   prep_vcfs_script
        File   relatedness_flagged_samples
        File   phenotype_file
        Int    cpu
        Int    memory_gb
        String docker
        Int    disk_gb
        Int    preemptible
    }

    command <<<
        set -euo pipefail

        ## Build comma-separated SV subset from bonferroni_hits_tsv (col 1, no header)
        sv_subset=$(awk '{print $1}' "~{bonferroni_hits_tsv}" | paste -sd,)

        ## Generate per-SV VCFs using the shared script.
        ## --min-genotyped-samples 0: skip count threshold (SVs already validated upstream).
        python3 "~{prep_vcfs_script}" \
            --geno-tsv                    "~{sv_geno_tsv}" \
            --pheno-tsv                   "~{phenotype_file}" \
            --relatedness-flagged-samples "~{relatedness_flagged_samples}" \
            --sv-subset                   "$sv_subset" \
            --min-genotyped-samples       0 \
            --output-dir                  . \
            --sv-list-out                 sv_list.txt

        ## Group each successfully-generated SV's phecodes (col 2, space-separated) into
        ## a JSON array of arrays, in the same alphabetical order glob("*_GT.vcf") returns.
        python3 - <<'PYEOF'
import json

bonferroni_hits = "~{bonferroni_hits_tsv}"

with open("sv_list.txt") as f:
    successful = set(line.strip().replace("_GT.vcf", "") for line in f if line.strip())

sv_to_phenos = {}
with open(bonferroni_hits) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2 and parts[0] in successful:
            sv_to_phenos[parts[0]] = parts[1].split()

unique_svs = sorted(sv_to_phenos)
with open("phenotypes_grouped.json", "w") as f:
    json.dump([sv_to_phenos[sv] for sv in unique_svs], f)
PYEOF
    >>>

    output {
        Array[File]           vcf_files          = glob("*_GT.vcf")
        Array[Array[String]]  phenotypes_grouped = read_json("phenotypes_grouped.json")
    }

    runtime {
        docker:      docker
        cpu:         cpu
        memory:      "~{memory_gb} GB"
        disks:       "local-disk ~{disk_gb} SSD"
        preemptible: preemptible
    }
}

# ---------------------------------------------------------------------------
task PrepMergedPgen {
    input {
        String sv_id
        File   sv_gt_vcf   # per-SV VCF from PrepVcfs
        String google_project
        Int    flank_bp
        Int    min_ac
        Int    cpu
        Int    memory_gb
        String docker
        Int    disk_gb
        Int    preemptible
    }

    command <<<
        set -euo pipefail

        sv="~{sv_id}"

        export GCS_REQUESTER_PAYS_PROJECT="~{google_project}"

        ## Extract sample list from VCF header for bcftools -S
        bcftools query -l "~{sv_gt_vcf}" > "${sv}_samples.txt"

        ## Parse chrom / start / type / len from sv_id (chrN-POS-TYPE-LEN)
        IFS="-" read -ra _parts <<< "$sv"
        chrom="${_parts[0]}"
        pos="${_parts[1]}"
        sv_type="${_parts[2]}"
        sv_len="${_parts[3]}"
        start=$(( pos - ~{flank_bp} ))
        if [[ "$sv_type" == "DEL" ]]; then
            end=$(( pos + sv_len + ~{flank_bp} ))
        else
            end=$(( pos + ~{flank_bp} ))
        fi
        region="${chrom}:${start}-${end}"

        ## 1. Fetch + filter flanking SNPs from GCS phased VCF
        snp_vcf_path="gs://vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/phasing/${chrom}_AOU_v8.2_allsamples_phased.vcf.gz"

        ## bcftools' GCS reads occasionally fail mid-stream with a libcurl
        ## HTTP/2 framing error (curl error 92). Retry with new token
        fetch_snps() {
            export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
            bcftools view --threads ~{cpu} -r "$region" -S "${sv}_samples.txt" "$snp_vcf_path" \
                | bcftools norm --threads ~{cpu} -m - \
                | bcftools view --threads ~{cpu} -m2 -M2 -v snps \
                | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
                | bcftools filter -i "AC > ~{min_ac}" \
                -Oz -o "${sv}_SNP.vcf.gz"
        }

        max_attempts=5
        for attempt in $(seq 1 "$max_attempts"); do
            if fetch_snps; then
                break
            fi
            if [[ "$attempt" -eq "$max_attempts" ]]; then
                echo "GCS fetch failed after ${max_attempts} attempts" >&2
                exit 1
            fi
            echo "GCS fetch attempt ${attempt} failed, retrying in 10s..." >&2
            sleep 10
        done
        bcftools index --threads ~{cpu} -t "${sv}_SNP.vcf.gz"

        ## 2. SNP summary table (MAF / AC)
        bcftools +fill-tags "${sv}_SNP.vcf.gz" \
            | bcftools query -f "%ID\t%INFO/MAF\t%INFO/AC\n" \
            > "${sv}_SNP.tsv"

        ## 3. Bgzip + index SV GT VCF
        bgzip -@ ~{cpu} -c "~{sv_gt_vcf}" > "${sv}_GT.vcf.gz"
        tabix -p vcf "${sv}_GT.vcf.gz"

        ## 4. Merge SV + flanking SNPs → single pgen for association and susieR
        bcftools concat -a "${sv}_GT.vcf.gz" "${sv}_SNP.vcf.gz" -Ou \
            | bcftools sort -Oz -o "${sv}_merged.vcf.gz"
        tabix -p vcf "${sv}_merged.vcf.gz"

        plink2 --vcf "${sv}_merged.vcf.gz" --make-pgen --out "${sv}_merged"

        ## 5. Unphased pairwise LD (r^2) across the merged SV + flanking-SNP set
        plink2 --pfile "${sv}_merged" --r2-unphased --ld-window-r2 0 --out "${sv}_ld"
    >>>

    output {
        File snp_tsv     = "~{sv_id}_SNP.tsv"
        File merged_pgen = "~{sv_id}_merged.pgen"
        File merged_pvar = "~{sv_id}_merged.pvar"
        File merged_psam = "~{sv_id}_merged.psam"
        File merged_ld   = "~{sv_id}_ld.vcor"
    }

    runtime {
        docker:      docker
        cpu:         cpu
        memory:      "~{memory_gb} GB"
        disks:       "local-disk ~{disk_gb} SSD"
        preemptible: preemptible
    }
}

# ---------------------------------------------------------------------------
task RunSusie {
    input {
        String sv_id
        String phenotype
        File   merged_pgen
        File   merged_pvar
        File   merged_psam
        File   phenotype_file
        File   covariate_file
        String covar_names = "age,SEX,principal_component_1,principal_component_2,principal_component_3,principal_component_4,principal_component_5,principal_component_6,principal_component_7,principal_component_8,principal_component_9,principal_component_10"
        Int    min_case_carriers
        Int    susie_L
        Float  susie_coverage
        Int    cpu
        Int    memory_gb
        String docker
        Int    disk_gb
        Int    preemptible
    }

    command <<<
        set -euo pipefail

        ln -sf "~{merged_pgen}" merged.pgen
        ln -sf "~{merged_pvar}" merged.pvar
        ln -sf "~{merged_psam}" merged.psam

        Rscript - <<'REOF'
        library(pgenlibr)
        library(susieR)

        sv        <- "~{sv_id}"
        phecode   <- "~{phenotype}"
        covar_vec <- strsplit("~{covar_names}", ",")[[1]]

        ## ── 1. Read genotypes from merged pgen (SV + flanking SNPs) ──────────
        pvar        <- NewPvar("merged.pvar")
        pgen        <- NewPgen("merged.pgen", pvar = pvar)
        n_variants  <- GetVariantCt(pvar)
        variant_ids <- sapply(seq_len(n_variants), function(i) GetVariantId(pvar, i))
        X_full      <- ReadList(pgen, variant_subset = seq_len(n_variants), meanimpute = TRUE)
        colnames(X_full) <- variant_ids

        ## plink2 files use a leading '#' comment char on the header line —
        ## read.table's defaults (comment.char = "#", check.names = TRUE)
        ## mangle it, so read manually and strip the '#'.
        read_plink_table <- function(path) {
            df <- read.table(path, header = TRUE, sep = "\t",
                             comment.char = "", check.names = FALSE)
            colnames(df)[1] <- sub("^#", "", colnames(df)[1])
            df
        }

        psam       <- read_plink_table("merged.psam")
        sample_ids <- as.character(psam$IID)

        ## ── 2. Merge phenotype + covariates, align to pgen sample order ───────
        pheno_dat <- read_plink_table("~{phenotype_file}")
        covar_dat <- read_plink_table("~{covariate_file}")

        missing_cols <- setdiff(c("IID", phecode), colnames(pheno_dat))
        if (length(missing_cols) > 0) {
            stop(sprintf(
                "phenotype_file is missing expected column(s): %s\nphecode requested: [%s]\ncolumns found: %s",
                paste(missing_cols, collapse = ", "),
                phecode,
                paste(sprintf("[%s]", colnames(pheno_dat)), collapse = ", ")
            ))
        }

        pheno_sub <- pheno_dat[, c("IID", phecode)]
        dat       <- merge(pheno_sub, covar_dat, by = "IID")
        dat       <- dat[!is.na(dat[[phecode]]), ]

        ## phecode columns use the standard 1 = control / 2 = case coding —
        ## recode to 0/1 for glm(family = binomial).
        dat[[phecode]] <- dat[[phecode]] - 1

        ## Restrict to samples that actually have genotype data before indexing,
        ## otherwise unmatched IIDs introduce NA rows into X.
        dat <- dat[dat$IID %in% sample_ids, ]

        row_idx <- match(dat$IID, sample_ids)
        X_raw   <- X_full[row_idx, , drop = FALSE]

        ## ── 3. Drop variants with too few case carriers, or with no variance
        ## in this subset (susie's scale() would turn those into NaN columns).
        case_status      <- dat[[phecode]]
        is_carrier_full  <- round(X_raw) >= 1
        n_case_with_full <- colSums(is_carrier_full[case_status == 1, , drop = FALSE], na.rm = TRUE)
        has_variance     <- apply(X_raw, 2, function(col) length(unique(col[!is.na(col)])) > 1)
        keep             <- (n_case_with_full >= ~{min_case_carriers}) & has_variance

        message(sprintf(
            "dropping %d/%d variant(s): %d with fewer than %d case carrier(s), %d monomorphic in this subset: %s",
            sum(!keep), length(keep),
            sum(n_case_with_full < ~{min_case_carriers}), ~{min_case_carriers},
            sum(!has_variance), paste(variant_ids[!has_variance], collapse = ", ")
        ))

        variant_ids <- variant_ids[keep]
        X_raw       <- X_raw[, keep, drop = FALSE]
        n_variants  <- length(variant_ids)
        X           <- scale(X_raw)

        ## ── 4. Null model residuals ───────────────────────────────────────────
        null_formula <- as.formula(paste(phecode, "~", paste(covar_vec, collapse = " + ")))
        null <- glm(null_formula, family = binomial, data = dat)
        y    <- residuals(null, type = "response")

        ## ── 5. Run susie ──────────────────────────────────────────────────────
        fit <- susie(X, y, L = ~{susie_L}, coverage = ~{susie_coverage})

        message(sprintf(
            "susie found %d credible set(s); target SV '%s' %s in variant_ids (n_variants=%d, max_pip=%.4g)",
            length(fit$sets$cs), sv,
            ifelse(sv %in% variant_ids, "IS present", "is NOT present"),
            n_variants, max(fit$pip)
        ))

        ## CS: credible set number a variant belongs to (NA if none).
        ## Lead_in_CS: TRUE for the highest-PIP variant within its own credible set.
        cs_membership <- rep(NA_integer_, n_variants)
        lead_in_cs    <- rep(FALSE, n_variants)
        for (j in seq_along(fit$sets$cs)) {
            idx <- fit$sets$cs[[j]]
            cs_membership[idx] <- j
            lead_in_cs[idx[which.max(fit$pip[idx])]] <- TRUE
        }

        ## Per-variant case/control counts by carrier status (>= 1 copy of the
        ## alt allele, rounding the mean-imputed dosage to the nearest hard call).
        is_carrier        <- round(X_raw) >= 1
        n_case_with       <- colSums(is_carrier[case_status == 1, , drop = FALSE])
        n_case_without    <- colSums(!is_carrier[case_status == 1, , drop = FALSE])
        n_control_with    <- colSums(is_carrier[case_status == 0, , drop = FALSE])
        n_control_without <- colSums(!is_carrier[case_status == 0, , drop = FALSE])

        final_df <- data.frame(
            Variant           = variant_ids,
            PIP               = fit$pip,
            CS                = cs_membership,
            Lead_in_CS        = lead_in_cs,
            Cases_W_Var       = n_case_with,
            Cases_WO_Var    = n_case_without,
            Ctrls_W_Var    = n_control_with,
            Ctrls_WO_Var = n_control_without,
            stringsAsFactors = FALSE
        )
        final_df <- final_df[order(-final_df$PIP), ]

        ## Phenotype suffix keeps filenames unique when an SV has more than one phenotype.
        write.table(final_df, paste0(sv, "_", phecode, "_susie.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        REOF
    >>>

    output {
        File susie_out = "~{sv_id}_~{phenotype}_susie.tsv"
    }

    runtime {
        docker:      docker
        cpu:         cpu
        memory:      "~{memory_gb} GB"
        disks:       "local-disk ~{disk_gb} SSD"
        preemptible: preemptible
    }
}

# ---------------------------------------------------------------------------
task MergeSusieResults {
    input {
        Array[File] susie_files
        String      docker
        Int         cpu         = 1
        Int         memory_gb   = 4
        Int         disk_gb     = 20
        Int         preemptible = 1
    }

    command <<<
        set -euo pipefail

        ## Packaged as a tarball since WDL/Cromwell task outputs must be files, not directories.
        mkdir -p susie_results
        while read -r f; do
            cp "$f" susie_results/
        done < "~{write_lines(susie_files)}"

        tar -czf susie_results.tar.gz susie_results
    >>>

    output {
        File merged_dir = "susie_results.tar.gz"
    }

    runtime {
        docker:      docker
        cpu:         cpu
        memory:      "~{memory_gb} GB"
        disks:       "local-disk ~{disk_gb} SSD"
        preemptible: preemptible
    }
}
