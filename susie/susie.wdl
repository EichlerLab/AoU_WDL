version 1.0

## Fine-mapping pipeline: per-SV scatter for susieR.
##
## Input TSV (no header): sv_id <TAB> phenotype
## Workflow:
##   1. PrepVcfs        — single task; splits geno TSV into one VCF per SV
##   2. PrepMergedPgen  — scattered; fetches flanking SNPs, merges with SV GT, builds pgen
##   3. RunAssociation  — scattered; plink2 Firth → z-scores for all variants
##   4. RunSusie        — scattered; susie_rss(z, R, n) → PIPs + credible sets

workflow susieR_finemap_prep {
    input {
        ## TSV of significant SV–phenotype pairs (no header): sv_id <TAB> phenotype
        File bonferroni_hits_tsv

        ## Long-format genotype TSV (header: SVID\tSample\tGT)
        File sv_geno_tsv

        ## Scripts + supporting files reused from the phewas pipeline
        File prep_vcfs_script
        File relatedness_flagged_samples

        ## Phenotype + covariate tables
        File phenotype_file
        File covariate_file

        ## Covariates passed to --covar-name in plink2
        String covar_names = "age,SEX,principal_component_1,principal_component_2,principal_component_3,principal_component_4,principal_component_5,principal_component_6,principal_component_7,principal_component_8,principal_component_9,principal_component_10"

        ## Flanking window around SV start (bp)
        Int flank_bp = 100000

        ## Minimum AC for SNP inclusion
        Int min_ac = 10

        ## susieR: maximum number of causal signals per locus (L=10 is standard for GWAS)
        Int susie_L = 10

        ## Resources — PrepVcfs needs more memory (full geno TSV in RAM)
        Int    prep_cpu        = 1
        Int    prep_memory_gb  = 32
        Int    prep_disk_gb    = 50
        Int    prep_preemptible = 1

        Int    cpu         = 2
        Int    memory_gb   = 32
        String docker      = "eichlerlab/plink-phewas:0.2"
        Int    disk_gb     = 50
        Int    preemptible = 1
    }

    ## ── Step 1: split geno TSV into per-SV VCFs (single task) ───────────────
    call PrepVcfs {
        input:
            bonferroni_hits_tsv        = bonferroni_hits_tsv,
            sv_geno_tsv                = sv_geno_tsv,
            prep_vcfs_script           = prep_vcfs_script,
            relatedness_flagged_samples = relatedness_flagged_samples,
            phenotype_file             = phenotype_file,
            cpu                        = prep_cpu,
            memory_gb                  = prep_memory_gb,
            docker                     = docker,
            disk_gb                    = prep_disk_gb,
            preemptible                = prep_preemptible
    }

    ## ── Step 2–4: per-SV scatter ─────────────────────────────────────────────
    scatter (i in range(length(PrepVcfs.vcf_files))) {
        File   sv_gt_vcf = PrepVcfs.vcf_files[i]
        String sv_id     = PrepVcfs.sv_ids[i]
        String phenotype = PrepVcfs.phenotypes[i]

        call PrepMergedPgen {
            input:
                sv_id     = sv_id,
                sv_gt_vcf = sv_gt_vcf,
                flank_bp  = flank_bp,
                min_ac    = min_ac,
                cpu       = cpu,
                memory_gb = memory_gb,
                docker    = docker,
                disk_gb   = disk_gb,
                preemptible = preemptible
        }

        call RunAssociation {
            input:
                sv_id          = sv_id,
                phenotype      = phenotype,
                merged_pgen    = PrepMergedPgen.merged_pgen,
                merged_pvar    = PrepMergedPgen.merged_pvar,
                merged_psam    = PrepMergedPgen.merged_psam,
                phenotype_file = phenotype_file,
                covariate_file = covariate_file,
                covar_names    = covar_names,
                cpu            = cpu,
                memory_gb      = memory_gb,
                docker         = docker,
                disk_gb        = disk_gb,
                preemptible    = preemptible
        }

        call RunSusie {
            input:
                sv_id       = sv_id,
                merged_pgen = PrepMergedPgen.merged_pgen,
                merged_pvar = PrepMergedPgen.merged_pvar,
                merged_psam = PrepMergedPgen.merged_psam,
                assoc_tsv   = RunAssociation.assoc_results,
                susie_L     = susie_L,
                cpu         = cpu,
                memory_gb   = memory_gb,
                docker      = docker,
                disk_gb     = disk_gb,
                preemptible = preemptible
        }
    }

    output {
        Array[File] snp_tsv       = PrepMergedPgen.snp_tsv
        Array[File] assoc_results = RunAssociation.assoc_results
        Array[File] susie_pips    = RunSusie.susie_pips
        Array[File] susie_cs      = RunSusie.susie_cs
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
            --sv-list-out                 sv_list.txt \
            --sv-metadata-out             sv_metadata.tsv

        ## Write sv_ids and phenotypes in alphabetical order (matching glob("*_GT.vcf"))
        python3 - <<'PYEOF'
bonferroni_hits = "~{bonferroni_hits_tsv}"

with open("sv_list.txt") as f:
    successful = set(line.strip().replace("_GT.vcf", "") for line in f if line.strip())

sv_to_pheno = {}
with open(bonferroni_hits) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2 and parts[0] in successful:
            sv_to_pheno[parts[0]] = parts[1]

ordered = sorted(sv_to_pheno)
with open("sv_ids_ordered.txt", "w") as f:
    f.write("\n".join(ordered) + "\n")
with open("phenotypes_ordered.txt", "w") as f:
    f.write("\n".join(sv_to_pheno[sv] for sv in ordered) + "\n")
PYEOF
    >>>

    output {
        ## glob returns files in alphabetical order; sv_ids/phenotypes are written
        ## in the same alphabetical order so the three arrays are aligned.
        Array[File]   vcf_files  = glob("*_GT.vcf")
        Array[String] sv_ids     = read_lines("sv_ids_ordered.txt")
        Array[String] phenotypes = read_lines("phenotypes_ordered.txt")
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

        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)

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

        bcftools view --threads ~{cpu} -r "$region" -S "${sv}_samples.txt" "$snp_vcf_path" \
            | bcftools norm --threads ~{cpu} -m - \
            | bcftools view --threads ~{cpu} -m2 -M2 -v snps \
            | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
            | bcftools filter -i "AC > ~{min_ac}" \
            -Oz -o "${sv}_SNP.vcf.gz"
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
    >>>

    output {
        File snp_tsv     = "~{sv_id}_SNP.tsv"
        File merged_pgen = "~{sv_id}_merged.pgen"
        File merged_pvar = "~{sv_id}_merged.pvar"
        File merged_psam = "~{sv_id}_merged.psam"
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
task RunAssociation {
    input {
        String sv_id
        String phenotype
        File   merged_pgen
        File   merged_pvar
        File   merged_psam
        File   phenotype_file
        File   covariate_file
        String covar_names
        Int    cpu
        Int    memory_gb
        String docker
        Int    disk_gb
        Int    preemptible
    }

    command <<<
        set -euo pipefail

        sv="~{sv_id}"

        ln -sf "~{merged_pgen}" merged.pgen
        ln -sf "~{merged_pvar}" merged.pvar
        ln -sf "~{merged_psam}" merged.psam

        plink2 \
            --pfile merged \
            --pheno "~{phenotype_file}" \
            --pheno-name ~{phenotype} \
            --covar "~{covariate_file}" \
            --covar-name ~{covar_names} \
            --covar-variance-standardize \
            --glm firth hide-covar \
            --threads ~{cpu} \
            --out "${sv}"

        mv "${sv}.~{phenotype}.glm.firth" "${sv}_assoc.tsv"
    >>>

    output {
        File assoc_results = "~{sv_id}_assoc.tsv"
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
        File   merged_pgen
        File   merged_pvar
        File   merged_psam
        File   assoc_tsv
        Int    susie_L
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

        sv <- "~{sv_id}"

        ## ── 1. Compute in-sample LD correlation matrix from merged pgen ───
        pvar        <- NewPvar("merged.pvar")
        pgen        <- NewPgen("merged.pgen", pvar = pvar)
        n_variants  <- GetVariantCt(pvar)
        variant_ids <- sapply(seq_len(n_variants), function(i) GetVariantId(pvar, i))

        X <- ReadList(pgen, variant_subset = seq_len(n_variants), meanimpute = TRUE)
        R <- cor(X)
        rownames(R) <- colnames(R) <- variant_ids

        ## ── 2. Parse z-scores and sample size from plink2 Firth output ────
        ## Columns: #CHROM POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P
        assoc <- read.table("~{assoc_tsv}", header = TRUE, sep = "\t",
                            comment.char = "", check.names = FALSE, na.strings = "NA")
        colnames(assoc)[1] <- "CHROM"
        assoc <- assoc[!is.na(assoc$BETA) & !is.na(assoc$SE), ]

        n <- as.integer(median(assoc$OBS_CT))
        z <- setNames(assoc$BETA / assoc$SE, assoc$ID)

        ## ── 3. Align variant order between LD matrix and z-scores ─────────
        shared <- intersect(variant_ids, names(z))
        R <- R[shared, shared]
        z <- z[shared]

        ## ── 4. Run susie_rss ──────────────────────────────────────────────
        fit <- susie_rss(z = z, R = R, n = n,
                         L = ~{susie_L},
                         estimate_residual_variance = TRUE)

        ## ── 5. Write PIPs ─────────────────────────────────────────────────
        pip_df <- data.frame(variant_id = shared,
                             pip        = susie_get_pip(fit),
                             stringsAsFactors = FALSE)
        pip_df <- pip_df[order(-pip_df$pip), ]
        write.table(pip_df, paste0(sv, "_susie_pips.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)

        ## ── 6. Write 95% credible sets ────────────────────────────────────
        cs <- susie_get_cs(fit, coverage = 0.95)
        if (length(cs$cs) > 0) {
            cs_rows <- lapply(names(cs$cs), function(nm) {
                idx <- cs$cs[[nm]]
                data.frame(cs_id      = nm,
                           variant_id = shared[idx],
                           pip        = susie_get_pip(fit)[idx],
                           stringsAsFactors = FALSE)
            })
            cs_df <- do.call(rbind, cs_rows)
        } else {
            cs_df <- data.frame(cs_id      = character(0),
                                variant_id = character(0),
                                pip        = numeric(0))
        }
        write.table(cs_df, paste0(sv, "_susie_cs.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)

        cat(sprintf("susieR done: %d variants, %d credible set(s)\n",
                    length(shared), length(cs$cs)))
        REOF
    >>>

    output {
        File susie_pips = "~{sv_id}_susie_pips.tsv"
        File susie_cs   = "~{sv_id}_susie_cs.tsv"
    }

    runtime {
        docker:      docker
        cpu:         cpu
        memory:      "~{memory_gb} GB"
        disks:       "local-disk ~{disk_gb} SSD"
        preemptible: preemptible
    }
}
