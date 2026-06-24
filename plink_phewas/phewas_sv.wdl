version 1.0

## PheWAS across structural variants using plink2 Firth regression.
##
## Workflow:
##   1. PrepVcfs   – one task; writes one VCF per passing SV
##   2. RunPhewas  – scattered over SVs; plink2 GLM → format → Manhattan plot
##

workflow PhewasAcrossSVs {

  input {
    # ── Input data files (gs:// paths) ──────────────────────────────────────
    File geno_tsv                     # long-format genotypes:  SVID \t Sample \t GT
    File relatedness_flagged_samples  # one sample per line (first line = header)
    File pheno_tsv                    # plink2 phenotype file (also used to validate sample set)
    File covars_tsv                   # plink2 covariates file
    File sex_file                     # plink2 sex file (FID IID SEX)
    File phecode_definitions          # phecode_definitions1.2.csv

    # ── Python helper scripts (gs:// paths or local) ─────────────────────────
    File prep_vcfs_script
    File format_results_script
    File plot_phewas_script

    # ── Analysis parameters ──────────────────────────────────────────────────
    # Optional comma-separated list of SV IDs to process (e.g. "chr1-100-DEL-1,chr2-200-INS-2").
    # Leave empty ("") to process all passing SVs.
    String sv_subset = ""

    Array[String] covar_names = [
      "age", "SEX",
      "principal_component_1", "principal_component_2", "principal_component_3",
      "principal_component_4",  "principal_component_5",  "principal_component_6",
      "principal_component_7",  "principal_component_8",  "principal_component_9",
      "principal_component_10"
    ]
    Int     min_genotyped_samples  = 200000
    Int     plink_threads          = 4
    Boolean label_only_bonferroni  = true  # if false, label top manhattan_label_count instead
    Int     manhattan_label_count  = 15

    # ── Docker image ─────────────────────────────────────────────────────────
    String docker = "eichlerlab/plink-phewas:0.1"

    # ── Runtime overrides ────────────────────────────────────────────────────
    Int prep_memory_gb   = 32
    Int prep_cpu         = 1
    Int prep_disk_gb     = 50
    Int prep_preemptible = 1

    Int phewas_memory_gb   = 16
    Int phewas_cpu         = 4
    Int phewas_disk_gb     = 30
    Int phewas_preemptible = 2
  }

  # ── Step 1: create all per-SV VCFs (single task) ──────────────────────────
  call PrepVcfs {
    input:
      geno_tsv                    = geno_tsv,
      pheno_tsv                   = pheno_tsv,
      relatedness_flagged_samples = relatedness_flagged_samples,
      script                      = prep_vcfs_script,
      min_genotyped_samples       = min_genotyped_samples,
      sv_subset                   = sv_subset,
      docker                      = docker,
      memory_gb                   = prep_memory_gb,
      cpu                         = prep_cpu,
      disk_gb                     = prep_disk_gb,
      preemptible                 = prep_preemptible,
  }

  # ── Step 2: plink2 PheWAS + plotting, one shard per SV ────────────────────
  scatter (vcf_file in PrepVcfs.vcf_files) {
    call RunPhewasPerSv {
      input:
        vcf_file               = vcf_file,
        pheno_tsv              = pheno_tsv,
        covars_tsv             = covars_tsv,
        sex_file               = sex_file,
        phecode_definitions    = phecode_definitions,
        format_results_script  = format_results_script,
        plot_phewas_script     = plot_phewas_script,
        covar_names            = covar_names,
        plink_threads          = plink_threads,
        label_only_bonferroni = label_only_bonferroni,
        manhattan_label_count = manhattan_label_count,
        docker                 = docker,
        memory_gb              = phewas_memory_gb,
        cpu                    = phewas_cpu,
        disk_gb                = phewas_disk_gb,
        preemptible            = phewas_preemptible,
    }
  }

  # ── Step 3: collect all per-SV plots and results into one archive ────────────
  call GatherOutputs {
    input:
      results_tsvs    = flatten(RunPhewasPerSv.results_tsv),
      manhattan_plots = flatten(RunPhewasPerSv.manhattan_plot),
      docker          = docker,
  }

  output {
    File        sv_list        = PrepVcfs.sv_list
    Array[File] phewas_results = flatten(RunPhewasPerSv.results_tsv)
    Array[File] phewas_plots   = flatten(RunPhewasPerSv.manhattan_plot)
    File        results_tar    = GatherOutputs.results_tar
    File        plots_tar      = GatherOutputs.plots_tar
  }
}


# ── Task: PrepVcfs ────────────────────────────────────────────────────────────
task PrepVcfs {
  input {
    File   geno_tsv
    File   pheno_tsv
    File   relatedness_flagged_samples
    File   script
    Int    min_genotyped_samples
    String sv_subset
    String docker
    Int    memory_gb
    Int    cpu
    Int    disk_gb
    Int    preemptible
  }

  command <<<
    set -euo pipefail

    python ~{script} \
      --geno-tsv                    ~{geno_tsv} \
      --pheno-tsv                   ~{pheno_tsv} \
      --relatedness-flagged-samples ~{relatedness_flagged_samples} \
      --min-genotyped-samples       ~{min_genotyped_samples} \
      ~{if sv_subset != "" then "--sv-subset '" + sv_subset + "'" else ""} \
      --output-dir    plink_gt \
      --sv-list-out   SV_VCF_list.txt
  >>>

  output {
    Array[File] vcf_files = glob("plink_gt/*.vcf")
    File        sv_list   = "SV_VCF_list.txt"
  }

  runtime {
    docker:      docker
    memory:      "~{memory_gb} GB"
    cpu:         cpu
    disks:       "local-disk ~{disk_gb} SSD"
    preemptible: preemptible
    maxRetries:  preemptible
  }
}


# ── Task: RunPhewasPerSv ──────────────────────────────────────────────────────
task RunPhewasPerSv {
  input {
    File          vcf_file
    File          pheno_tsv
    File          covars_tsv
    File          sex_file
    File          phecode_definitions
    File          format_results_script
    File          plot_phewas_script
    Array[String] covar_names
    Int           plink_threads
    Boolean       label_only_bonferroni
    Int           manhattan_label_count
    String        docker
    Int           memory_gb
    Int           cpu
    Int           disk_gb
    Int           preemptible
  }

  # Derive SV ID from VCF filename (e.g. "chr19-48697404-DEL-10077_GT.vcf" → "chr19-48697404-DEL-10077")
  String sv_id = basename(vcf_file, "_GT.vcf")

  command <<<
    set -uo pipefail

    sv="~{sv_id}"

    mkdir -p plink_pgen plink_results_init plink_results_final

    run_phewas() {
      # ── Build plink2 pgen from VCF ─────────────────────────────────────────
      plink2 \
        --vcf        ~{vcf_file} \
        --make-pgen \
        --update-sex ~{sex_file} \
        --out        plink_pgen/"${sv}"

      # ── Firth regression against all phenotypes ────────────────────────────
      plink2 \
        --threads               ~{plink_threads} \
        --geno-counts \
        --pfile                 plink_pgen/"${sv}" \
        --out                   plink_results_init/"${sv}" \
        --pheno                 ~{pheno_tsv} \
        --covar                 ~{covars_tsv} \
        --covar-variance-standardize \
        --covar-name            ~{sep=" " covar_names} \
        --glm firth hide-covar cols=+nobs,+a1countcc,+gcountcc

      # ── Aggregate per-phenotype result files into one TSV ──────────────────
      # Each plink2 output file is named <sv>.<phenotype>.glm.firth;
      # this awk block prepends the phenotype name as a leading column.
      awk '
        BEGIN { OFS="\t" }
        FNR==1 {
          file = FILENAME
          sub(/^.*\//, "", file)    # strip directory
          sub(/^[^.]+\./, "", file) # strip "<sv>." prefix
          sub(/\.glm\.firth$/, "", file)
          pheno = file
          if (NR == 1) print "PHENOTYPE", $0
          next
        }
        { print pheno, $0 }
      ' plink_results_init/"${sv}".*.glm.firth > temp_aggregated.txt

      # Sort by p-value (find its column index dynamically)
      p_col=$(head -1 temp_aggregated.txt | tr '\t' '\n' | grep -n "^P$" | cut -d: -f1)
      {
        head -1 temp_aggregated.txt
        tail -n +2 temp_aggregated.txt | sort -gk"${p_col}"
      } > plink_results_final/"${sv}_results_raw.tsv"

      rm -f temp_aggregated.txt plink_results_init/"${sv}".*.glm.firth

      # ── Format results for PheTK ───────────────────────────────────────────
      python ~{format_results_script} \
        --plink-results       plink_results_final/"${sv}_results_raw.tsv" \
        --phecode-definitions ~{phecode_definitions} \
        --output              plink_results_final/"${sv}_phewas_results.tsv"

      # ── Manhattan plot ─────────────────────────────────────────────────────
      python ~{plot_phewas_script} \
        --phewas-results plink_results_final/"${sv}_phewas_results.tsv" \
        --output         "${sv}_manhattan.png" \
        --title          "PheWAS: ${sv}" \
        --manhattan-label-count  ~{manhattan_label_count} \
        ~{if label_only_bonferroni then "--label-only-bonferroni" else ""}
    }

    if ! run_phewas; then
      echo "ERROR: phewas failed for ${sv}" >&2
    fi
  >>>

  output {
    Array[File] results_tsv    = glob("plink_results_final/*_phewas_results.tsv")
    Array[File] manhattan_plot = glob("*_manhattan.png")
  }

  runtime {
    docker:      docker
    memory:      "~{memory_gb} GB"
    cpu:         cpu
    disks:       "local-disk ~{disk_gb} HDD"
    preemptible: preemptible
    maxRetries:  preemptible
  }
}


# ── Task: GatherOutputs ───────────────────────────────────────────────────────
task GatherOutputs {
  input {
    Array[File] results_tsvs
    Array[File] manhattan_plots
    String      docker
    Int         disk_gb     = 50
    Int         preemptible = 1
  }

  command <<<
    set -euo pipefail
    mkdir -p phewas_results phewas_plots

    results=(~{sep=" " results_tsvs})
    plots=(~{sep=" " manhattan_plots})

    if [ ${#results[@]} -gt 0 ]; then
      cp "${results[@]}" phewas_results/
    fi
    if [ ${#plots[@]} -gt 0 ]; then
      cp "${plots[@]}" phewas_plots/
    fi

    tar -czf phewas_results.tar.gz phewas_results/
    tar -czf phewas_plots.tar.gz   phewas_plots/
  >>>

  output {
    File results_tar = "phewas_results.tar.gz"
    File plots_tar   = "phewas_plots.tar.gz"
  }

  runtime {
    docker:      docker
    memory:      "4 GB"
    cpu:         1
    disks:       "local-disk ~{disk_gb} HDD"
    preemptible: preemptible
    maxRetries:  preemptible
  }
}
