# MR

library(data.table)
library(dplyr)
library(stringr)
library(MungeSumstats)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(TwoSampleMR)
library(ieugwasr)
library(ggplot2)
library(ggrepel)

Sys.setenv(OPENGWAS_JWT = "")
user()

user_temp_dir <- "D:/temp_R_files"
if (!dir.exists(user_temp_dir)) {
  dir.create(user_temp_dir, recursive = TRUE)
}

Sys.setenv(TMP = user_temp_dir)
Sys.setenv(TEMP = user_temp_dir)
Sys.setenv(TMPDIR = user_temp_dir)

tempdir_path <- file.path(user_temp_dir, "R_temp")
if (!dir.exists(tempdir_path)) {
  dir.create(tempdir_path, recursive = TRUE)
}

old_tempdir <- tempdir()
if (!identical(tempdir_path, old_tempdir)) {
  tryCatch({
    .libPaths(c(tempdir_path, .libPaths()))
    assign(".TempDir", tempdir_path, envir = baseenv())
    cat("Temporary directory set to:", tempdir_path, "\n")
  }, error = function(e) {
    cat("Cannot change temporary directory, using default:", old_tempdir, "\n")
  })
}

protein_dir <- "D:/Project/4.15_protein"
output_base_dir <- "D:/Project/results"
protein_map_file <- "D:/Project/olink_protein_map_3k_v1.tsv"

outcomes <- list(
  list(
    name = "phenoage_acceleration",
    file = "D:/phenoage_acceleration.tsv.gz",
    samplesize = 311471,
    outcome_name = "phenoage_acceleration"
  )
)

if (!dir.exists(output_base_dir)) dir.create(output_base_dir, recursive = TRUE)
for (outcome in outcomes) {
  outcome_dir <- file.path(output_base_dir, outcome$name)
  if (!dir.exists(outcome_dir)) dir.create(outcome_dir, recursive = TRUE)
}

create_chromosome_mapper <- function() {
  chr_mapping <- c(
    "1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5, "6" = 6, "7" = 7, "8" = 8,
    "9" = 9, "10" = 10, "11" = 11, "12" = 12, "13" = 13, "14" = 14, "15" = 15,
    "16" = 16, "17" = 17, "18" = 18, "19" = 19, "20" = 20, "21" = 21, "22" = 22,
    "X" = 23, "Y" = 24, "MT" = 25, "M" = 25
  )
  
  function(chr_value) {
    if (is.numeric(chr_value)) return(chr_value)
    if (chr_value %in% names(chr_mapping)) return(chr_mapping[chr_value])
    suppressWarnings(as.numeric(chr_value))
  }
}

chr_mapper <- create_chromosome_mapper()

rename_first_match <- function(dt, candidates, new_name) {
  existing_names <- colnames(dt)
  found <- intersect(candidates, existing_names)
  
  if (length(found) == 0) return(dt)
  
  if (new_name %in% existing_names) {
    return(dt)
  } else {
    setnames(dt, found[1], new_name)
    return(dt)
  }
}

safe_extract_result <- function(mr_results, method_pattern, col_name, default = NA) {
  if (is.null(mr_results) || nrow(mr_results) == 0) return(default)
  matched_rows <- mr_results[grep(method_pattern, mr_results$method), , drop = FALSE]
  if (nrow(matched_rows) > 0 && col_name %in% colnames(matched_rows)) {
    return(matched_rows[[col_name]][1])
  }
  default
}

safe_try <- function(expr, default = NULL) {
  tryCatch(expr, error = function(e) default)
}

cat("Reading protein mapping file...\n")
protein_map <- fread(protein_map_file)
cat("Protein mapping file read completed, total proteins:", nrow(protein_map), "\n")

cat("Processing protein annotation file...\n")
gene_positions <- protein_map[
  , .(UKBPPP_ProteinID, HGNC.symbol, UniProt, chr, gene_start, gene_end)
][
  , chr := sapply(chr, chr_mapper)
][
  , gene_start := as.numeric(gene_start)
][
  , gene_end := as.numeric(gene_end)
][
  !is.na(chr) & !is.na(gene_start) & !is.na(gene_end) & !is.na(UniProt)
]

cat("Successfully extracted", nrow(gene_positions), "gene positions\n")

chr_stats <- gene_positions %>%
  count(chr) %>%
  arrange(chr)
cat("Chromosome distribution statistics:\n")
print(chr_stats)

tar_files <- list.files(protein_dir, pattern = "\\.tar$", full.names = TRUE)
cat("Found", length(tar_files), "protein files\n")

if (length(tar_files) == 0) {
  stop("No .tar files found in protein_dir, please check the path.")
}

preprocess_exposure_data <- function(tar_file, protein_name, gene_positions, tempdir_path) {
  cat("Preprocessing protein:", protein_name, "\n")
  
  tryCatch({
    file_parts <- str_split(protein_name, "_")[[1]]
    uniprot_id <- if (length(file_parts) >= 2) file_parts[2] else NA_character_
    
    if (!is.na(uniprot_id)) {
      cat("Extracted UniProt ID:", uniprot_id, "\n")
    } else {
      cat("Cannot extract UniProt ID from filename\n")
    }
    
    extract_dir <- tempfile(paste0(protein_name, "_extract_"), tmpdir = tempdir_path)
    dir.create(extract_dir, recursive = TRUE)
    
    untar(tar_file, exdir = extract_dir)
    
    all_files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
    gz_files <- all_files[str_detect(all_files, "\\.gz$")]
    
    if (length(gz_files) == 0) {
      cat(".gz files not found\n")
      return(NULL)
    }
    
    read_gz_file <- function(file_path) {
      tryCatch({
        fread(file_path, showProgress = FALSE)
      }, error = function(e) {
        cat("Cannot read file:", basename(file_path), "\n")
        NULL
      })
    }
    
    data_list <- lapply(gz_files, read_gz_file)
    data_list <- Filter(Negate(is.null), data_list)
    
    if (length(data_list) == 0) {
      cat("All .gz files failed to read\n")
      return(NULL)
    }
    
    combined_data <- rbindlist(data_list, use.names = TRUE, fill = TRUE)
    cat("Merge completed, total rows:", nrow(combined_data), "\n")
    
    if (!"P" %in% colnames(combined_data) && "LOG10P" %in% colnames(combined_data)) {
      combined_data[, P := 10^(-LOG10P)]
    }
    
    if (!"P" %in% colnames(combined_data)) {
      cat("P or LOG10P column not found in exposure file\n")
      return(NULL)
    }
    
    dt_sig <- combined_data[P < 5e-8]
    n_sig_snps <- nrow(dt_sig)
    cat("Significant SNPs (p < 5e-8):", n_sig_snps, "\n")
    
    if (n_sig_snps == 0) {
      cat("No SNPs satisfying p < 5e-8\n")
      return(NULL)
    }
    
    map_names <- list(
      chr  = c("CHROM", "chr", "CHR"),
      bp   = c("GENPOS", "pos", "BP", "bp"),
      a1   = c("ALLELE1", "A1", "a1", "effect_allele"),
      a2   = c("ALLELE0", "A2", "a2", "other_allele"),
      eaf  = c("A1FREQ", "EAF", "eaf"),
      beta = c("BETA", "beta", "effect"),
      se   = c("SE", "se"),
      n    = c("N", "n", "samplesize"),
      p    = c("P", "p", "pval", "PVALUE")
    )
    
    dt2 <- copy(dt_sig)
    new_names <- list()
    for (target in names(map_names)) {
      found <- intersect(map_names[[target]], colnames(dt2))
      if (length(found) > 0) new_names[[found[1]]] <- target
    }
    if (length(new_names) > 0) {
      setnames(dt2, old = names(new_names), new = unlist(new_names))
    }
    
    dt_out <- dt2 %>%
      dplyr::select(any_of(c("chr", "bp", "a1", "a2", "eaf", "beta", "se", "n", "p")))
    
    empty_cols <- names(which(sapply(dt_out, function(x) all(is.na(x) | x == ""))))
    if (length(empty_cols) > 0) dt_out[, (empty_cols) := NULL]
    
    cat("MungeSumstats formatting...\n")
    tmpf <- tempfile(fileext = ".tsv", tmpdir = tempdir_path)
    fwrite(dt_out, tmpf, sep = "\t", na = "NA", quote = FALSE)
    
    formatted_file <- tempfile(fileext = ".tsv.gz", tmpdir = tempdir_path)
    
    format_sumstats(
      path = tmpf,
      ref_genome = "GRCh38",
      dbSNP = "155",
      nThread = 4,
      indel = FALSE,
      save_path = formatted_file,
      return_data = FALSE,
      force_new = TRUE
    )
    
    formatted_data <- fread(formatted_file)
    cat("SNPs after formatting:", nrow(formatted_data), "\n")
    
    cat("Performing cis screening...\n")
    cis_screened <- FALSE
    
    gene_info <- NULL
    if (!is.na(uniprot_id)) {
      gene_info <- gene_positions %>% filter(UniProt == uniprot_id)
    }
    
    if (is.null(gene_info) || nrow(gene_info) == 0) {
      hgnc_symbol <- str_split(protein_name, "_")[[1]][1]
      cat("UniProt not matched, trying HGNC.symbol:", hgnc_symbol, "\n")
      gene_info <- gene_positions %>% filter(HGNC.symbol == hgnc_symbol)
    }
    
    if (nrow(gene_info) == 0) {
      cat("Cannot find gene position information\n")
      return(NULL)
    }
    
    gene_chr <- gene_info$chr[1]
    gene_start <- gene_info$gene_start[1]
    gene_end <- gene_info$gene_end[1]
    hgnc_symbol <- gene_info$HGNC.symbol[1]
    
    cis_start <- gene_start - 1000000
    cis_end <- gene_end + 1000000
    
    cat("Gene:", hgnc_symbol, "Chromosome:", gene_chr,
        "Start:", gene_start, "End:", gene_end, "\n")
    cat("Cis region:", cis_start, "-", cis_end, "\n")
    
    cis_data <- formatted_data %>%
      filter(CHR == gene_chr & BP >= cis_start & BP <= cis_end)
    
    n_cis_snps <- nrow(cis_data)
    cat("SNPs after cis screening:", n_cis_snps, "\n")
    
    if (n_cis_snps == 0) {
      cat("No SNPs retained after cis screening\n")
      return(NULL)
    }
    
    cis_screened <- TRUE
    
    exposure_dat <- cis_data %>%
      mutate(
        SNP = SNP,
        effect_allele.exposure = A1,
        other_allele.exposure = A2,
        eaf.exposure = FRQ,
        beta.exposure = BETA,
        se.exposure = SE,
        pval.exposure = P,
        samplesize.exposure = N,
        id.exposure = protein_name,
        exposure = protein_name
      ) %>%
      filter(!is.na(SNP) & !is.na(beta.exposure) & !is.na(se.exposure))
    
    if ("eaf.exposure" %in% colnames(exposure_dat)) {
      exposure_dat <- exposure_dat %>%
        filter(!is.na(eaf.exposure) & eaf.exposure >= 0.01 & eaf.exposure <= 0.99)
    }
    
    exposure_dat$f_stat <- (exposure_dat$beta.exposure / exposure_dat$se.exposure)^2
    
    exposure_dat_f <- exposure_dat %>% filter(!is.na(f_stat) & f_stat > 10)
    n_f_gt10 <- nrow(exposure_dat_f)
    
    cat("SNPs with F > 10:", n_f_gt10, "\n")
    
    if (n_f_gt10 == 0) {
      cat("No SNPs satisfying F > 10\n")
      return(NULL)
    }
    
    cat("Starting clumping...\n")
    exposure_dat_clumped <- clump_data(
      exposure_dat_f,
      clump_kb = 10000,
      clump_r2 = 0.001,
      clump_p1 = 1,
      pop = "EUR"
    )
    
    if (nrow(exposure_dat_clumped) == 0) {
      cat("No SNPs retained after clumping\n")
      return(NULL)
    }
    
    if (!"f_stat" %in% colnames(exposure_dat_clumped)) {
      exposure_dat_clumped$f_stat <- (exposure_dat_clumped$beta.exposure / exposure_dat_clumped$se.exposure)^2
    }
    
    final_iv_n <- nrow(exposure_dat_clumped)
    
    iv_summary <- data.table(
      Protein = protein_name,
      UniProt_ID = ifelse(is.na(uniprot_id), NA_character_, uniprot_id),
      HGNC_Symbol = hgnc_symbol,
      Initial_Significant_SNPs = n_sig_snps,
      Cis_SNPs = n_cis_snps,
      F_gt10_SNPs = n_f_gt10,
      Final_IVs = final_iv_n,
      F_Mean = mean(exposure_dat_clumped$f_stat, na.rm = TRUE),
      F_Median = median(exposure_dat_clumped$f_stat, na.rm = TRUE),
      F_Min = min(exposure_dat_clumped$f_stat, na.rm = TRUE),
      F_Max = max(exposure_dat_clumped$f_stat, na.rm = TRUE)
    )
    
    cat("Final IV count:", final_iv_n, "\n")
    
    return(list(
      exposure_dat = exposure_dat_clumped,
      protein_name = protein_name,
      uniprot_id = uniprot_id,
      hgnc_symbol = hgnc_symbol,
      cis_screened = cis_screened,
      n_snps = final_iv_n,
      iv_summary = iv_summary
    ))
    
  }, error = function(e) {
    cat("Error preprocessing protein", protein_name, ":", e$message, "\n")
    return(NULL)
  })
}

outcome_cache <- new.env(parent = emptyenv())

load_outcome_data <- function(outcome_info) {
  if (exists(outcome_info$name, envir = outcome_cache, inherits = FALSE)) {
    return(get(outcome_info$name, envir = outcome_cache))
  }
  
  cat("First-time reading outcome file into cache:", outcome_info$name, "\n")
  dt <- fread(outcome_info$file)
  
  dt <- rename_first_match(dt, c("CHR", "chr", "chrom", "chromosome", "#chrom"), "CHR")
  dt <- rename_first_match(dt, c("SNP", "snp", "rsid", "rsids", "variant_id", "ID", "markername"), "SNP")
  dt <- rename_first_match(dt, c("BP", "bp", "pos", "position", "base_pair_location", "GENPOS"), "BP")
  dt <- rename_first_match(dt, c("A1", "effect_allele", "EA", "ALLELE1", "alt"), "A1")
  dt <- rename_first_match(dt, c("A2", "other_allele", "NEA", "ALLELE0", "ref"), "A2")
  dt <- rename_first_match(dt, c("MAF", "maf"), "MAF")
  dt <- rename_first_match(dt, c("BETA", "beta", "effect"), "BETA")
  dt <- rename_first_match(dt, c("SE", "se", "sebeta", "standard_error"), "SE")
  dt <- rename_first_match(dt, c("P", "p", "pval", "p_value", "PVALUE"), "P")
  dt <- rename_first_match(dt, c("N", "n", "samplesize"), "N")
  
  required_cols <- c("CHR", "SNP", "BP", "A1", "A2", "BETA", "SE", "P")
  missing_cols <- setdiff(required_cols, colnames(dt))
  if (length(missing_cols) > 0) {
    stop("Outcome file missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  dt[, CHR := suppressWarnings(as.numeric(CHR))]
  dt[, BP := suppressWarnings(as.numeric(BP))]
  dt[, BETA := suppressWarnings(as.numeric(BETA))]
  dt[, SE := suppressWarnings(as.numeric(SE))]
  dt[, P := suppressWarnings(as.numeric(P))]
  
  if ("N" %in% colnames(dt)) {
    dt[, N := suppressWarnings(as.numeric(N))]
  } else {
    dt[, N := outcome_info$samplesize]
  }
  
  if ("MAF" %in% colnames(dt)) {
    dt[, MAF := suppressWarnings(as.numeric(MAF))]
  }
  
  assign(outcome_info$name, dt, envir = outcome_cache)
  dt
}

preprocess_outcome_data <- function(outcome_info, exposure_snps) {
  cat("Preprocessing outcome data:", outcome_info$name, "\n")
  
  tryCatch({
    outcome_data <- load_outcome_data(outcome_info)
    
    cat("Total rows in cached outcome:", nrow(outcome_data), "\n")
    cat("Matching exposure SNPs...\n")
    
    outcome_sub <- outcome_data[SNP %in% exposure_snps]
    
    if (nrow(outcome_sub) == 0) {
      cat("No exposure SNPs matched in outcome\n")
      return(NULL)
    }
    
    outcome_sub <- outcome_sub[!is.na(SE) & SE > 0]
    
    if (nrow(outcome_sub) == 0) {
      cat("No valid SE in matched SNPs\n")
      return(NULL)
    }
    
    if ("MAF" %in% colnames(outcome_sub)) {
      outcome_sub <- outcome_sub[is.na(MAF) | MAF >= 0.01]
    }
    
    if (nrow(outcome_sub) == 0) {
      cat("No SNPs retained after MAF filtering\n")
      return(NULL)
    }
    
    outcome_dat <- outcome_sub %>%
      mutate(
        id.outcome = outcome_info$name,
        outcome = outcome_info$outcome_name,
        effect_allele.outcome = A1,
        other_allele.outcome = A2,
        beta.outcome = BETA,
        se.outcome = SE,
        pval.outcome = P,
        samplesize.outcome = ifelse(!is.na(N), N, outcome_info$samplesize)
      ) %>%
      filter(!is.na(SNP) & !is.na(beta.outcome) & !is.na(se.outcome))
    
    if ("MAF" %in% colnames(outcome_sub)) {
      outcome_dat$maf.outcome <- outcome_sub$MAF
    }
    
    cat("SNPs retained after outcome matching:", nrow(outcome_dat), "\n")
    
    return(list(
      outcome_dat = outcome_dat,
      outcome_name = outcome_info$name,
      n_snps = nrow(outcome_dat)
    ))
    
  }, error = function(e) {
    cat("Error preprocessing outcome data:", e$message, "\n")
    return(NULL)
  })
}

perform_mr_analysis <- function(exposure_data, outcome_data, output_dir) {
  cat("Starting MR analysis:", exposure_data$protein_name, "->", outcome_data$outcome_name, "\n")
  
  tryCatch({
    protein_output_dir <- file.path(output_dir, exposure_data$protein_name)
    if (!dir.exists(protein_output_dir)) {
      dir.create(protein_output_dir, recursive = TRUE)
    }
    
    iv_summary_to_save <- copy(exposure_data$iv_summary)
    iv_summary_to_save[, Outcome := outcome_data$outcome_name]
    fwrite(iv_summary_to_save, file.path(protein_output_dir, "iv_strength_summary.csv"))
    
    exposure_dat <- exposure_data$exposure_dat %>%
      filter(SNP %in% outcome_data$outcome_dat$SNP)
    
    if (nrow(exposure_dat) == 0) {
      cat("No matched SNPs for MR analysis\n")
      return(FALSE)
    }
    
    outcome_dat <- outcome_data$outcome_dat %>%
      filter(SNP %in% exposure_dat$SNP)
    
    if (nrow(outcome_dat) == 0) {
      cat("No matched SNPs in outcome data\n")
      return(FALSE)
    }
    
    cat("Matched SNPs count:", nrow(exposure_dat), "\n")
    
    required_exposure_cols <- c(
      "SNP", "effect_allele.exposure", "other_allele.exposure",
      "beta.exposure", "se.exposure", "exposure", "id.exposure"
    )
    missing_exposure <- setdiff(required_exposure_cols, colnames(exposure_dat))
    if (length(missing_exposure) > 0) {
      cat("Exposure data missing columns:", paste(missing_exposure, collapse = ", "), "\n")
      for (col in missing_exposure) {
        if (col %in% c("exposure", "id.exposure")) {
          exposure_dat[[col]] <- exposure_data$protein_name
        } else {
          exposure_dat[[col]] <- NA
        }
      }
    }
    
    cat("Harmonising data...\n")
    harmonised_dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat,
      action = 2
    )
    
    if (nrow(harmonised_dat) == 0) {
      cat("Harmonisation failed\n")
      return(FALSE)
    }
    
    n_harmonised_before_steiger <- nrow(harmonised_dat)
    cat("SNPs after harmonisation:", n_harmonised_before_steiger, "\n")
    
    cat("Performing Steiger filtering...\n")
    
    steiger_filtered_dat <- harmonised_dat
    steiger_status <- "Skipped (no eaf)"
    steiger_removed <- 0
    
    if ("eaf.exposure" %in% colnames(harmonised_dat) && "eaf.outcome" %in% colnames(harmonised_dat)) {
      steiger_res <- tryCatch({
        directionality_test(harmonised_dat)
      }, error = function(e) {
        cat("Steiger test error:", e$message, "\n")
        NULL
      })
      
      if (!is.null(steiger_res) && nrow(steiger_res) > 0) {
        correct_dir <- steiger_res$correct_causal_direction[1]
        
        if (isTRUE(correct_dir)) {
          steiger_status <- "Passed"
          steiger_filtered_dat <- harmonised_dat
          steiger_removed <- 0
        } else {
          temp_dat <- harmonised_dat %>% filter(mr_keep)
          steiger_removed <- nrow(harmonised_dat) - nrow(temp_dat)
          
          if (nrow(temp_dat) > 0) {
            steiger_filtered_dat <- temp_dat
            steiger_status <- ifelse(
              steiger_removed > 0,
              paste0("Filtered (removed ", steiger_removed, " SNPs)"),
              "No SNPs removed"
            )
          } else {
            steiger_filtered_dat <- harmonised_dat
            steiger_status <- "Failed but no filtered set available"
            steiger_removed <- 0
          }
        }
      } else {
        steiger_status <- "Skipped (null result)"
        steiger_removed <- 0
      }
    }
    
    harmonised_dat <- steiger_filtered_dat
    
    steiger_summary <- data.table(
      Protein = exposure_data$protein_name,
      Outcome = outcome_data$outcome_name,
      Initial_Harmonised_SNPs = n_harmonised_before_steiger,
      After_Steiger_SNPs = nrow(harmonised_dat),
      Steiger_Removed = steiger_removed,
      Steiger_Status = steiger_status
    )
    fwrite(steiger_summary, file.path(protein_output_dir, "steiger_filter_summary.csv"))
    
    if (nrow(harmonised_dat) == 0) {
      cat("No SNPs retained after Steiger filtering\n")
      return(FALSE)
    }
    
    cat("SNPs after Steiger:", nrow(harmonised_dat), "\n")
    
    cat("Starting MR main analysis...\n")
    
    if (nrow(harmonised_dat) == 1) {
      cat("Single SNP: using Wald ratio\n")
      
      mr_res <- suppressWarnings(
        mr_wald_ratio(
          b_exp = harmonised_dat$beta.exposure,
          b_out = harmonised_dat$beta.outcome,
          se_exp = harmonised_dat$se.exposure,
          se_out = harmonised_dat$se.outcome
        )
      )
      
      mr_results <- data.frame(
        id.exposure = harmonised_dat$id.exposure[1],
        id.outcome = harmonised_dat$id.outcome[1],
        outcome = unique(harmonised_dat$outcome),
        exposure = unique(harmonised_dat$exposure),
        method = "Wald ratio",
        nsnp = 1,
        b = mr_res$b,
        se = mr_res$se,
        pval = mr_res$pval,
        lci95 = mr_res$b - 1.96 * mr_res$se,
        uci95 = mr_res$b + 1.96 * mr_res$se,
        stringsAsFactors = FALSE
      )
      
    } else {
      mr_results <- mr(
        harmonised_dat,
        method_list = c(
          "mr_ivw",
          "mr_egger_regression",
          "mr_weighted_median",
          "mr_weighted_mode",
          "mr_simple_mode"
        )
      )
      
      if (!is.null(mr_results) && nrow(mr_results) > 0) {
        mr_results$lci95 <- mr_results$b - 1.96 * mr_results$se
        mr_results$uci95 <- mr_results$b + 1.96 * mr_results$se
      }
    }
    
    if (is.null(mr_results) || nrow(mr_results) == 0) {
      cat("MR main analysis results empty\n")
      return(FALSE)
    }
    
    cat("Sensitivity analysis...\n")
    valid_snps <- if ("mr_keep" %in% colnames(harmonised_dat)) {
      sum(harmonised_dat$mr_keep, na.rm = TRUE)
    } else {
      nrow(harmonised_dat)
    }
    
    heterogeneity <- safe_try(
      if (valid_snps >= 2) mr_heterogeneity(harmonised_dat) else NULL,
      default = NULL
    )
    
    pleiotropy <- safe_try(
      if (valid_snps >= 3) mr_pleiotropy_test(harmonised_dat) else NULL,
      default = NULL
    )
    
    leaveoneout <- safe_try(
      if (valid_snps >= 2) mr_leaveoneout(harmonised_dat) else NULL,
      default = NULL
    )
    
    res_single <- safe_try(
      if (valid_snps >= 2) mr_singlesnp(harmonised_dat) else NULL,
      default = NULL
    )
    
    cat("Saving results...\n")
    fwrite(as.data.table(mr_results), file.path(protein_output_dir, "mr_main_results.csv"))
    fwrite(as.data.table(harmonised_dat), file.path(protein_output_dir, "harmonised_data.csv"))
    
    if (!is.null(heterogeneity) && nrow(heterogeneity) > 0) {
      fwrite(as.data.table(heterogeneity), file.path(protein_output_dir, "heterogeneity_test.csv"))
    } else {
      fwrite(data.table(Note = "Insufficient SNPs for heterogeneity analysis"),
             file.path(protein_output_dir, "heterogeneity_test.csv"))
    }
    
    if (!is.null(pleiotropy) && nrow(pleiotropy) > 0) {
      fwrite(as.data.table(pleiotropy), file.path(protein_output_dir, "pleiotropy_test.csv"))
    } else {
      fwrite(data.table(Note = "Insufficient SNPs for pleiotropy analysis"),
             file.path(protein_output_dir, "pleiotropy_test.csv"))
    }
    
    if (!is.null(leaveoneout) && nrow(leaveoneout) > 0) {
      fwrite(as.data.table(leaveoneout), file.path(protein_output_dir, "leaveoneout_analysis.csv"))
    } else {
      fwrite(data.table(Note = "Insufficient SNPs for leave-one-out analysis"),
             file.path(protein_output_dir, "leaveoneout_analysis.csv"))
    }
    
    cat("Generating single-protein plots...\n")
    
    if (valid_snps >= 2) {
      p_scatter <- safe_try(mr_scatter_plot(mr_results, harmonised_dat), default = NULL)
      if (!is.null(p_scatter) && length(p_scatter) > 0) {
        ggsave(
          file.path(protein_output_dir, "scatter_plot.png"),
          p_scatter[[1]], width = 8, height = 6, dpi = 300
        )
      }
      
      if (!is.null(res_single) && nrow(res_single) > 0) {
        p_funnel <- safe_try(mr_funnel_plot(res_single), default = NULL)
        if (!is.null(p_funnel) && length(p_funnel) > 0) {
          ggsave(
            file.path(protein_output_dir, "funnel_plot.png"),
            p_funnel[[1]], width = 8, height = 6, dpi = 300
          )
        }
        
        p_forest <- safe_try(mr_forest_plot(res_single), default = NULL)
        if (!is.null(p_forest) && length(p_forest) > 0) {
          ggsave(
            file.path(protein_output_dir, "forest_plot.png"),
            p_forest[[1]], width = 10, height = 8, dpi = 300
          )
        }
      }
      
      if (!is.null(leaveoneout) && nrow(leaveoneout) > 0) {
        p_loo <- safe_try(mr_leaveoneout_plot(leaveoneout), default = NULL)
        if (!is.null(p_loo) && length(p_loo) > 0) {
          ggsave(
            file.path(protein_output_dir, "leaveoneout_plot.png"),
            p_loo[[1]], width = 10, height = 8, dpi = 300
          )
        }
      }
    }
    
    cat("Generating single-protein summary report...\n")
    
    summary_results <- data.table(
      Protein = exposure_data$protein_name,
      UniProt_ID = ifelse(is.na(exposure_data$uniprot_id), NA_character_, exposure_data$uniprot_id),
      HGNC_Symbol = exposure_data$hgnc_symbol,
      Outcome = outcome_data$outcome_name,
      
      Initial_Significant_SNPs = exposure_data$iv_summary$Initial_Significant_SNPs[1],
      Cis_SNPs = exposure_data$iv_summary$Cis_SNPs[1],
      F_gt10_SNPs = exposure_data$iv_summary$F_gt10_SNPs[1],
      Final_IVs = exposure_data$iv_summary$Final_IVs[1],
      F_Mean = exposure_data$iv_summary$F_Mean[1],
      F_Median = exposure_data$iv_summary$F_Median[1],
      F_Min = exposure_data$iv_summary$F_Min[1],
      F_Max = exposure_data$iv_summary$F_Max[1],
      
      NSnps_After_Harmonise = nrow(harmonised_dat),
      Valid_Snps = valid_snps,
      Steiger_Status = steiger_status,
      Steiger_Removed = steiger_removed,
      
      Main_beta = safe_extract_result(mr_results, "Inverse variance weighted|Wald ratio", "b", NA),
      Main_se = safe_extract_result(mr_results, "Inverse variance weighted|Wald ratio", "se", NA),
      Main_lci95 = safe_extract_result(mr_results, "Inverse variance weighted|Wald ratio", "lci95", NA),
      Main_uci95 = safe_extract_result(mr_results, "Inverse variance weighted|Wald ratio", "uci95", NA),
      Main_Pvalue = safe_extract_result(mr_results, "Inverse variance weighted|Wald ratio", "pval", NA),
      
      Egger_beta = safe_extract_result(mr_results, "Egger", "b", NA),
      Egger_Pvalue = safe_extract_result(mr_results, "Egger", "pval", NA),
      
      Method_Used = ifelse(
        nrow(harmonised_dat) == 1,
        "Wald ratio",
        "Inverse variance weighted"
      ),
      Cis_Screened = ifelse(exposure_data$cis_screened, "Yes", "No"),
      Analysis_Status = ifelse(valid_snps > 0, "Completed", "Incomplete"),
      Notes = ifelse(
        nrow(harmonised_dat) == 1,
        "Single SNP analysis - interpret with caution",
        "Continuous outcome: interpret beta change in phenoage acceleration"
      )
    )
    
    fwrite(summary_results, file.path(protein_output_dir, "summary_results.csv"))
    
    cat("MR analysis completed\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("MR analysis error:", e$message, "\n")
    return(FALSE)
  })
}

cat("\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("Starting MR analysis workflow (UKB pQTL exposure x phenoage_acceleration outcome)\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")

overall_stats <- data.table(
  Protein = character(),
  Preprocess_Status = character(),
  MR_Analysis_Status = character(),
  Final_IVs = integer(),
  F_Mean = numeric(),
  F_Median = numeric(),
  F_Min = numeric(),
  F_Max = numeric(),
  MR_SNPs = integer(),
  Steiger_Status = character(),
  Steiger_Removed = integer()
)

successful_proteins <- 0
failed_proteins <- 0
skipped_proteins <- 0

for (outcome in outcomes) {
  invisible(load_outcome_data(outcome))
}

cat("\nStarting per-protein processing...\n")
for (i in seq_along(tar_files)) {
  tar_file <- tar_files[i]
  protein_name <- tools::file_path_sans_ext(basename(tar_file))
  
  cat("\n")
  cat(paste0(rep("=", 70), collapse = ""), "\n")
  cat("Processing protein", i, "/", length(tar_files), ":", protein_name, "\n")
  cat(paste0(rep("=", 70), collapse = ""), "\n")
  
  cat("Preprocessing protein exposure data...\n")
  exposure_data <- preprocess_exposure_data(tar_file, protein_name, gene_positions, tempdir_path)
  
  if (is.null(exposure_data)) {
    cat("Protein", protein_name, "preprocessing failed, skipping\n")
    
    overall_stats <- rbind(
      overall_stats,
      data.table(
        Protein = protein_name,
        Preprocess_Status = "Failed",
        MR_Analysis_Status = "Skipped",
        Final_IVs = 0,
        F_Mean = NA_real_,
        F_Median = NA_real_,
        F_Min = NA_real_,
        F_Max = NA_real_,
        MR_SNPs = 0,
        Steiger_Status = "Skipped",
        Steiger_Removed = 0
      ),
      fill = TRUE
    )
    
    skipped_proteins <- skipped_proteins + 1
    next
  }
  
  cat("Exposure preprocessing successful, final IV count:", exposure_data$n_snps, "\n")
  
  protein_stats <- data.table(
    Protein = protein_name,
    Preprocess_Status = "Success",
    MR_Analysis_Status = NA_character_,
    Final_IVs = exposure_data$iv_summary$Final_IVs[1],
    F_Mean = exposure_data$iv_summary$F_Mean[1],
    F_Median = exposure_data$iv_summary$F_Median[1],
    F_Min = exposure_data$iv_summary$F_Min[1],
    F_Max = exposure_data$iv_summary$F_Max[1],
    MR_SNPs = 0,
    Steiger_Status = NA_character_,
    Steiger_Removed = 0
  )
  
  for (outcome in outcomes) {
    cat("\n")
    cat(paste0(rep("-", 55), collapse = ""), "\n")
    cat("Analyzing outcome:", outcome$name, "\n")
    cat(paste0(rep("-", 55), collapse = ""), "\n")
    
    exposure_snps <- exposure_data$exposure_dat$SNP
    outcome_data <- preprocess_outcome_data(outcome, exposure_snps)
    
    if (is.null(outcome_data)) {
      cat("Outcome data preprocessing failed\n")
      protein_stats$MR_Analysis_Status <- "Failed"
      protein_stats$MR_SNPs <- 0
      protein_stats$Steiger_Status <- "Skipped"
      failed_proteins <- failed_proteins + 1
      next
    }
    
    cat("Outcome data preprocessing successful, matched SNPs:", outcome_data$n_snps, "\n")
    
    outcome_dir <- file.path(output_base_dir, outcome$name)
    success <- perform_mr_analysis(exposure_data, outcome_data, outcome_dir)
    
    if (success) {
      protein_stats$MR_Analysis_Status <- "Success"
      protein_stats$MR_SNPs <- outcome_data$n_snps
      
      steiger_file <- file.path(outcome_dir, exposure_data$protein_name, "steiger_filter_summary.csv")
      if (file.exists(steiger_file)) {
        steiger_info <- fread(steiger_file)
        protein_stats$Steiger_Status <- steiger_info$Steiger_Status[1]
        protein_stats$Steiger_Removed <- steiger_info$Steiger_Removed[1]
      }
      
      successful_proteins <- successful_proteins + 1
      cat("Protein", protein_name, "analysis completed\n")
    } else {
      protein_stats$MR_Analysis_Status <- "Failed"
      protein_stats$MR_SNPs <- 0
      protein_stats$Steiger_Status <- "Failed"
      failed_proteins <- failed_proteins + 1
      cat("Protein", protein_name, "analysis failed\n")
    }
  }
  
  overall_stats <- rbind(overall_stats, protein_stats, fill = TRUE)
}

fwrite(overall_stats, file.path(output_base_dir, "overall_analysis_statistics.csv"))
cat("\nOverall analysis statistics saved\n")

cat("\nAnalysis progress:\n")
cat("   Successfully analyzed proteins:", successful_proteins, "/", length(tar_files), "\n")
cat("   Failed proteins:", failed_proteins, "/", length(tar_files), "\n")
cat("   Skipped proteins:", skipped_proteins, "/", length(tar_files), "\n")

cat("\nGenerating overall summary results...\n")

for (outcome in outcomes) {
  cat("\nGenerating overall summary for outcome", outcome$name, "...\n")
  
  outcome_dir <- file.path(output_base_dir, outcome$name)
  all_summaries <- list()
  all_iv_summaries <- list()
  
  sub_dirs <- list.dirs(outcome_dir, recursive = FALSE)
  
  if (length(sub_dirs) > 0) {
    for (protein_dir_i in sub_dirs) {
      summary_file <- file.path(protein_dir_i, "summary_results.csv")
      iv_file <- file.path(protein_dir_i, "iv_strength_summary.csv")
      
      if (file.exists(summary_file)) {
        all_summaries[[basename(protein_dir_i)]] <- fread(summary_file)
      }
      if (file.exists(iv_file)) {
        all_iv_summaries[[basename(protein_dir_i)]] <- fread(iv_file)
      }
    }
  }
  
  if (length(all_summaries) > 0) {
    combined_summary <- rbindlist(all_summaries, fill = TRUE)
    
    required_cols <- c(
      "Main_beta", "Main_lci95", "Main_uci95", "Protein",
      "Valid_Snps", "Main_Pvalue", "Steiger_Status"
    )
    for (col in required_cols) {
      if (!col %in% colnames(combined_summary)) {
        combined_summary[, (col) := NA]
      }
    }
    
    valid_results <- combined_summary[
      !is.na(Main_beta) & !is.na(Main_lci95) & !is.na(Main_uci95)
    ]
    
    if (nrow(valid_results) > 0) {
      result_table <- valid_results %>%
        dplyr::select(
          Protein, UniProt_ID, HGNC_Symbol,
          Initial_Significant_SNPs, Cis_SNPs, F_gt10_SNPs, Final_IVs,
          F_Mean, F_Median, F_Min, F_Max,
          NSnps_After_Harmonise, Valid_Snps,
          Steiger_Status, Steiger_Removed,
          Main_beta, Main_se, Main_lci95, Main_uci95, Main_Pvalue,
          Egger_beta, Egger_Pvalue,
          Method_Used, Analysis_Status, Notes
        ) %>%
        arrange(Main_Pvalue)
      
      fwrite(result_table, file.path(outcome_dir, paste0("mr_results_table_", outcome$name, ".csv")))
      cat("Results table generated\n")
    } else {
      cat("Insufficient data to generate results table\n")
    }
    
    fwrite(combined_summary, file.path(outcome_dir, paste0("all_proteins_summary_", outcome$name, ".csv")))
    cat("Outcome overall summary saved\n")
  } else {
    cat("Single protein summary_results.csv not found\n")
  }
  
  if (length(all_iv_summaries) > 0) {
    combined_iv_summary <- rbindlist(all_iv_summaries, fill = TRUE)
    combined_iv_summary <- combined_iv_summary[order(Final_IVs, -F_Mean)]
    fwrite(combined_iv_summary, file.path(outcome_dir, paste0("iv_strength_summary_", outcome$name, ".csv")))
    cat("IV strength overall summary saved\n")
  } else {
    cat("Single protein iv_strength_summary.csv not found\n")
  }
}

cat("\nGenerating root directory summary report...\n")

for (outcome in outcomes) {
  outcome_dir <- file.path(output_base_dir, outcome$name)
  
  protein_dirs <- list.dirs(outcome_dir, recursive = FALSE)
  summary_list <- list()
  
  if (length(protein_dirs) > 0) {
    for (protein_dir_i in protein_dirs) {
      summary_file <- file.path(protein_dir_i, "summary_results.csv")
      if (file.exists(summary_file)) {
        summary_list[[basename(protein_dir_i)]] <- fread(summary_file)
      }
    }
  }
  
  if (length(summary_list) > 0) {
    phenoage_summary_report <- rbindlist(summary_list, fill = TRUE)
    fwrite(
      phenoage_summary_report,
      file.path(output_base_dir, paste0(outcome$name, "_summary_report.csv"))
    )
    cat(paste0(outcome$name, "_summary_report.csv saved\n"))
  } else {
    cat("No summary_results.csv files found for summarization\n")
  }
}

all_outcomes_iv_summary <- list()

for (outcome in outcomes) {
  outcome_dir <- file.path(output_base_dir, outcome$name)
  iv_file <- file.path(outcome_dir, paste0("iv_strength_summary_", outcome$name, ".csv"))
  
  if (file.exists(iv_file)) {
    tmp_iv <- fread(iv_file)
    tmp_iv$Outcome_Group <- outcome$name
    all_outcomes_iv_summary[[outcome$name]] <- tmp_iv
  }
}

if (length(all_outcomes_iv_summary) > 0) {
  cross_iv_summary <- rbindlist(all_outcomes_iv_summary, fill = TRUE)
  fwrite(cross_iv_summary, file.path(output_base_dir, "ukbpqtl_iv_strength_report.csv"))
  cat("ukbpqtl_iv_strength_report.csv saved\n")
}

cat("\nAnalysis statistics:\n")
cat("   Total protein files:", length(tar_files), "\n")
cat("   Successfully preprocessed proteins:", nrow(overall_stats[Preprocess_Status == "Success"]), "\n")
cat("   Successfully completed MR analysis proteins:", nrow(overall_stats[MR_Analysis_Status == "Success"]), "\n")

if (nrow(overall_stats[MR_Analysis_Status == "Success"]) > 0) {
  cat("   Mean final IV count:", round(mean(overall_stats[MR_Analysis_Status == "Success", Final_IVs], na.rm = TRUE), 2), "\n")
  cat("   Mean F statistic:", round(mean(overall_stats[MR_Analysis_Status == "Success", F_Mean], na.rm = TRUE), 2), "\n")
  cat("   Mean MR matched SNPs:", round(mean(overall_stats[MR_Analysis_Status == "Success", MR_SNPs], na.rm = TRUE), 2), "\n")
}

final_report <- data.table(
  Analysis_Item = c(
    "Total protein files",
    "Successfully preprocessed proteins",
    "Successfully completed MR analysis proteins",
    "Mean final IV count",
    "Mean F statistic",
    "Results save path"
  ),
  Result = c(
    as.character(length(tar_files)),
    as.character(nrow(overall_stats[Preprocess_Status == "Success"])),
    as.character(nrow(overall_stats[MR_Analysis_Status == "Success"])),
    ifelse(
      nrow(overall_stats[MR_Analysis_Status == "Success"]) > 0,
      as.character(round(mean(overall_stats[MR_Analysis_Status == "Success", Final_IVs], na.rm = TRUE), 2)),
      "0"
    ),
    ifelse(
      nrow(overall_stats[MR_Analysis_Status == "Success"]) > 0,
      as.character(round(mean(overall_stats[MR_Analysis_Status == "Success", F_Mean], na.rm = TRUE), 2)),
      "0"
    ),
    output_base_dir
  )
)

fwrite(final_report, file.path(output_base_dir, "analysis_summary.csv"))
print(final_report)

cat("\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("All analyses completed! Results saved in:", output_base_dir, "\n")
cat("Main output file descriptions:\n")
cat("   1. Per protein directory:\n")
cat("      - mr_main_results.csv\n")
cat("      - harmonised_data.csv\n")
cat("      - heterogeneity_test.csv\n")
cat("      - pleiotropy_test.csv\n")
cat("      - leaveoneout_analysis.csv\n")
cat("      - steiger_filter_summary.csv\n")
cat("      - iv_strength_summary.csv\n")
cat("      - summary_results.csv\n")
cat("   2. phenoage_acceleration directory:\n")
cat("      - all_proteins_summary_phenoage_acceleration.csv\n")
cat("      - mr_results_table_phenoage_acceleration.csv\n")
cat("      - iv_strength_summary_phenoage_acceleration.csv\n")
cat("   3. Root directory:\n")
cat("      - overall_analysis_statistics.csv\n")
cat("      - phenoage_acceleration_summary_report.csv\n")
cat("      - ukbpqtl_iv_strength_report.csv\n")
cat("      - analysis_summary.csv\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")
