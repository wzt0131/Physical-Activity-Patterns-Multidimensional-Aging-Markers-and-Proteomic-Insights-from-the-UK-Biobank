# colocalization

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(coloc)
})

exposure_dir <- "D:/Project/4.15_protein"
outcome_file <- "D:/phenoage_acceleration.tsv.gz"
output_dir   <- "D:/Project/5.2_colocalization"

protein_map_file <- "D:/Project/olink_protein_map_3k_v1.tsv"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

window_bp <- 1000000L
min_snps  <- 50L

p1_val  <- 1e-4
p2_val  <- 1e-4
p12_val <- 1e-5

outcome_default_n <- 311471

save_per_protein_detail <- TRUE

run_only_mr_hits <- FALSE
mr_summary_file  <- "D:/Project/5.2_results/phenoage_acceleration_summary_report.csv"
mr_p_col         <- "Main_Pvalue"
mr_p_threshold   <- 0.05

tmp_root <- file.path(output_dir, "tmp_tar_extract")
if (!dir.exists(tmp_root)) dir.create(tmp_root, recursive = TRUE)

chr_mapper <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x <- toupper(x)
  
  map <- c(
    "1"=1,"2"=2,"3"=3,"4"=4,"5"=5,"6"=6,"7"=7,"8"=8,"9"=9,"10"=10,"11"=11,
    "12"=12,"13"=13,"14"=14,"15"=15,"16"=16,"17"=17,"18"=18,"19"=19,"20"=20,
    "21"=21,"22"=22,"X"=23,"Y"=24,"MT"=25,"M"=25
  )
  
  out <- unname(map[x])
  suppressWarnings(out_num <- as.numeric(x))
  out[is.na(out)] <- out_num[is.na(out)]
  out
}

safe_upper <- function(x) {
  x <- as.character(x)
  x <- toupper(trimws(x))
  x
}

pick_first_col <- function(dt_names, candidates) {
  idx <- match(tolower(candidates), tolower(dt_names))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NA_character_)
  dt_names[idx[1]]
}

make_variant_id <- function(chr, pos, oa, ea) {
  paste(chr, as.integer(pos), oa, ea, sep = "_")
}

parse_protein_name <- function(tar_path) {
  nm <- basename(tar_path)
  nm <- sub("\\.tar$", "", nm, ignore.case = TRUE)
  parts <- strsplit(nm, "_")[[1]]
  
  gene_symbol <- if (length(parts) >= 1) parts[1] else NA_character_
  uniprot_id  <- if (length(parts) >= 2) parts[2] else NA_character_
  
  list(
    protein_name = nm,
    gene_symbol  = gene_symbol,
    uniprot_id   = uniprot_id
  )
}

coloc_support_label <- function(pp_h3, pp_h4) {
  if (is.na(pp_h4)) return("not_available")
  if (pp_h4 >= 0.80) return("strong_H4_shared_signal")
  if (pp_h4 >= 0.50) return("moderate_H4_shared_signal")
  if (!is.na(pp_h3) && pp_h3 >= 0.80) return("H3_distinct_signals_not_colocalized")
  return("weak_or_no_H4")
}

cat("Reading annotation file...\n")
protein_map <- fread(protein_map_file)

needed_map_cols <- c("UKBPPP_ProteinID", "HGNC.symbol", "UniProt", "chr", "gene_start", "gene_end")
missing_map_cols <- setdiff(needed_map_cols, names(protein_map))
if (length(missing_map_cols) > 0) {
  stop("Annotation file missing columns: ", paste(missing_map_cols, collapse = ", "))
}

gene_positions <- copy(protein_map)[
  , .(UKBPPP_ProteinID, HGNC.symbol, UniProt, chr, gene_start, gene_end)
]

gene_positions[, chr := chr_mapper(chr)]
gene_positions[, gene_start := as.numeric(gene_start)]
gene_positions[, gene_end   := as.numeric(gene_end)]

gene_positions <- gene_positions[
  !is.na(chr) & !is.na(gene_start) & !is.na(gene_end)
]

cat("Annotation file read completed, available protein annotations:", nrow(gene_positions), "\n")

read_outcome_data <- function(outcome_file) {
  cat("Reading outcome file header...\n")
  header_dt <- fread(outcome_file, nrows = 0)
  nm <- names(header_dt)
  
  col_chr  <- pick_first_col(nm, c("CHR","chr","chrom","chromosome","#chrom"))
  col_pos  <- pick_first_col(nm, c("BP","bp","GENPOS","pos","POS","base_pair_location","position"))
  col_ea   <- pick_first_col(nm, c("A1","effect_allele","EA","ea","ALLELE1","alt"))
  col_oa   <- pick_first_col(nm, c("A2","other_allele","NEA","nea","OA","oa","ALLELE0","ref"))
  col_beta <- pick_first_col(nm, c("BETA","beta","effect"))
  col_se   <- pick_first_col(nm, c("SE","se","standard_error","sebeta"))
  col_p    <- pick_first_col(nm, c("P","p","PVALUE","pval","p_value"))
  col_n    <- pick_first_col(nm, c("N","n","samplesize","sample_size"))
  col_maf  <- pick_first_col(nm, c("MAF","maf"))
  col_eaf  <- pick_first_col(nm, c("EAF","eaf","effect_allele_frequency","A1FREQ","AF","af"))
  col_rsid <- pick_first_col(nm, c("SNP","snp","rsid","rsids","ID","markername","variant_id"))
  
  required <- c(col_chr, col_pos, col_ea, col_oa, col_beta, col_se)
  if (any(is.na(required))) {
    stop(
      paste0(
        "Outcome file missing required columns. Identified results:",
        "\nchr=", col_chr,
        "\npos=", col_pos,
        "\nea=", col_ea,
        "\noa=", col_oa,
        "\nbeta=", col_beta,
        "\nse=", col_se,
        "\np=", col_p,
        "\nn=", col_n,
        "\nmaf=", col_maf,
        "\neaf=", col_eaf,
        "\nrsid=", col_rsid
      )
    )
  }
  
  select_cols <- unique(c(
    col_chr, col_pos, col_ea, col_oa, col_beta, col_se,
    col_p, col_n, col_maf, col_eaf, col_rsid
  ))
  select_cols <- select_cols[!is.na(select_cols)]
  
  cat("Reading outcome file required columns...\n")
  raw <- fread(outcome_file, select = select_cols, showProgress = TRUE)
  
  outcome_dt <- data.table(
    chr  = raw[[col_chr]],
    pos  = raw[[col_pos]],
    ea   = raw[[col_ea]],
    oa   = raw[[col_oa]],
    beta = raw[[col_beta]],
    se   = raw[[col_se]]
  )
  
  if (!is.na(col_p)) {
    outcome_dt[, p := raw[[col_p]]]
  } else {
    outcome_dt[, p := NA_real_]
  }
  
  if (!is.na(col_n)) {
    outcome_dt[, n := raw[[col_n]]]
  } else {
    outcome_dt[, n := outcome_default_n]
  }
  
  if (!is.na(col_maf)) {
    outcome_dt[, maf := raw[[col_maf]]]
  } else if (!is.na(col_eaf)) {
    outcome_dt[, maf := pmin(as.numeric(raw[[col_eaf]]), 1 - as.numeric(raw[[col_eaf]]))]
  } else {
    stop("Outcome file has neither MAF nor EAF/effect_allele_frequency, coloc cannot construct MAF.")
  }
  
  if (!is.na(col_rsid)) {
    outcome_dt[, rsid := as.character(raw[[col_rsid]])]
  } else {
    outcome_dt[, rsid := NA_character_]
  }
  
  outcome_dt[, chr  := chr_mapper(chr)]
  outcome_dt[, pos  := as.numeric(pos)]
  outcome_dt[, ea   := safe_upper(ea)]
  outcome_dt[, oa   := safe_upper(oa)]
  outcome_dt[, beta := as.numeric(beta)]
  outcome_dt[, se   := as.numeric(se)]
  outcome_dt[, p    := as.numeric(p)]
  outcome_dt[, n    := as.numeric(n)]
  outcome_dt[, maf  := as.numeric(maf)]
  
  outcome_dt <- outcome_dt[
    !is.na(chr) & !is.na(pos) &
      !is.na(ea) & !is.na(oa) &
      ea != "" & oa != "" &
      nchar(ea) == 1 & nchar(oa) == 1 &
      !is.na(beta) & !is.na(se) & se > 0 &
      !is.na(n) & n > 0 &
      !is.na(maf) & maf > 0 & maf <= 0.5
  ]
  
  outcome_dt[, variant_id := make_variant_id(chr, pos, oa, ea)]
  
  setorderv(outcome_dt, c("variant_id", "p"), c(1, 1), na.last = TRUE)
  outcome_dt <- outcome_dt[!duplicated(variant_id)]
  
  cat("Outcome file read completed, retained SNPs:", nrow(outcome_dt), "\n")
  outcome_dt
}

outcome_dt <- read_outcome_data(outcome_file)

read_and_prepare_exposure <- function(tar_file) {
  info <- parse_protein_name(tar_file)
  protein_name <- info$protein_name
  gene_symbol  <- info$gene_symbol
  uniprot_id   <- info$uniprot_id
  
  cat("\n====================================================\n")
  cat("Processing protein:", protein_name, "\n")
  cat("====================================================\n")
  
  gene_info <- gene_positions[UniProt == uniprot_id]
  if (nrow(gene_info) == 0 && !is.na(gene_symbol)) {
    gene_info <- gene_positions[HGNC.symbol == gene_symbol]
  }
  if (nrow(gene_info) == 0) {
    stop("Cannot find UniProt/HGNC mapping for this protein in annotation file")
  }
  
  gene_chr   <- gene_info$chr[1]
  gene_start <- gene_info$gene_start[1]
  gene_end   <- gene_info$gene_end[1]
  gene_hgnc  <- gene_info$HGNC.symbol[1]
  
  cis_start <- max(1, as.integer(gene_start - window_bp))
  cis_end   <- as.integer(gene_end + window_bp)
  
  cat("Gene:", gene_hgnc,
      " | UniProt:", uniprot_id,
      " | chr:", gene_chr,
      " | region:", cis_start, "-", cis_end, "\n")
  
  extract_dir <- file.path(tmp_root, protein_name)
  if (dir.exists(extract_dir)) unlink(extract_dir, recursive = TRUE, force = TRUE)
  dir.create(extract_dir, recursive = TRUE)
  
  untar(tar_file, exdir = extract_dir)
  
  all_files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
  gz_files  <- all_files[grepl("\\.gz$", all_files, ignore.case = TRUE)]
  
  if (length(gz_files) == 0) {
    stop("No .gz files found in tar")
  }
  
  cat("Number of chunk files:", length(gz_files), "\n")
  
  dt_list <- lapply(gz_files, function(f) {
    tryCatch(
      fread(f, showProgress = FALSE),
      error = function(e) NULL
    )
  })
  dt_list <- Filter(Negate(is.null), dt_list)
  
  if (length(dt_list) == 0) {
    stop("All .gz chunks failed to read")
  }
  
  exp_raw <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  cat("Exposure raw total rows:", nrow(exp_raw), "\n")
  
  nm <- names(exp_raw)
  
  col_chr  <- pick_first_col(nm, c("CHROM","CHR","chr","chrom","chromosome"))
  col_pos  <- pick_first_col(nm, c("GENPOS","BP","bp","POS","pos","base_pair_location","position"))
  col_ea   <- pick_first_col(nm, c("ALLELE1","A1","a1","effect_allele","EA","ea","alt"))
  col_oa   <- pick_first_col(nm, c("ALLELE0","A2","a2","other_allele","OA","oa","ref"))
  col_eaf  <- pick_first_col(nm, c("A1FREQ","EAF","eaf","effect_allele_frequency","AF","af"))
  col_beta <- pick_first_col(nm, c("BETA","beta","effect"))
  col_se   <- pick_first_col(nm, c("SE","se","standard_error"))
  col_n    <- pick_first_col(nm, c("N","n","samplesize","sample_size"))
  col_p    <- pick_first_col(nm, c("P","p","PVALUE","pval","p_value"))
  col_log10p <- pick_first_col(nm, c("LOG10P","log10p","LP"))
  col_rsid <- pick_first_col(nm, c("SNP","snp","rsid","RSID","ID","id","variant"))
  
  required_found <- c(col_chr, col_pos, col_ea, col_oa, col_eaf, col_beta, col_se, col_n)
  if (any(is.na(required_found))) {
    stop(
      paste0(
        "Exposure missing required columns. Identified results:",
        "\nchr=", col_chr,
        "\npos=", col_pos,
        "\nea=", col_ea,
        "\noa=", col_oa,
        "\neaf=", col_eaf,
        "\nbeta=", col_beta,
        "\nse=", col_se,
        "\nn=", col_n,
        "\np=", col_p,
        "\nlog10p=", col_log10p
      )
    )
  }
  
  exp_dt <- data.table(
    chr  = exp_raw[[col_chr]],
    pos  = exp_raw[[col_pos]],
    ea   = exp_raw[[col_ea]],
    oa   = exp_raw[[col_oa]],
    eaf  = exp_raw[[col_eaf]],
    beta = exp_raw[[col_beta]],
    se   = exp_raw[[col_se]],
    n    = exp_raw[[col_n]]
  )
  
  if (!is.na(col_p)) {
    exp_dt[, p := exp_raw[[col_p]]]
  } else if (!is.na(col_log10p)) {
    exp_dt[, p := 10^(-as.numeric(exp_raw[[col_log10p]]))]
  } else {
    exp_dt[, p := NA_real_]
  }
  
  if (!is.na(col_rsid)) {
    exp_dt[, rsid := as.character(exp_raw[[col_rsid]])]
  } else {
    exp_dt[, rsid := NA_character_]
  }
  
  exp_dt[, chr  := chr_mapper(chr)]
  exp_dt[, pos  := as.numeric(pos)]
  exp_dt[, ea   := safe_upper(ea)]
  exp_dt[, oa   := safe_upper(oa)]
  exp_dt[, eaf  := as.numeric(eaf)]
  exp_dt[, beta := as.numeric(beta)]
  exp_dt[, se   := as.numeric(se)]
  exp_dt[, n    := as.numeric(n)]
  exp_dt[, p    := as.numeric(p)]
  
  exp_dt <- exp_dt[
    !is.na(chr) & !is.na(pos) &
      !is.na(ea) & !is.na(oa) &
      ea != "" & oa != "" &
      nchar(ea) == 1 & nchar(oa) == 1 &
      !is.na(beta) & !is.na(se) & se > 0 &
      !is.na(n) & n > 0 &
      !is.na(eaf) & eaf > 0 & eaf < 1
  ]
  
  exp_dt <- exp_dt[
    chr == gene_chr & pos >= cis_start & pos <= cis_end
  ]
  
  if (nrow(exp_dt) == 0) {
    stop("No exposure SNPs in cis region")
  }
  
  exp_dt[, variant_id_forward := make_variant_id(chr, pos, oa, ea)]
  exp_dt[, variant_id_flip    := make_variant_id(chr, pos, ea, oa)]
  
  setorderv(exp_dt, c("variant_id_forward", "p"), c(1, 1), na.last = TRUE)
  exp_dt <- exp_dt[!duplicated(variant_id_forward)]
  
  list(
    protein_name = protein_name,
    gene_symbol  = gene_hgnc,
    uniprot_id   = uniprot_id,
    gene_chr     = gene_chr,
    cis_start    = cis_start,
    cis_end      = cis_end,
    exposure_dt  = exp_dt
  )
}

match_exposure_outcome <- function(exp_obj) {
  gene_chr  <- exp_obj$gene_chr
  cis_start <- exp_obj$cis_start
  cis_end   <- exp_obj$cis_end
  exp_dt    <- copy(exp_obj$exposure_dt)
  
  out_reg <- outcome_dt[
    chr == gene_chr & pos >= cis_start & pos <= cis_end
  ]
  
  if (nrow(out_reg) == 0) {
    stop("No outcome SNPs in cis region")
  }
  
  exp_dt2 <- copy(exp_dt)
  setnames(
    exp_dt2,
    old = c("chr","pos","ea","oa","eaf","beta","se","p","n","rsid"),
    new = c("chr.exposure","pos.exposure","ea.exposure","oa.exposure","eaf.exposure",
            "beta.exposure","se.exposure","p.exposure","n.exposure","rsid.exposure")
  )
  
  out_dt2 <- copy(out_reg)
  setnames(
    out_dt2,
    old = c("chr","pos","ea","oa","maf","beta","se","p","n","rsid"),
    new = c("chr.outcome","pos.outcome","ea.outcome","oa.outcome","maf.outcome",
            "beta.outcome","se.outcome","p.outcome","n.outcome","rsid.outcome")
  )
  
  exp_direct <- copy(exp_dt2)
  exp_direct[, match_id := variant_id_forward]
  out_match <- copy(out_dt2)
  out_match[, match_id := variant_id]
  
  direct <- merge(exp_direct, out_match, by = "match_id", all = FALSE)
  if (nrow(direct) > 0) {
    direct[, match_type := "direct"]
    direct[, variant_id := match_id]
  }
  
  exp_flip <- copy(exp_dt2)
  exp_flip[, match_id := variant_id_flip]
  
  flip <- merge(exp_flip, out_match, by = "match_id", all = FALSE)
  if (nrow(flip) > 0) {
    flip[, beta.exposure := -beta.exposure]
    flip[, eaf.exposure  := 1 - eaf.exposure]
    
    old_ea <- flip$ea.exposure
    old_oa <- flip$oa.exposure
    flip[, ea.exposure := old_oa]
    flip[, oa.exposure := old_ea]
    
    flip[, match_type := "flip"]
    flip[, variant_id := match_id]
  }
  
  matched <- rbindlist(list(direct, flip), use.names = TRUE, fill = TRUE)
  
  if (nrow(matched) == 0) {
    stop("No matching SNPs between exposure and outcome in cis region")
  }
  
  matched[, maf.exposure := pmin(eaf.exposure, 1 - eaf.exposure)]
  
  matched <- matched[
    !is.na(beta.exposure) & !is.na(se.exposure) & se.exposure > 0 &
      !is.na(beta.outcome)  & !is.na(se.outcome)  & se.outcome  > 0 &
      !is.na(maf.exposure)  & maf.exposure > 0 & maf.exposure <= 0.5 &
      !is.na(maf.outcome)   & maf.outcome  > 0 & maf.outcome  <= 0.5 &
      !is.na(n.exposure) & n.exposure > 0 &
      !is.na(n.outcome)  & n.outcome  > 0
  ]
  
  setorderv(matched, c("variant_id", "p.exposure", "p.outcome"), c(1, 1, 1), na.last = TRUE)
  matched <- matched[!duplicated(variant_id)]
  
  if (nrow(matched) < min_snps) {
    stop(paste0("Matching SNPs insufficient: ", nrow(matched), " < ", min_snps))
  }
  
  matched
}

run_coloc_for_one_protein <- function(tar_file) {
  info <- parse_protein_name(tar_file)
  protein_name <- info$protein_name
  protein_outdir <- file.path(output_dir, protein_name)
  if (!dir.exists(protein_outdir)) dir.create(protein_outdir, recursive = TRUE)
  
  res_row <- data.table(
    protein_name = protein_name,
    gene_symbol  = NA_character_,
    uniprot_id   = NA_character_,
    chr          = NA_real_,
    cis_start    = NA_real_,
    cis_end      = NA_real_,
    n_exposure_cis = NA_integer_,
    n_outcome_cis  = NA_integer_,
    n_matched      = NA_integer_,
    PP.H0 = NA_real_,
    PP.H1 = NA_real_,
    PP.H2 = NA_real_,
    PP.H3 = NA_real_,
    PP.H4 = NA_real_,
    top_snp = NA_character_,
    top_snp_pph4 = NA_real_,
    support = NA_character_,
    status = "INIT",
    error_message = NA_character_
  )
  
  tryCatch({
    exp_obj <- read_and_prepare_exposure(tar_file)
    
    res_row[, gene_symbol := exp_obj$gene_symbol]
    res_row[, uniprot_id  := exp_obj$uniprot_id]
    res_row[, chr         := exp_obj$gene_chr]
    res_row[, cis_start   := exp_obj$cis_start]
    res_row[, cis_end     := exp_obj$cis_end]
    res_row[, n_exposure_cis := nrow(exp_obj$exposure_dt)]
    
    out_reg_n <- outcome_dt[
      chr == exp_obj$gene_chr &
        pos >= exp_obj$cis_start &
        pos <= exp_obj$cis_end,
      .N
    ]
    res_row[, n_outcome_cis := out_reg_n]
    
    matched <- match_exposure_outcome(exp_obj)
    res_row[, n_matched := nrow(matched)]
    
    if (save_per_protein_detail) {
      fwrite(
        matched,
        file = file.path(protein_outdir, "matched_snps.tsv.gz"),
        sep = "\t"
      )
    }
    
    N1 <- as.numeric(stats::median(matched$n.exposure, na.rm = TRUE))
    N2 <- as.numeric(stats::median(matched$n.outcome,  na.rm = TRUE))
    
    d1 <- list(
      snp     = matched$variant_id,
      beta    = matched$beta.exposure,
      varbeta = matched$se.exposure^2,
      MAF     = matched$maf.exposure,
      N       = N1,
      type    = "quant"
    )
    
    d2 <- list(
      snp     = matched$variant_id,
      beta    = matched$beta.outcome,
      varbeta = matched$se.outcome^2,
      MAF     = matched$maf.outcome,
      N       = N2,
      type    = "quant"
    )
    
    coloc_res <- coloc.abf(
      dataset1 = d1,
      dataset2 = d2,
      p1 = p1_val,
      p2 = p2_val,
      p12 = p12_val
    )
    
    sm <- as.list(coloc_res$summary)
    res_row[, PP.H0 := as.numeric(sm[["PP.H0.abf"]])]
    res_row[, PP.H1 := as.numeric(sm[["PP.H1.abf"]])]
    res_row[, PP.H2 := as.numeric(sm[["PP.H2.abf"]])]
    res_row[, PP.H3 := as.numeric(sm[["PP.H3.abf"]])]
    res_row[, PP.H4 := as.numeric(sm[["PP.H4.abf"]])]
    
    snp_res <- as.data.table(coloc_res$results)
    if ("SNP.PP.H4" %in% names(snp_res)) {
      setorderv(snp_res, "SNP.PP.H4", -1, na.last = TRUE)
      res_row[, top_snp := snp_res$snp[1]]
      res_row[, top_snp_pph4 := snp_res$SNP.PP.H4[1]]
    }
    
    res_row[, support := coloc_support_label(PP.H3, PP.H4)]
    
    fwrite(
      as.data.table(coloc_res$summary),
      file = file.path(protein_outdir, "coloc_summary.tsv"),
      sep = "\t"
    )
    
    fwrite(
      snp_res,
      file = file.path(protein_outdir, "coloc_snp_results.tsv.gz"),
      sep = "\t"
    )
    
    gene_info_dt <- data.table(
      protein_name = exp_obj$protein_name,
      gene_symbol = exp_obj$gene_symbol,
      uniprot_id = exp_obj$uniprot_id,
      chr = exp_obj$gene_chr,
      cis_start = exp_obj$cis_start,
      cis_end = exp_obj$cis_end
    )
    fwrite(
      gene_info_dt,
      file = file.path(protein_outdir, "gene_info.tsv"),
      sep = "\t"
    )
    
    res_row[, status := "OK"]
    
    cat(
      "Analysis completed:", protein_name,
      " | matched SNPs=", nrow(matched),
      " | PP.H3=", signif(res_row$PP.H3, 4),
      " | PP.H4=", signif(res_row$PP.H4, 4),
      " | conclusion=", res_row$support,
      "\n"
    )
    
    extract_dir <- file.path(tmp_root, protein_name)
    if (dir.exists(extract_dir)) unlink(extract_dir, recursive = TRUE, force = TRUE)
    
    return(res_row)
    
  }, error = function(e) {
    res_row[, status := "ERROR"]
    res_row[, error_message := as.character(e$message)]
    cat("Failed:", protein_name, " | ", e$message, "\n")
    
    extract_dir <- file.path(tmp_root, protein_name)
    if (dir.exists(extract_dir)) unlink(extract_dir, recursive = TRUE, force = TRUE)
    
    fwrite(
      data.table(error_message = as.character(e$message)),
      file = file.path(protein_outdir, "error_log.tsv"),
      sep = "\t"
    )
    
    return(res_row)
  })
}

tar_files <- list.files(
  exposure_dir,
  pattern = "\\.tar$",
  full.names = TRUE,
  recursive = FALSE
)

cat("Number of tar files found:", length(tar_files), "\n")

if (length(tar_files) == 0) {
  stop("No .tar files found in exposure directory")
}

if (run_only_mr_hits) {
  if (!file.exists(mr_summary_file)) {
    stop("run_only_mr_hits=TRUE, but cannot find MR summary file:", mr_summary_file)
  }
  
  mr_dt <- fread(mr_summary_file)
  if (!mr_p_col %in% names(mr_dt)) {
    stop("MR summary file missing p-value column:", mr_p_col)
  }
  if (!"Protein" %in% names(mr_dt)) {
    stop("MR summary file missing Protein column.")
  }
  
  mr_hits <- mr_dt[!is.na(get(mr_p_col)) & get(mr_p_col) < mr_p_threshold, unique(Protein)]
  cat("MR significant proteins count:", length(mr_hits), "\n")
  
  tar_protein_names <- sub("\\.tar$", "", basename(tar_files), ignore.case = TRUE)
  tar_files <- tar_files[tar_protein_names %in% mr_hits]
  
  cat("Proteins to colocalize after matching with tar files:", length(tar_files), "\n")
  
  if (length(tar_files) == 0) {
    stop("No matches between MR significant proteins and tar files.")
  }
}

manifest <- data.table(
  tar_file = tar_files,
  protein_name = sub("\\.tar$", "", basename(tar_files), ignore.case = TRUE)
)
fwrite(manifest, file.path(output_dir, "tar_manifest.tsv"), sep = "\t")

all_results <- vector("list", length(tar_files))

start_time <- Sys.time()
cat("Starting batch colocalization:", format(start_time), "\n")

for (i in seq_along(tar_files)) {
  cat("\n###############################\n")
  cat("Progress:", i, "/", length(tar_files), "\n")
  cat("###############################\n")
  
  all_results[[i]] <- run_coloc_for_one_protein(tar_files[i])
  
  tmp_res <- rbindlist(all_results[!sapply(all_results, is.null)], fill = TRUE)
  fwrite(tmp_res, file.path(output_dir, "coloc_master_results_running.tsv"), sep = "\t")
}

end_time <- Sys.time()
cat("\nAll completed:", format(end_time), "\n")
cat("Total time:", difftime(end_time, start_time, units = "mins"), "minutes\n")

final_results <- rbindlist(all_results, fill = TRUE)

final_results[, status_order := ifelse(status == "OK", 0, 1)]
setorderv(final_results, c("status_order", "PP.H4"), c(1, -1), na.last = TRUE)
final_results[, status_order := NULL]

fwrite(final_results, file.path(output_dir, "coloc_master_results.tsv"), sep = "\t")

strong_hits <- final_results[status == "OK" & !is.na(PP.H4) & PP.H4 >= 0.80]
moderate_hits <- final_results[status == "OK" & !is.na(PP.H4) & PP.H4 >= 0.50]
distinct_h3_hits <- final_results[status == "OK" & !is.na(PP.H3) & PP.H3 >= 0.80 & (is.na(PP.H4) | PP.H4 < 0.50)]

fwrite(strong_hits,      file.path(output_dir, "coloc_strong_hits_PP.H4_ge_0.80.tsv"), sep = "\t")
fwrite(moderate_hits,    file.path(output_dir, "coloc_moderate_hits_PP.H4_ge_0.50.tsv"), sep = "\t")
fwrite(distinct_h3_hits, file.path(output_dir, "coloc_distinct_signals_PP.H3_ge_0.80.tsv"), sep = "\t")

analysis_summary <- data.table(
  item = c(
    "total_tar_files",
    "status_OK",
    "status_ERROR",
    "strong_H4_PP.H4_ge_0.80",
    "moderate_H4_PP.H4_ge_0.50",
    "distinct_H3_PP.H3_ge_0.80_and_PP.H4_lt_0.50",
    "output_dir"
  ),
  value = c(
    as.character(length(tar_files)),
    as.character(nrow(final_results[status == "OK"])),
    as.character(nrow(final_results[status == "ERROR"])),
    as.character(nrow(strong_hits)),
    as.character(nrow(moderate_hits)),
    as.character(nrow(distinct_h3_hits)),
    output_dir
  )
)
fwrite(analysis_summary, file.path(output_dir, "analysis_summary.tsv"), sep = "\t")

cat("\n============================================\n")
cat("Result files:\n")
cat("1) coloc_master_results.tsv\n")
cat("2) coloc_strong_hits_PP.H4_ge_0.80.tsv\n")
cat("3) coloc_moderate_hits_PP.H4_ge_0.50.tsv\n")
cat("4) coloc_distinct_signals_PP.H3_ge_0.80.tsv\n")
cat("5) analysis_summary.tsv\n")
cat("6) Per protein subdirectory: coloc_summary.tsv / coloc_snp_results.tsv.gz / matched_snps.tsv.gz / gene_info.tsv\n")
cat("============================================\n")