library(optparse)
library(readr)
library(bigreadr)
library(bigsnpr)
library(data.table)
library(dplyr)
library(stringr)

option_list = list(
  make_option(c("--trait"), type = "character", default = NULL,
              help = "Trait ID", metavar = "character")
)
opt_parser = OptionParser(option_list = option_list)
option = parse_args(opt_parser)

if (is.null(option$trait)) {
  stop("Please provide --trait.")
}

if (utils::packageVersion("bigsnpr") < package_version("1.11.4")) {
  stop("This revised script requires bigsnpr >= 1.11.4 for LDpred2-auto.")
}

opt = list(
  LDrefpanel = "1kg",
  partitions = "0.8,0.2",
  delta = "0.001,0.01,0.1,1",
  nlambda = 30,
  lambda.min.ratio = 0.01,
  alpha = "0.7,1.0,1.4",
  p_seq = "1.0e-05,3.2e-05,1.0e-04,3.2e-04,1.0e-03,3.2e-03,1.0e-02,3.2e-02,1.0e-01,3.2e-01,1.0e+00",
  sparse = FALSE,
  kb = 500,
  Pvalthr = "5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01",
  R2 = "0.1",
  ensemble = FALSE,
  type = "logical",
  verbose = 1
)


# Example manual input (kept close to the original script): ---------
userID = "user1"
submissionID = "single_ans"
method = "LDpred2_auto"
trait = option$trait
# trait = 23
race = "EUR"
LDrefpanel = "1kg"
ensemble = FALSE

# ------------------------------------------------------------------

# Optional input parameters retained from the original script.
partitions <- opt$partitions
homedir = "/path/to/PennPRS/Files/"
type = "single-ancestry"
ld_path <- paste0(homedir)
threads = 1

trait_name = paste0(race, "_", trait)
ld_path0 <- paste0(ld_path, "LD_1kg/")
if (LDrefpanel == "1kg") {
  eval_ld_ref_path <- paste0(ld_path, "/1KGref_plinkfile/", race, "/")
  path_precalLD <- paste0(ld_path, "/LDpred2_lassosum2_corr_1kg/")
}

jobID = paste(c(trait, race, method, userID, submissionID), collapse = "_")
tempdir = "/path/to/summary_training/results_DXA/summery/"
workdir = paste0(tempdir, jobID, "/")
suppressWarnings(dir.create(workdir, recursive = TRUE))
setwd(workdir)

h2.ratio = as.numeric(str_split(opt$alpha, ",")[[1]])
p_init_seq <- as.numeric(str_split(opt$p_seq, ",")[[1]])
sp.temp = toupper(str_split(as.character(opt$sparse), ",")[[1]])
sparse.option = sp.temp == "TRUE"
if (length(sparse.option) != 1) {
  stop("LDpred2-auto requires a single sparse setting. Please set opt$sparse to either TRUE or FALSE.")
}

gwas_path <- paste0(workdir, "sumdata/")
PennPRS_finalresults_path <- paste0(workdir, "PennPRS_results/")
prsdir0 = paste0(workdir, "PRS_model_training/")
prsdir = paste0(prsdir0, method, "/")
dir.create(gwas_path, showWarnings = FALSE, recursive = TRUE)
dir.create(PennPRS_finalresults_path, showWarnings = FALSE, recursive = TRUE)
dir.create(prsdir, showWarnings = FALSE, recursive = TRUE)

source_sumstats = paste0("/path/to/summary_training/results_DXA/summery/", trait_name, ".txt")
target_sumstats = paste0(gwas_path, trait_name, ".txt")
if (!file.exists(source_sumstats)) {
  stop(paste0("Input GWAS summary file not found: ", source_sumstats))
}
file.copy(source_sumstats, target_sumstats, overwrite = TRUE)

make_h2_init_seq <- function(ldsc_h2_est, h2_ratio) {
  h2_seq <- round(ldsc_h2_est * h2_ratio, 5)
  h2_seq[h2_seq == 0] <- 1e-5
  dup_idx <- duplicated(h2_seq)
  if (any(dup_idx)) h2_seq[dup_idx] <- h2_seq[dup_idx] * 1.01
  n.inflated <- sum(h2_seq > 1)
  if (n.inflated > 0) {
    h2_seq[h2_seq > 1] <- 0.95 + seq(0, 0.01 * (n.inflated - 1), by = 0.01)
  }
  h2_seq
}

extract_chain_beta <- function(auto_fit, sparse_requested) {
  if (sparse_requested) {
    if (!is.null(auto_fit$beta_est_sparse)) return(as.numeric(auto_fit$beta_est_sparse))
  }
  as.numeric(auto_fit$beta_est)
}

standardize_sumstats <- function(sumraw) {
  colnames_upper <- toupper(colnames(sumraw))

  find_col <- function(candidates, required = TRUE) {
    hit <- match(candidates, colnames_upper)
    hit <- hit[!is.na(hit)]
    if (length(hit) > 0) return(colnames(sumraw)[hit[1]])
    if (required) {
      stop(paste0(
        "Input summary statistics are missing required column(s): ",
        paste(candidates, collapse = "/"),
        ". Available columns are: ",
        paste(colnames(sumraw), collapse = ", ")
      ))
    }
    NA_character_
  }

  chr_col <- find_col(c("CHR", "#CHROM", "CHROM"))
  snp_col <- find_col(c("SNP", "RSID", "ID", "MARKERNAME"))
  a1_col <- find_col(c("A1", "EA", "EFFECT_ALLELE", "ALLELE1"))
  a2_col <- find_col(c("A2", "NEA", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "ALLELE2"))
  beta_col <- find_col(c("BETA", "EFFECT", "ESTIMATE"))
  se_col <- find_col(c("SE", "STDERR", "SEBETA"))
  p_col <- find_col(c("P", "PVAL", "PVALUE", "P_VALUE"))
  n_col <- find_col(c("N", "N_EFF", "NEFF", "TOTAL_N", "OBS_CT"))
  af_col <- find_col(
    c("AF1", "A1FREQ", "A1_AF", "EAF", "EAF1", "ALT_FREQ", "FRQ", "MAF"),
    required = FALSE
  )

  if (is.na(af_col)) {
    stop(paste0(
      "Input summary statistics must contain an allele-frequency column such as AF1 or MAF. ",
      "Available columns are: ", paste(colnames(sumraw), collapse = ", ")
    ))
  }

  standardized <- data.frame(
    CHR = sumraw[[chr_col]],
    SNP = sumraw[[snp_col]],
    A1 = sumraw[[a1_col]],
    A2 = sumraw[[a2_col]],
    BETA = sumraw[[beta_col]],
    SE = sumraw[[se_col]],
    P = sumraw[[p_col]],
    N = sumraw[[n_col]],
    MAF = sumraw[[af_col]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  standardized
}

compute_corr_range <- function(auto_fit) {
  corr_est <- auto_fit$corr_est
  corr_est <- corr_est[is.finite(corr_est)]
  if (length(corr_est) == 0) return(NA_real_)
  diff(range(corr_est))
}

build_ld_matrix <- function(df_beta, map_ldref, path_precalLD, race, prsdir) {
  ld_tmpdir <- paste0(prsdir, "temporary_LDpred2_auto/")
  if (!dir.exists(ld_tmpdir)) dir.create(ld_tmpdir, recursive = TRUE)
  tmp_sfbm <- tempfile(tmpdir = ld_tmpdir)
  corr <- NULL
  ld <- numeric()

  for (chr in 1:22) {
    ind.chr <- which(df_beta$chr == chr)
    if (length(ind.chr) == 0) next

    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
    ind.chr3 <- ind.chr3[!is.na(ind.chr3)]
    if (length(ind.chr3) == 0) next

    corr0 <- readRDS(paste0(path_precalLD, race, "/LD_ref_chr", chr, ".rds"))[ind.chr3, ind.chr3]
    if (length(corr0) == 1) corr0 <- as(1, "sparseMatrix")

    if (is.null(corr)) {
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, tmp_sfbm, compact = TRUE)
    } else {
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }

    if (opt$verbose >= 1) {
      print(paste0("Complete calculating LD for CHR ", chr))
    }
    rm(corr0)
  }

  if (is.null(corr) || length(ld) == 0) {
    stop("Unable to build the LD correlation matrix for the matched variants.")
  }

  list(corr = corr, ld = ld, tmp_sfbm = tmp_sfbm)
}

run_ldpred2_auto_suite <- function(corr, df_beta, h2_init_seq, p_init_seq, sparse_requested, ncores) {
  coef_shrink <- 0.95
  burn_in <- 500
  num_iter <- 500
  use_MLE <- TRUE
  allow_jump_sign <- FALSE

  all_chains <- list()
  chain_records <- list()
  chain_counter <- 0L

  set.seed(2023)
  for (h2_idx in seq_along(h2_init_seq)) {
    h2_init <- h2_init_seq[h2_idx]
    if (opt$verbose >= 1) {
      print(paste0("Starting LDpred2-auto chains for h2_init = ", signif(h2_init, 4)))
    }

    auto_block <- snp_ldpred2_auto(
      corr = corr,
      df_beta = df_beta,
      h2_init = h2_init,
      vec_p_init = p_init_seq,
      burn_in = burn_in,
      num_iter = num_iter,
      sparse = sparse_requested,
      allow_jump_sign = allow_jump_sign,
      shrink_corr = coef_shrink,
      use_MLE = use_MLE,
      ncores = min(ncores, length(p_init_seq))
    )

    for (p_idx in seq_along(auto_block)) {
      chain_counter <- chain_counter + 1L
      auto_fit <- auto_block[[p_idx]]
      beta_used <- extract_chain_beta(auto_fit, sparse_requested)
      chain_records[[chain_counter]] <- data.frame(
        chain_id = chain_counter,
        h2_init_index = h2_idx,
        p_init_index = p_idx,
        h2_init = as.numeric(auto_fit$h2_init),
        p_init = as.numeric(auto_fit$p_init),
        sparse = sparse_requested,
        h2_est = as.numeric(auto_fit$h2_est),
        p_est = as.numeric(auto_fit$p_est),
        alpha_est = as.numeric(auto_fit$alpha_est),
        range_corr = compute_corr_range(auto_fit),
        n_nonzero_beta = sum(abs(beta_used) > 1e-10, na.rm = TRUE),
        has_missing_beta = anyNA(beta_used),
        stringsAsFactors = FALSE
      )
      all_chains[[chain_counter]] <- auto_fit
    }
  }

  chain_summary <- bind_rows(chain_records)
  chain_summary$valid_chain <- !chain_summary$has_missing_beta &
    is.finite(chain_summary$h2_est) &
    is.finite(chain_summary$p_est) &
    is.finite(chain_summary$alpha_est) &
    is.finite(chain_summary$range_corr) &
    (chain_summary$n_nonzero_beta > 0)

  if (!any(chain_summary$valid_chain)) {
    stop("All LDpred2-auto chains were invalid (missing effects or invalid corr_est range).")
  }

  range_threshold <- 0.95 * as.numeric(quantile(
    chain_summary$range_corr[chain_summary$valid_chain],
    probs = 0.95,
    na.rm = TRUE,
    names = FALSE
  ))
  keep_chain <- chain_summary$valid_chain & (chain_summary$range_corr > range_threshold)

  fallback_used <- FALSE
  if (!any(keep_chain)) {
    fallback_used <- TRUE
    best_chain <- which.max(ifelse(chain_summary$valid_chain, chain_summary$range_corr, -Inf))
    keep_chain[best_chain] <- TRUE
  }

  chain_summary$kept_for_final <- keep_chain
  chain_summary$selection_rule <- "range_corr > 0.95 * quantile(range_corr, 0.95)"

  list(
    auto_fits = all_chains,
    chain_summary = chain_summary,
    keep_chain = keep_chain,
    range_threshold = range_threshold,
    coef_shrink = coef_shrink,
    burn_in = burn_in,
    num_iter = num_iter,
    use_MLE = use_MLE,
    allow_jump_sign = allow_jump_sign,
    fallback_used = fallback_used
  )
}

cat(paste0("\n********************************************"))
cat(paste0("\n**** Step 0: QC for the input GWAS data ****"))
cat(paste0("\n********************************************\n"))

sumraw = bigreadr::fread2(target_sumstats)
sumraw = standardize_sumstats(sumraw)
sumraw$BETA = as.numeric(sumraw$BETA)
sumraw$SE = as.numeric(sumraw$SE)
sumraw$MAF = as.numeric(sumraw$MAF)
sumraw$P = as.numeric(sumraw$P)
sumraw$N = as.numeric(sumraw$N)

chi2_thr = 30
remaining.SNPs = which(abs(sumraw$BETA / sumraw$SE) < sqrt(chi2_thr))
if (length(remaining.SNPs) < 5) {
  stop("[Terminated] Job is terminated because less than 5 SNPs have z-score < sqrt(30), suggesting issues with the input GWAS data.")
}

n.na = sum(!complete.cases(sumraw))
if (n.na > 0) {
  sumraw = sumraw[complete.cases(sumraw), ]
  if (n.na == 1) print("* 1 SNP has missing GWAS summary-level information and is removed.")
  if (n.na > 1) print(paste0("* ", n.na, " SNPs have missing GWAS summary-level information and are removed."))
}

beta.thr = 1e3
rm.indx1 = which(abs(sumraw$BETA) > beta.thr)
if (length(rm.indx1) > 0) {
  if (length(rm.indx1) == 1) print(paste0("* 1 SNP has problematic GWAS summary statistic with abs(BETA) > ", beta.thr, " and is removed."))
  if (length(rm.indx1) > 1) print(paste0("* ", length(rm.indx1), " SNPs have problematic GWAS summary statistics with abs(BETA) > ", beta.thr, " and are removed."))
}

rm.indx2 = which((sumraw$P > 1) | (sumraw$P < 0))
if (length(rm.indx2) > 0) {
  if (length(rm.indx2) == 1) print("* 1 SNP has p-value > 1 or < 0 and is removed.")
  if (length(rm.indx2) > 1) print(paste0("* ", length(rm.indx2), " SNPs have p-value > 1 or < 0 and are removed."))
}

rm.indx3 = numeric()

chi2.thr = 1e3
rm.indx4 = which((sumraw$BETA / sumraw$SE)^2 > chi2.thr)
if (length(rm.indx4) > 0) {
  if (length(rm.indx4) == 1) print(paste0("* 1 SNP has an extremely large effect size  (z-score^2 > ", chi2.thr, ") and is removed."))
  if (length(rm.indx4) > 1) print(paste0("* ", length(rm.indx4), " SNPs have extremely large effect sizes  (z-score^2 > ", chi2.thr, ") and are removed."))
}

rm.indx5 = which(sumraw$SE == 0)
if (length(rm.indx5) > 0) {
  if (length(rm.indx5) == 1) print("* 1 SNP has SE = 0 and is removed.")
  if (length(rm.indx5) > 1) print(paste0("* ", length(rm.indx5), " SNPs have SE = 0 and are removed."))
}

rm.indx = unique(c(rm.indx1, rm.indx2, rm.indx3, rm.indx4, rm.indx5))
if (length(rm.indx) > 0) {
  sumraw = sumraw[-rm.indx, ]
  if (nrow(sumraw) == 0) {
    stop("[Terminated] 0 SNPs remaining after QC. Job terminated.\n * Please check the quality of the input GWAS summary data and make sure the columns are in correct format.")
  }
  write_delim(sumraw, path = target_sumstats, delim = "\t")
  if (length(rm.indx) == 1) print("* 1 problematic SNP removed. QC step completed.")
  if (length(rm.indx) > 1) print(paste0("* QC step completed. ", nrow(sumraw), " SNPs remaining. ", length(rm.indx), " problematic SNPs removed."))
}
if (length(rm.indx) == 0) print(paste0("* QC step completed. ", nrow(sumraw), " SNPs remaining. No SNP was removed."))

cat(paste0("\n********************************************************"))
cat(paste0("\n**** Step 1: Prepare matched data for LDpred2-auto ****"))
cat(paste0("\n********************************************************\n"))

map_ldref <- readRDS(paste0(ld_path, "map_", LDrefpanel, "_ldref_", race, ".rds"))
sumstats = sumraw[, c("CHR", "SNP", "A1", "A2", "BETA", "SE", "P", "N", "MAF")]
names(sumstats) <- c("chr", "rsid", "a1", "a0", "beta", "beta_se", "p", "n_eff", "a1_sumdata_af")

info_snp <- snp_match(sumstats, map_ldref, strand_flip = TRUE, join_by_pos = FALSE)
info_snp <- as.data.frame(info_snp)
info_snp <- info_snp[complete.cases(info_snp), ]
sd_ldref <- with(info_snp, sqrt(2 * a1_af * (1 - a1_af)))
sd_ss <- with(info_snp, sqrt(2 * a1_sumdata_af * (1 - a1_sumdata_af)))
is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
df_beta <- info_snp[!is_bad, ]
if (nrow(df_beta) == 0) {
  stop("No SNPs remained after matching and LDpred2 QC against the LD reference.")
}
write_delim(as.data.frame(df_beta), path = paste0(prsdir, trait_name, ".", method, ".matched_sumstats.txt"), delim = "\t")

ld_obj <- build_ld_matrix(df_beta, map_ldref, path_precalLD, race, prsdir)
corr <- ld_obj$corr
ld <- ld_obj$ld

cat(paste0("\n******************************************************"))
cat(paste0("\n**** Step 2: Run LD score regression for h2 init ****"))
cat(paste0("\n******************************************************\n"))

ldsc <- with(df_beta, snp_ldsc(
  ld,
  length(ld),
  chi2 = (beta / beta_se)^2,
  sample_size = n_eff,
  blocks = NULL
))
ldsc_h2_est <- abs(ldsc[["h2"]])
h2_init_seq <- make_h2_init_seq(ldsc_h2_est, h2.ratio)

cat(paste0("Heritability estimate based on LD score regression: ", signif(ldsc[["h2"]], 4), "\n"))
cat(paste0("LDpred2-auto h2 initialization values: ", paste(signif(h2_init_seq, 4), collapse = ", "), "\n"))
cat(paste0("LDpred2-auto p initialization values: ", paste(signif(p_init_seq, 4), collapse = ", "), "\n"))
cat(paste0("LDpred2-auto sparse output requested: ", sparse.option, "\n"))

cat(paste0("\n*******************************************"))
cat(paste0("\n**** Step 3: Run LDpred2-auto chains ****"))
cat(paste0("\n*******************************************\n"))

auto_run <- run_ldpred2_auto_suite(
  corr = corr,
  df_beta = df_beta,
  h2_init_seq = h2_init_seq,
  p_init_seq = p_init_seq,
  sparse_requested = sparse.option,
  ncores = 17
)

chain_summary <- auto_run$chain_summary
keep_chain <- auto_run$keep_chain

beta_matrix_list <- lapply(auto_run$auto_fits, extract_chain_beta, sparse_requested = sparse.option)
beta_matrix_list <- lapply(beta_matrix_list, function(beta_vec) {
  beta_vec[is.na(beta_vec)] <- 0
  beta_vec
})
beta_matrix <- do.call(cbind, beta_matrix_list)
if (is.null(dim(beta_matrix))) beta_matrix <- matrix(beta_matrix, ncol = 1)
colnames(beta_matrix) <- paste0(method, "_", seq_len(ncol(beta_matrix)))

ldpred2.full <- data.frame(
  CHR = df_beta$chr,
  SNP = df_beta$rsid,
  A1 = df_beta$a1,
  A2 = df_beta$a0,
  beta_matrix,
  check.names = FALSE
)
nonzero.rows <- which(rowSums(abs(beta_matrix) > 0, na.rm = TRUE) > 0)
if (length(nonzero.rows) > 0) {
  ldpred2.full <- ldpred2.full[nonzero.rows, ]
}
write_delim(ldpred2.full, path = paste0(prsdir, trait_name, ".", method, ".full.txt"), delim = "\t")

beta_keep <- beta_matrix[, keep_chain, drop = FALSE]
beta_final <- rowMeans(beta_keep)
final_prs <- data.frame(
  CHR = df_beta$chr,
  SNP = df_beta$rsid,
  A1 = df_beta$a1,
  A2 = df_beta$a0,
  BETA = beta_final
)
final_prs <- final_prs[abs(final_prs$BETA) > 1e-10, ]
if (nrow(final_prs) == 0) {
  stop("The retained LDpred2-auto chains produced an all-zero final PRS.")
}
write_delim(final_prs, path = paste0(PennPRS_finalresults_path, trait_name, ".", method, ".PRS.txt"), delim = "\t")

chain_summary_out <- chain_summary[, c(
  "chain_id", "h2_init_index", "p_init_index", "h2_init", "p_init", "sparse",
  "h2_est", "p_est", "alpha_est", "range_corr", "n_nonzero_beta",
  "has_missing_beta", "valid_chain", "kept_for_final", "selection_rule"
)]
write_delim(chain_summary_out,
            path = paste0(PennPRS_finalresults_path, trait_name, ".", method, ".chain_estimates.txt"),
            delim = "\t")

optimal_params_out <- chain_summary_out[chain_summary_out$kept_for_final, ]
write_delim(optimal_params_out,
            path = paste0(PennPRS_finalresults_path, trait_name, ".", method, ".optimal_params.txt"),
            delim = "\t")

selection_summary <- data.frame(
  method = method,
  trait = trait,
  race = race,
  ldrefpanel = LDrefpanel,
  ldsc_h2_est = as.numeric(ldsc[["h2"]]),
  h2_init_values = paste(signif(h2_init_seq, 6), collapse = ","),
  p_init_values = paste(signif(p_init_seq, 6), collapse = ","),
  sparse = sparse.option,
  burn_in = auto_run$burn_in,
  num_iter = auto_run$num_iter,
  allow_jump_sign = auto_run$allow_jump_sign,
  shrink_corr = auto_run$coef_shrink,
  use_MLE = auto_run$use_MLE,
  range_threshold = auto_run$range_threshold,
  fallback_used = auto_run$fallback_used,
  n_total_chains = nrow(chain_summary_out),
  n_kept_chains = sum(keep_chain),
  kept_chain_ids = paste(chain_summary_out$chain_id[keep_chain], collapse = ","),
  mean_h2_est_kept = mean(chain_summary_out$h2_est[keep_chain], na.rm = TRUE),
  mean_p_est_kept = mean(chain_summary_out$p_est[keep_chain], na.rm = TRUE),
  mean_alpha_est_kept = mean(chain_summary_out$alpha_est[keep_chain], na.rm = TRUE),
  stringsAsFactors = FALSE
)
write_delim(selection_summary,
            path = paste0(PennPRS_finalresults_path, trait_name, ".", method, ".selection_summary.txt"),
            delim = "\t")

kept_chain_ids <- chain_summary_out$chain_id[keep_chain]
save(
  chain_summary,
  kept_chain_ids,
  selection_summary,
  ldsc,
  ldsc_h2_est,
  file = paste0(PennPRS_finalresults_path, "step1.RData")
)

if (file.exists(paste0(ld_obj$tmp_sfbm, ".sbk"))) {
  file.remove(paste0(ld_obj$tmp_sfbm, ".sbk"))
}

cat(paste0("\n********************************************************"))
cat(paste0("\n**** LDpred2-auto training completed successfully ****"))
cat(paste0("\n********************************************************\n"))
