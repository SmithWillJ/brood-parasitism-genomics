### LFMM - one het filter per run (index 1-3)  -- alpha = 0.05
### Reorders sample info to VCF/.geno order, runs observed LFMM, GENERATES site-level
### JOINT permutations and writes all results to LFMM_out_alpha05_<filter>.
### "joint_support" = FDR candidate AND joint permutation p < alpha.

library(LEA)
library(raster)
library(tidyverse)
library(lfmm)
library(qqman)

#### Settings ####
n_perms_joint <- 5000                 # number of null permutations
fdr     <- 0.1                        # FDR threshold (Benjamini-Hochberg)
alpha   <- 0.05                       # permutation p-value threshold
predNames_analysis <- c("parasit", "lat", "bio4", "bio8", "bio13", "bio15")  # predictors tested
joint_seed <- 100000L                 # base seed for permutations (reproducibility)

#### Choose filter ####
args <- commandArgs(trailingOnly = TRUE)
idx <- as.integer(args[1])
filters <- c("0.5", "0.6", "0.75")    # heterozygosity thresholds
if (is.na(idx) || idx < 1 || idx > length(filters)) stop("Give a filter index 1-3")
het_filter <- filters[idx]

#### Paths ####
base      <- "/project/filepath/Reed_Warblers_2026"
run_dir   <- file.path(base, "filepath/populations_runs", paste0("GEA_auto_", het_filter, "_run3"))
prefix    <- "populations.snps.chromfilt.ABfilt.maxmiss0.8.recode"
geno_file <- file.path(run_dir, paste0(prefix, ".geno"))
vcf_file  <- file.path(run_dir, paste0(prefix, ".vcf"))
pos_file  <- file.path(base, "scripts", paste0("POS_FIL_EDITED_", het_filter, ".012.pos"))
out_dir   <- file.path(base, "scripts", paste0("LFMM_out_alpha05_", het_filter))
csv_file  <- file.path(base, "filepath/populations_runs/RW_LFMM_sampleinfo2_prec.csv")
bio_dir   <- file.path(base, "BioClim/")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("##### LFMM, het filter", het_filter, "| alpha =", alpha, "#####\n\n")

# Standard BH: rejects everything up to the largest rank that passes the threshold.
# Returns ORIGINAL SNP indices, ordered by p-value.
get_candidates <- function(pval, fdr) {
  ok <- which(!is.na(pval) & is.finite(pval))
  if (length(ok) == 0) return(integer())
  ord <- order(pval[ok])
  sorted_p <- pval[ok][ord]
  passing <- which(sorted_p <= fdr * seq_along(sorted_p) / length(sorted_p))
  if (length(passing) == 0) return(integer())
  ok[ord[seq_len(max(passing))]]
}

# Permutation p-value: how often the null minimum p is <= observed p
# (+1 correction avoids p = 0)
perm_p <- function(p, null_minp) {
  ns <- sort(null_minp)
  (findInterval(p, ns) + 1) / (length(ns) + 1)
}

# md5 of a serialised R object; used as the permutation cache signature
object_md5 <- function(object) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  con <- file(tmp, open = "wb")
  tryCatch(writeBin(serialize(object, NULL, version = 2), con),  # version 2 = stable across R versions
           finally = close(con))
  unname(tools::md5sum(tmp))
}

#### Genotypes + position check ####
if (!file.exists(geno_file)) {                         # convert VCF -> .geno if missing
  if (!file.exists(vcf_file)) stop("Neither geno nor vcf found for ", het_filter)
  vcf2geno(vcf_file)
}
snp_matrix <- read.geno(geno_file)                     # individuals x SNPs (0/1/2, 9 = missing)
n_snps <- ncol(snp_matrix); n_inds <- nrow(snp_matrix)
cat("SNP matrix:", n_inds, "individuals x", n_snps, "SNPs\n")

# Impute missing genotypes with the most common genotype per SNP
n_miss <- sum(snp_matrix == 9)
if (n_miss > 0) {
  snp_matrix[snp_matrix == 9] <- NA
  col_mode <- apply(snp_matrix, 2, function(x) {
    t <- table(x); if (length(t) == 0) 0 else as.numeric(names(t)[which.max(t)]) })
  na_idx <- which(is.na(snp_matrix), arr.ind = TRUE)
  snp_matrix[na_idx] <- col_mode[na_idx[, 2]]
  cat("Imputed", n_miss, "missing genotypes with the per-SNP most common genotype\n")
}

snp_positions <- read.table(pos_file, header = FALSE, stringsAsFactors = FALSE)
if (nrow(snp_positions) != n_snps)                     # pos file must match geno
  stop("Pos file rows (", nrow(snp_positions), ") != SNPs in geno (", n_snps, ")")
snp_positions <- snp_positions[, 1:3]
colnames(snp_positions) <- c("V1", "V2", "V3")         # scaffold, bp, chromosome

#### Sample info and climate ####
sample_info <- read.csv(csv_file, header = TRUE, sep = ";", dec = ".")
sample_info$parasit <- as.numeric(gsub(",", ".", sample_info$parasit))  # decimal comma -> point
if (nrow(sample_info) != n_inds)
  stop("Sample info rows (", nrow(sample_info), ") != individuals in geno (", n_inds, ")")

#### REORDER sample_info TO MATCH VCF (= .geno) SAMPLE ORDER ####
# Read sample IDs from the VCF #CHROM header (columns 10 onward)
vcf_ids <- system(paste0("grep -m1 '^#CHROM' ", shQuote(vcf_file), " | cut -f10- | tr '\\t' '\\n'"),
                  intern = TRUE)
vcf_ids <- trimws(vcf_ids); vcf_ids <- vcf_ids[nzchar(vcf_ids)]
csv_ids <- trimws(as.character(sample_info[[1]]))
# Safety checks: counts, duplicates, identical ID sets
if (length(vcf_ids) != n_inds) stop("VCF samples (", length(vcf_ids), ") != individuals in geno (", n_inds, ")")
if (anyDuplicated(vcf_ids) || anyDuplicated(csv_ids)) stop("Duplicate sample IDs")
if (!setequal(vcf_ids, csv_ids)) stop("Sample IDs differ between VCF and CSV: ",
    paste(head(setdiff(vcf_ids, csv_ids), 5), collapse = ","), " | ",
    paste(head(setdiff(csv_ids, vcf_ids), 5), collapse = ","))
n_moved <- sum(vcf_ids != csv_ids)                     # how many rows change position
sample_info <- sample_info[match(vcf_ids, csv_ids), ]  # actual reordering
rownames(sample_info) <- NULL
stopifnot(identical(trimws(as.character(sample_info[[1]])), vcf_ids))
cat("Reordered sample_info to VCF order; rows that moved:", n_moved, "\n")
write.csv(data.frame(vcf_order = seq_along(vcf_ids), id = vcf_ids, moved = vcf_ids != csv_ids),
          file.path(out_dir, "sample_order_check.csv"), row.names = FALSE)  # audit trail

# Extract the bioclimatic variables at the sample coordinates
predNames <- paste0("bio", 1:19)
presClim <- stack(paste0(bio_dir, "wc2.1_30s_bio_", 1:19, ".tif"))
names(presClim) <- predNames
climate_data <- raster::extract(presClim, y = sample_info[, c("lon", "lat")])
sample_info2 <- cbind(sample_info, climate_data) %>% mutate(bio4 = bio4 / 100)  # WorldClim bio4 is stored x100

missing_preds <- setdiff(predNames_analysis, colnames(sample_info2))
if (length(missing_preds) > 0) stop("Missing predictors: ", paste(missing_preds, collapse = ", "))
if (anyNA(sample_info2[, predNames_analysis])) stop("NAs in predictors")

# A "site" = unique lon_lat pair
pop <- paste(sample_info2$lon, sample_info2$lat, sep = "_")
n_pop <- length(unique(pop))

#### Scale predictors (mean 0, sd 1) ####
sample_info2[, predNames_analysis] <- scale(sample_info2[, predNames_analysis])
X <- as.matrix(sample_info2[, predNames_analysis])

#### LFMM function ####
# Ridge LFMM with K = 1; p-values calibrated with GIF
run_lfmm <- function(Xmat) {
  m  <- lfmm_ridge(Y = snp_matrix, X = Xmat, K = 1)
  mt <- lfmm_test(Y = snp_matrix, X = Xmat, lfmm = m)
  p <- mt$calibrated.pvalue
  colnames(p) <- colnames(Xmat)
  list(p = p, gif = mt$gif)
}

#### Observed LFMM ####
cat("Running observed LFMM...\n")
obs_lfmm <- run_lfmm(X)
adjusted_p <- obs_lfmm$p                               # SNPs x predictors
candidates <- lapply(predNames_analysis, function(p) get_candidates(adjusted_p[, p], fdr))
names(candidates) <- predNames_analysis
obs_n    <- sapply(candidates, length)                 # observed number of FDR candidates
obs_minp <- apply(adjusted_p, 2, min, na.rm = TRUE)    # observed minimum p per predictor
for (pred in predNames_analysis) cat(sprintf("  LFMM %-8s: %4d candidates\n", pred, obs_n[[pred]]))

# Sanity check: the first candidate must be the SNP with the minimum p
for (pred in predNames_analysis) {
  if (length(candidates[[pred]]) > 0)
    stopifnot(adjusted_p[candidates[[pred]][1], pred] == obs_minp[[pred]])
}

#### Site-level JOINT permutations (generated here, from the REORDERED data) ####
# Site-level predictor table (predictors are constant within site; site mean used)
pop_f    <- factor(pop, levels = unique(pop))
site_idx <- as.integer(pop_f)                          # site of each individual
site_X   <- apply(X, 2, function(v) tapply(v, pop_f, mean))   # sites x predictors
if (is.null(dim(site_X))) stop("Site table failed")
within_sd <- apply(X, 2, function(v) max(tapply(v, pop_f, sd), na.rm = TRUE))
if (any(within_sd > 1e-8)) cat("NOTE: predictors vary within site (site means used):",
                               paste(names(within_sd)[within_sd > 1e-8], collapse = ", "), "\n")

# Joint: shuffle whole site-level predictor sets together (preserves correlation
# among predictors), then expand back to individuals
perm_X_joint <- function() site_X[sample(nrow(site_X)), , drop = FALSE][site_idx, , drop = FALSE]

# Cache is reused only if signature + CSV checksum + structure all validate
run_perms <- function(n_perm, perm_fun, seed0, file, label, signature) {
  manifest_file <- paste0(file, ".rds")                # manifest holding signature and CSV md5
  expected_cols <- c(paste0("n_", predNames_analysis), paste0("minp_", predNames_analysis))

  if (file.exists(file) && file.exists(manifest_file)) {
    manifest <- tryCatch(readRDS(manifest_file), error = function(e) {
      cat("Cannot read permutation cache manifest:", conditionMessage(e), "\n"); NULL })
    csv_md5 <- unname(tools::md5sum(file))
    cache_matches <- !is.null(manifest) &&
      identical(manifest$signature, signature) &&      # same inputs/parameters
      identical(manifest$csv_md5, csv_md5)             # CSV not modified

    if (cache_matches) {
      old <- tryCatch(read.csv(file, check.names = FALSE), error = function(e) {
        cat("Cannot read permutation cache:", conditionMessage(e), "\n"); NULL })
      # Structural validation: rows, columns, numeric type, finite values
      valid_data <- !is.null(old) &&
        nrow(old) == n_perm &&
        identical(names(old), expected_cols) &&
        all(vapply(old, is.numeric, logical(1))) &&
        all(is.finite(as.matrix(old)))
      if (valid_data) {
        cat("Reusing validated", label, "permutations from", file, "\n")
        return(old)
      }
    }
    cat("Existing", label, "permutation cache failed validation; regenerating.\n")
  } else if (file.exists(file) || file.exists(manifest_file)) {
    cat("Incomplete", label, "permutation cache; regenerating.\n")
  }

  cat("Running", n_perm, label, "permutations...\n")
  res <- matrix(NA_real_, nrow = n_perm, ncol = 2 * length(predNames_analysis))
  colnames(res) <- expected_cols
  t0 <- Sys.time()
  for (i in seq_len(n_perm)) {
    set.seed(seed0 + i)                                # per-permutation seed
    Xp <- perm_fun(); colnames(Xp) <- predNames_analysis
    pp <- run_lfmm(Xp)$p
    for (pred in predNames_analysis) {
      res[i, paste0("n_", pred)]    <- length(get_candidates(pp[, pred], fdr))  # null candidate count
      res[i, paste0("minp_", pred)] <- min(pp[, pred], na.rm = TRUE)            # null minimum p
    }
    if (i %% 100 == 0)
      cat("  ", label, i, "/", n_perm, "(", round(difftime(Sys.time(), t0, units = "mins"), 1), "min )\n")
  }
  df <- as.data.frame(res)
  write.csv(df, file, row.names = FALSE)
  saveRDS(list(signature = signature, csv_md5 = unname(tools::md5sum(file))), manifest_file)
  df
}

joint_file <- file.path(out_dir, "permutation_joint_site_lfmm.csv")
# cache_version: increment to force regeneration if the permutation/LFMM procedure changes
perm_signature <- object_md5(list(
  cache_version = 1L,
  snp_matrix = snp_matrix,
  site_X = site_X,
  site_idx = site_idx,
  predictors = predNames_analysis,
  fdr = fdr,
  K = 1L,
  n_perms = n_perms_joint,
  seed = joint_seed
))
null_df <- run_perms(n_perms_joint, perm_X_joint, joint_seed, joint_file, "joint",
                     signature = perm_signature)
cat("\nPermutations ready (saved in", out_dir, ")\n")

#### Per-SNP permutation support ####
snp_tab <- list()
for (pred in predNames_analysis) {
  p <- adjusted_p[, pred]
  snp_tab[[pred]] <- data.frame(
    snp_index = seq_len(n_snps), scaffold = snp_positions$V1, bp = snp_positions$V2,
    chr_id = snp_positions$V3, predictor = pred, p_calibrated = p,
    fdr_candidate = seq_len(n_snps) %in% candidates[[pred]],
    p_perm_joint  = perm_p(p, null_df[[paste0("minp_", pred)]]))  # compared with null min-p (accounts for multiple testing)
}
snp_all <- do.call(rbind, snp_tab)
snp_all$joint_support <- snp_all$fdr_candidate & snp_all$p_perm_joint < alpha
snp_all <- snp_all[!is.na(snp_all$p_calibrated), ]

# Outputs: all FDR candidates, those with joint support, and one file per predictor
cand_tab <- snp_all[snp_all$fdr_candidate, ]
cand_tab <- cand_tab[order(cand_tab$predictor, cand_tab$p_calibrated), ]
write.csv(cand_tab, file.path(out_dir, "candidates_with_permutation_support.csv"), row.names = FALSE)
write.csv(cand_tab[cand_tab$joint_support, ],
          file.path(out_dir, "candidates_JOINT_support.csv"), row.names = FALSE)
for (pred in predNames_analysis)
  write.csv(cand_tab[cand_tab$predictor == pred, ],
            file.path(out_dir, paste0("candidates_", pred, ".csv")), row.names = FALSE)

cat("\nFDR candidates and supported SNPs (joint permutation p <", alpha, "):\n")
print(data.frame(Predictor = predNames_analysis,
                 FDR_candidates = sapply(predNames_analysis, function(p) sum(cand_tab$predictor == p)),
                 Support_joint  = sapply(predNames_analysis, function(p) sum(cand_tab$predictor == p & cand_tab$joint_support))),
      row.names = FALSE)

#### Manhattan data (for the RDA script) + plots ####
manhattan_df <- data.frame(snp_positions, SNP = as.character(seq_len(n_snps)), as.data.frame(adjusted_p))
write.csv(manhattan_df, file.path(out_dir, "manhattan_data.csv"), row.names = FALSE)

for (model in predNames_analysis) {
  tryCatch({                                           # a plot error won't stop the script
    md <- manhattan_df[!is.na(manhattan_df[[model]]), ]
    png(file.path(out_dir, paste0("manhattan_", model, ".png")), width = 1400, height = 600)
    hl <- as.character(candidates[[model]])            # highlighted candidate SNPs
    manhattan(md, chr = "V3", bp = "V2", snp = "SNP", p = model,
              suggestiveline = TRUE,
              # red line = alpha quantile of the null min-p distribution
              genomewideline = -log10(quantile(null_df[[paste0("minp_", model)]], alpha)),
              highlight = if (length(hl) > 0) hl else NULL,
              main = paste("Manhattan -", model, "(het", het_filter, "; red line =", alpha * 100, "% permutation threshold)"))
    dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE)
    cat("Manhattan failed for", model, ":", conditionMessage(e), "\n") })
}

#### Count-based summary (joint design) ####
# For each predictor, compare observed vs null:
#  P_count = empirical p on the number of candidates; P_minp = empirical p on the minimum p
mk <- function(df, design, nperm) do.call(rbind, lapply(predNames_analysis, function(pred) {
  nd <- df[[paste0("n_", pred)]]; pd <- df[[paste0("minp_", pred)]]
  data.frame(Filter = het_filter, Design = design, N_perms = nperm, N_SNPs = n_snps, N_sites = n_pop,
             Predictor = pred, Observed_candidates = obs_n[[pred]],
             Null_mean = round(mean(nd), 2), Null_95pct = round(unname(quantile(nd, 0.95)), 2),
             Null_max = max(nd),
             P_count = round((sum(nd >= obs_n[[pred]]) + 1) / (nperm + 1), 4),
             Observed_min_p = signif(obs_minp[[pred]], 4),
             P_minp = round((sum(pd <= obs_minp[[pred]]) + 1) / (nperm + 1), 4))
}))
perm_summary <- mk(null_df, "joint", n_perms_joint)
write.csv(perm_summary, file.path(out_dir, "permutation_test_summary.csv"), row.names = FALSE)

cat("\n===== DONE (alpha =", alpha, "): het filter", het_filter, "=====\n")
