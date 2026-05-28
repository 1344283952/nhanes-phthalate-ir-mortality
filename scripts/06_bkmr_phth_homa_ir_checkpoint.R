# ============================================
# 009 / scripts/06_bkmr_phth_homa_ir_checkpoint.R
# Checkpointed BKMR for Phthalate 8 metabolites → log-HOMA-IR continuous
# 基于 007/06_bkmr_pfas_checkpoint.R 改 outcome 为 HOMA-IR
# 含 rolling-3 cleanup (feedback_bkmr_history_rolling_3 强制)
#
# 实验设计 (Stack 1 主):
#   Y = log(HOMA-IR)  (continuous, primary)
#   Z = 8 Phthalate metabolites z-scores
#       (MEP/MBP/MIB/MZP/MHP/MHH/MOH/ECP)
#   X = age, sex, race(4 dummies), edu(2), pir, smoke_ever, cotinine_log, bmi, waist
#   n ≈ 2,239 (Stack 1 主)
#
# Checkpoint: 100 iter/block × 100 blocks = 10,000 iter/chain × 2 chains
# Wall time estimate: n=2.2k 比 007 n=8.6k 快 ~3-4x → ~3-5 min/block → ~5-8 hours/chain
# Total: 10-16 hours (~ 1 day) for both chains complete
# ============================================

RNGkind("L'Ecuyer-CMRG")

suppressPackageStartupMessages({
  library(dplyr); library(bkmr); library(coda); library(rstan); library(ggplot2)
})

cat("========================================\n")
cat("009 BKMR (checkpoint) — 8 Phthalate -> log-HOMA-IR (Stack 1 主)\n")
cat("========================================\n\n")

# ---- Constants ----
SENSITIVITY_MODE <- Sys.getenv("BKMR_SENSITIVITY_MODE", "primary")
if (!SENSITIVITY_MODE %in% c("primary","knots_K50")) {
  stop("Unknown SENSITIVITY_MODE: ", SENSITIVITY_MODE)
}
CHECKPOINT_FILE <- switch(SENSITIVITY_MODE,
  primary    = "data/processed/bkmr_phth_homa_009_checkpoint.rds",
  knots_K50  = "data/processed/bkmr_phth_homa_009_sens_K50.rds"
)
BLOCK_ITER  <- as.integer(Sys.getenv("BKMR_BLOCK_ITER", "100"))
TARGET_ITER <- as.integer(Sys.getenv("BKMR_TARGET_ITER", "10000"))
N_CHAINS    <- 2L
BURN_IN     <- as.integer(Sys.getenv("BKMR_BURN_IN", as.character(min(5000L, TARGET_ITER %/% 2L))))
SEED_BASE   <- 20260523L
KNOTS_K     <- as.integer(Sys.getenv("BKMR_KNOTS_K",
                                     if (SENSITIVITY_MODE == "primary") "100" else "50"))

cat(sprintf("Config: MODE=%s  BLOCK_ITER=%d  TARGET_ITER=%d  N_CHAINS=%d  BURN_IN=%d  KNOTS=%d\n",
            SENSITIVITY_MODE, BLOCK_ITER, TARGET_ITER, N_CHAINS, BURN_IN, KNOTS_K))
cat(sprintf("        CHECKPOINT=%s\n\n", CHECKPOINT_FILE))

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Load + prep ----
load("data/processed/nhanes_final.RData")

# Phthalate z-scores (03_clean 已生成 URX*_z)
phth_z_cols <- c("URXMEP_z","URXMBP_z","URXMIB_z","URXMZP_z",
                 "URXMHP_z","URXMHH_z","URXMOH_z","URXECP_z")
phth_z_cols <- intersect(phth_z_cols, names(nhanes_final))
phth_labels <- c("MEP","MnBP","MiBP","MBzP","MEHP","MEHHP","MEOHP","MECPP")[seq_along(phth_z_cols)]

if (length(phth_z_cols) < 4) {
  stop("Phthalate z-score 列不足 (n=", length(phth_z_cols), "), 检查 03_clean_data.R")
}

cat(sprintf("Phthalate 化合物: %s\n", paste(phth_z_cols, collapse = ", ")))

# Y = log(HOMA-IR)
nhanes_final$y_outcome <- nhanes_final$homa_ir_log

# X covariates
nhanes_final$sex_male_i   <- as.integer(nhanes_final$RIAGENDR == 1)
nhanes_final$race_nhw     <- as.integer(as.character(nhanes_final$race) == "Non-Hispanic White")
nhanes_final$race_nhb     <- as.integer(as.character(nhanes_final$race) == "Non-Hispanic Black")
nhanes_final$race_mex     <- as.integer(as.character(nhanes_final$race) == "Mexican American")
nhanes_final$race_othhisp <- as.integer(as.character(nhanes_final$race) == "Other Hispanic")
nhanes_final$edu_lths     <- as.integer(as.character(nhanes_final$education) == "Less than HS")
nhanes_final$edu_hs       <- as.integer(as.character(nhanes_final$education) == "High school")
nhanes_final$smoke_ever_i <- as.integer(as.character(nhanes_final$smoke) == "Ever")

X_covars <- c("age","sex_male_i",
              "race_nhw","race_nhb","race_mex","race_othhisp",
              "edu_lths","edu_hs","pir","bmi","waist","smoke_ever_i")

bkmr_vars_keep <- c(phth_z_cols, "y_outcome", X_covars, "wt_pooled")
df_bkmr <- nhanes_final %>%
  select(any_of(bkmr_vars_keep)) %>%
  filter(if_all(everything(), ~ !is.na(.)))

cat(sprintf("\n009 Phthalate BKMR analytic n = %d (mode=%s)\n", nrow(df_bkmr), SENSITIVITY_MODE))

Z_mat <- as.matrix(df_bkmr[, phth_z_cols])
colnames(Z_mat) <- phth_labels
y_vec <- df_bkmr$y_outcome
X_mat <- as.matrix(df_bkmr[, X_covars])

# ---- Load or init checkpoint ----
if (file.exists(CHECKPOINT_FILE)) {
  cp <- readRDS(CHECKPOINT_FILE)
  cat(sprintf("[RESUME] done iter/chain: %s / %d\n",
              paste(cp$done_iter_per_chain, collapse = ", "), TARGET_ITER))
  cat(sprintf("[RESUME] block count: %d  started_at: %s\n\n", cp$block_count, cp$started_at))
} else {
  cp <- list(
    done_iter_per_chain = rep(0L, N_CHAINS),
    last_states         = vector("list", N_CHAINS),
    accumulated_fits    = vector("list", N_CHAINS),
    knot_grids          = vector("list", N_CHAINS),
    block_count         = 0L,
    started_at          = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    target_iter         = TARGET_ITER,
    block_iter          = BLOCK_ITER,
    seed_base           = SEED_BASE,
    knots_k             = KNOTS_K
  )
  cat("[FRESH] No checkpoint found. Initializing new run.\n\n")
}

all_done <- all(cp$done_iter_per_chain >= TARGET_ITER)

extract_last_state <- function(fit) {
  n_keep <- if (!is.null(fit$beta))
              { if (is.matrix(fit$beta)) nrow(fit$beta) else length(fit$beta) }
            else length(fit$sigsq.eps)
  if (n_keep == 0) return(NULL)
  list(
    beta      = if (!is.null(fit$beta))
                  { if (is.matrix(fit$beta)) fit$beta[n_keep, ] else fit$beta[n_keep] } else NULL,
    lambda    = if (!is.null(fit$lambda))
                  { if (is.matrix(fit$lambda)) fit$lambda[n_keep, ] else fit$lambda[n_keep] } else NULL,
    r         = if (!is.null(fit$r))
                  { if (is.matrix(fit$r)) fit$r[n_keep, ] else fit$r[n_keep] } else NULL,
    sigsq.eps = fit$sigsq.eps[n_keep],
    delta     = if (!is.null(fit$delta))
                  { if (is.matrix(fit$delta)) fit$delta[n_keep, ] else fit$delta[n_keep] } else NULL
  )
}

concat_fits <- function(prev, new) {
  if (is.null(prev)) return(new)
  out <- new
  rbind_safe <- function(a, b) {
    if (is.null(a)) return(b); if (is.null(b)) return(a)
    if (is.matrix(a) && is.matrix(b)) rbind(a, b)
    else if (is.matrix(a)) rbind(a, matrix(b, nrow = 1))
    else if (is.matrix(b)) rbind(matrix(a, nrow = 1), b)
    else rbind(matrix(a, nrow = 1), matrix(b, nrow = 1))
  }
  for (p in c("beta","lambda","r","delta","h.hat","ystar")) {
    if (!is.null(prev[[p]]) || !is.null(new[[p]]))
      out[[p]] <- rbind_safe(prev[[p]], new[[p]])
  }
  out$sigsq.eps <- c(prev$sigsq.eps %||% numeric(0), new$sigsq.eps %||% numeric(0))
  out$iter <- (prev$iter %||% nrow(prev$r %||% prev$beta) %||% length(prev$sigsq.eps)) +
              (new$iter  %||% nrow(new$r  %||% new$beta)  %||% length(new$sigsq.eps))
  out
}

# ---- Run one block per chain ----
if (!all_done) {
  cp$block_count <- cp$block_count + 1L
  block_t0 <- Sys.time()
  cat(sprintf("[BLOCK %d] starting %s\n",
              cp$block_count, format(block_t0, "%Y-%m-%d %H:%M:%S")))

  for (c_i in seq_len(N_CHAINS)) {
    if (cp$done_iter_per_chain[c_i] >= TARGET_ITER) {
      cat(sprintf("  [chain %d] already at target, skip\n", c_i)); next
    }
    remaining <- TARGET_ITER - cp$done_iter_per_chain[c_i]
    this_iter <- min(BLOCK_ITER, remaining)
    t0 <- Sys.time()
    cat(sprintf("\n  [chain %d] running %d iter (cumulative %d -> %d / %d) at %s\n",
                c_i, this_iter, cp$done_iter_per_chain[c_i],
                cp$done_iter_per_chain[c_i] + this_iter, TARGET_ITER,
                format(t0, "%H:%M:%S")))
    flush(stdout())

    if (is.null(cp$knot_grids[[c_i]])) {
      set.seed(SEED_BASE + c_i * 10L)
      knot_idx <- sample(seq_len(nrow(Z_mat)), min(KNOTS_K, nrow(Z_mat)))
      cp$knot_grids[[c_i]] <- Z_mat[knot_idx, , drop = FALSE]
      cat(sprintf("  [chain %d] knot grid init (K=%d)\n", c_i, KNOTS_K))
    }
    knot_grid <- cp$knot_grids[[c_i]]
    set.seed(SEED_BASE + c_i * 1000000L + cp$done_iter_per_chain[c_i])

    fit_block <- bkmr::kmbayes(
      y = y_vec, Z = Z_mat, X = X_mat,
      iter = this_iter, family = "gaussian",
      verbose = TRUE, varsel = TRUE,
      knots = knot_grid,
      starting.values = cp$last_states[[c_i]],
      control.params = list(lambda.jump = 10, mu.r = 5, sigma.r = 25,
                            a.p0 = 1, b.p0 = 1, a.sigsq = 1e-3, b.sigsq = 1e-3),
      est.h = TRUE
    )

    t1 <- Sys.time()
    elapsed_min <- as.numeric(difftime(t1, t0, units = "mins"))
    cat(sprintf("  [chain %d] block done in %.1f min (%.1f iter/min)\n",
                c_i, elapsed_min, this_iter / elapsed_min))

    cp$last_states[[c_i]]       <- extract_last_state(fit_block)
    cp$accumulated_fits[[c_i]]  <- concat_fits(cp$accumulated_fits[[c_i]], fit_block)
    cp$done_iter_per_chain[c_i] <- cp$done_iter_per_chain[c_i] + this_iter

    tmp <- paste0(CHECKPOINT_FILE, ".tmp")
    saveRDS(cp, tmp); file.rename(tmp, CHECKPOINT_FILE)
    cat(sprintf("  [chain %d] checkpoint saved (%d / %d, %.1f%%)\n",
                c_i, cp$done_iter_per_chain[c_i], TARGET_ITER,
                100 * cp$done_iter_per_chain[c_i] / TARGET_ITER))

    # ---- History backup + ROLLING-3 cleanup (feedback_bkmr_history_rolling_3) ----
    tryCatch({
      hist_dir <- "_bkmr_checkpoint_history_phth"
      if (!dir.exists(hist_dir)) dir.create(hist_dir, recursive = TRUE)
      hist_path <- file.path(hist_dir,
        sprintf("block_%03d_chain%d_%s.rds", cp$block_count, c_i,
                format(Sys.time(), "%Y%m%d_%H%M%S")))
      file.copy(CHECKPOINT_FILE, hist_path, overwrite = FALSE)

      # Rolling-3: keep newest 3 history files, delete older (防 269 GB 累积)
      all_snaps <- list.files(hist_dir, pattern = "^block_.*\\.rds$", full.names = TRUE)
      if (length(all_snaps) > 3) {
        snap_info <- file.info(all_snaps)
        snap_info <- snap_info[order(snap_info$mtime, decreasing = TRUE), ]
        to_delete <- rownames(snap_info)[-(1:3)]
        n_deleted <- length(to_delete)
        file.remove(to_delete)
        cat(sprintf("  [chain %d] history rolling-3: kept newest 3, deleted %d old\n",
                    c_i, n_deleted))
      }
    }, error = function(e) {})

    # ---- Block history CSV ----
    tryCatch({
      blk_csv <- "_bkmr_block_history_009.csv"
      blk_row <- data.frame(
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        block_id = cp$block_count, chain_id = c_i,
        iter_in_block = this_iter,
        cumulative_iter = cp$done_iter_per_chain[c_i],
        target_iter = TARGET_ITER,
        block_wall_min = round(elapsed_min, 3),
        knots_K = KNOTS_K, burn_in = BURN_IN,
        stringsAsFactors = FALSE)
      write.table(blk_row, file = blk_csv,
                  append = file.exists(blk_csv),
                  col.names = !file.exists(blk_csv),
                  row.names = FALSE, sep = ",", quote = FALSE)
    }, error = function(e) {})

    rm(fit_block); invisible(gc(verbose = FALSE))
  }

  block_t1 <- Sys.time()
  cat(sprintf("\n[BLOCK %d DONE] total %.1f min\n",
              cp$block_count, as.numeric(difftime(block_t1, block_t0, units = "mins"))))
  cat("\nPer-chain progress:\n")
  for (c_i in seq_len(N_CHAINS)) {
    pct <- 100 * cp$done_iter_per_chain[c_i] / TARGET_ITER
    cat(sprintf("  Chain %d: %5d / %5d  (%.1f%%)\n",
                c_i, cp$done_iter_per_chain[c_i], TARGET_ITER, pct))
  }

  all_done <- all(cp$done_iter_per_chain >= TARGET_ITER)
}

# ---- Finalize ----
if (all_done && !file.exists("output/tables/bkmr_phth_results.RData")) {
  cat("\n========================================\n")
  cat("FINALIZE: convergence + PIP + figures\n")
  cat("========================================\n\n")

  bkmr_fit <- list(fits = cp$accumulated_fits)
  n_iter <- TARGET_ITER; n_burn <- BURN_IN; n_chain <- N_CHAINS

  if (!dir.exists("output/tables"))  dir.create("output/tables",  recursive = TRUE)
  if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

  # Convergence diag
  cat("[1/3] Convergence...\n")
  extract_chain_param <- function(fit_obj, pname, ci, after_burn = n_burn) {
    fit <- fit_obj$fits[[ci]]; vals <- fit[[pname]]
    if (is.null(vals)) return(NULL)
    if (is.matrix(vals)) vals[(after_burn + 1):n_iter, , drop = FALSE]
    else vals[(after_burn + 1):n_iter]
  }
  check_params <- list(sigsq.eps = "sigsq.eps", lambda = "lambda")
  for (j in seq_along(phth_labels))
    check_params[[paste0("r_", phth_labels[j])]] <- list(name = "r", col = j)

  diag_table <- data.frame()
  for (pname in names(check_params)) {
    p <- check_params[[pname]]
    chains_data <- if (is.list(p) && !is.null(p$col)) {
      lapply(seq_len(n_chain), function(ci) {
        m <- extract_chain_param(bkmr_fit, p$name, ci)
        if (is.null(m)) return(rep(NA, n_iter - n_burn)); m[, p$col]
      })
    } else {
      lapply(seq_len(n_chain), function(ci) extract_chain_param(bkmr_fit, p, ci))
    }
    chains_data <- chains_data[!sapply(chains_data, function(x) is.null(x) || all(is.na(x)))]
    if (length(chains_data) < 2) next
    M <- do.call(cbind, chains_data)
    rhat <- tryCatch(rstan::Rhat(M), error = function(e) NA)
    ess  <- tryCatch(coda::effectiveSize(coda::as.mcmc.list(
              lapply(chains_data, coda::as.mcmc)))[1], error = function(e) NA)
    pass <- !is.na(rhat) & !is.na(ess) & rhat < 1.1 & ess > 400
    diag_table <- rbind(diag_table,
      data.frame(param = pname, rhat = rhat, ess = as.numeric(ess), pass = pass))
    cat(sprintf("  %-20s rhat=%.3f  ESS=%.0f  %s\n",
                pname, rhat, ess, ifelse(isTRUE(pass), "OK", "FAIL")))
  }
  write.csv(diag_table, "output/tables/bkmr_phth_convergence.csv", row.names = FALSE)

  # PIP + summaries
  cat("\n[2/3] PIP + summaries...\n")
  fit_combined <- bkmr_fit$fits[[1]]
  pips <- tryCatch(bkmr::ExtractPIPs(fit_combined),
                   error = function(e) data.frame(variable = phth_labels, PIP = NA))
  write.csv(pips, "output/tables/bkmr_phth_pip.csv", row.names = FALSE)
  print(pips)

  univar <- tryCatch(
    bkmr::PredictorResponseUnivar(fit = fit_combined, q.fixed = 0.5,
      sel = seq(n_burn + 1, n_iter, by = 25), method = "approx"),
    error = function(e) NULL)
  bivar <- tryCatch(
    bkmr::PredictorResponseBivar(fit = fit_combined,
      z.pairs = data.frame(z1 = rep(1:7, each = 7), z2 = rep(2:8, times = 7)) %>%
                dplyr::filter(z1 < z2),
      q.fixed = 0.5, sel = seq(n_burn + 1, n_iter, by = 50), method = "approx"),
    error = function(e) NULL)
  overall <- tryCatch(
    bkmr::OverallRiskSummaries(fit = fit_combined, y = y_vec, Z = Z_mat, X = X_mat,
      qs = seq(0.1, 0.9, by = 0.1), q.fixed = 0.5, method = "approx",
      sel = seq(n_burn + 1, n_iter, by = 25)),
    error = function(e) NULL)

  cat("\n[3/3] Save bkmr_phth_results.RData\n")
  save(bkmr_fit, pips, univar, bivar, overall, diag_table, phth_labels,
       file = "output/tables/bkmr_phth_results.RData")
  cat("DONE.\n")
} else if (all_done) {
  cat("\n[ALREADY FINALIZED]\n")
}
