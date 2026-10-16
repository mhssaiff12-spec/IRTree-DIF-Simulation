# =====================================================================
# Stage 2 (fast, equivalent form): latent multiple-group LR DIF tests per IRTree node
# Because the three node traits are uncorrelated in the model and every pseudo-item
# loads on exactly one trait, the 3-D likelihood factorises into three 1-D node models.
# The LR test for a candidate pseudo-item therefore depends only on its own node model.
# Verified: 1-D LR = 8.9767 vs 3-D LR = 8.9767 (same data, same 15 quadrature points).
# NOTE: equivalence holds ONLY for uncorrelated node traits. If you add COV between
# traits (realistic for real data), go back to the full 3-D model in 02_latent_mirt_validation.R.
# =====================================================================
suppressMessages(library(mirt))
source("irtree_dif_sim.R")

fit_node <- function(Wn, groups, inv, quadpts) {
  multipleGroup(Wn, 1, group = groups, itemtype = "Rasch",
                invariance = c(inv, "free_means", "free_var"), method = "EM",
                quadpts = quadpts, verbose = FALSE, technical = list(NCYCLES = 1000))
}

# One case: all node tests for items in test_items, anchors = pure_anchors minus tested item
run_case_fast <- function(case_name, d, test_items = 1:6, pure_anchors = 7:20,
                          dif_node = "extreme", shift = 0.6, tree = "MDE", quadpts = 61) {
  P <- pseudo_items(d$Y, tree); J <- ncol(d$Y)
  groups <- factor(ifelse(d$group == 1, "focal", "reference"), levels = c("reference", "focal"))
  key <- c(middle = "M", direction = "D", extreme = "E")
  rows <- list()
  for (j in test_items) for (nd in names(key)) {
    Wn <- P[[nd]]; storage.mode(Wn) <- "numeric"; colnames(Wn) <- paste0(key[nd], seq_len(J))
    anc <- paste0(key[nd], setdiff(pure_anchors, j)); cand <- paste0(key[nd], j)
    z <- tryCatch({
      alt <- fit_node(Wn, groups, anc, quadpts); nul <- fit_node(Wn, groups, c(anc, cand), quadpts)
      conv <- isTRUE(extract.mirt(alt, "converged")) && isTRUE(extract.mirt(nul, "converged"))
      LR <- max(0, 2 * (as.numeric(logLik(alt)) - as.numeric(logLik(nul))))
      ddf <- extract.mirt(alt, "nest") - extract.mirt(nul, "nest")
      data.frame(item = j, node = key[nd], LR = LR, df = ddf, p = pchisq(LR, ddf, lower.tail = FALSE),
                 converged = conv, error = NA_character_)
    }, error = function(e) data.frame(item = j, node = key[nd], LR = NA, df = NA, p = NA,
                                       converged = FALSE, error = conditionMessage(e)))
    rows[[length(rows) + 1]] <- z
  }
  ans <- do.call(rbind, rows); rownames(ans) <- NULL
  ans$case <- case_name; ans$tree <- tree
  ans$true_dif_item <- shift != 0 & ans$item %in% 1:4
  ans$true_dif_node <- ans$true_dif_item & ans$node == key[dif_node]
  ans$p_BH_all <- p.adjust(ans$p, "BH"); ans$flag <- !is.na(ans$p_BH_all) & ans$p_BH_all < .05
  ans
}

CASES <- list(
  list(name = "A_noDIF_noImpact", seed = 20261014, dif_node = "direction", shift = 0,  it = 0,  ie = 0),
  list(name = "B_DIF_noImpact",   seed = 20261015, dif_node = "extreme",   shift = .6, it = 0,  ie = 0),
  list(name = "C_noDIF_impact",   seed = 20261016, dif_node = "direction", shift = 0,  it = .5, ie = .5),
  list(name = "D_DIF_impact",     seed = 20261017, dif_node = "extreme",   shift = .6, it = .5, ie = .5))

# legacy = TRUE reproduces the data actually used in the first stage-2 run (all cases on the
# seed-777 stream because the old make_items() reset the global seed). legacy = FALSE uses the
# declared case seeds (correct behaviour with the fixed make_items()).
run_all <- function(legacy = FALSE, quadpts = 61) {
  params <- make_items(20, seed = 777)
  do.call(rbind, lapply(CASES, function(q) {
    if (legacy) { set.seed(777); invisible(rnorm(60)) } else set.seed(q$seed)
    d <- sim_irtree(params, n = 1000, dif_items = 1:4, dif_node = q$dif_node, shift = q$shift,
                    impact_theta = q$it, impact_ers = q$ie)
    run_case_fast(q$name, d, dif_node = q$dif_node, shift = q$shift, quadpts = quadpts)
  }))
}

if (sys.nframe() == 0) {
  dir.create("stage2_output", showWarnings = FALSE)
  t0 <- Sys.time(); leg <- run_all(legacy = TRUE, quadpts = 15)
  write.csv(leg, "stage2_output/fast_legacy_q15.csv", row.names = FALSE)
  cat("legacy run time:", format(Sys.time() - t0), "\n")
  t0 <- Sys.time(); new <- run_all(legacy = FALSE, quadpts = 61)
  write.csv(new, "stage2_output/fast_newseeds_q61.csv", row.names = FALSE)
  cat("new-seed run time:", format(Sys.time() - t0), "\n")
}
