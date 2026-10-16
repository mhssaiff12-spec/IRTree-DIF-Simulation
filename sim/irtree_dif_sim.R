# =====================================================================
# IRTree DIF robustness study — Stages 0–4 pilot (Monte Carlo)
# Base R + MASS + parallel only. Seed-controlled and reproducible.
# Author of protocol: Maher; code revised/extended 2026-10-14.
# =====================================================================
suppressPackageStartupMessages({ library(MASS); library(parallel) })

## ---------- Stage 0: corrected data generator ------------------------
# True tree (MDE): Middle (3 vs rest) -> Direction (4/5 vs 1/2) -> Extremity (1/5 vs 2/4)
# theta drives direction; MRS drives middle; ERS drives extremity.
# Fixes relative to IRTree_pilot.R:
#  (1) item parameters fixed across replications (fixed truth);
#  (2) impact: focal group may differ in mean theta and mean ERS;
#  (3) anchor set built correctly: k DIF items deliberately included
#      in a nominated anchor set of fixed size (contamination);
#  (4) Monte Carlo loop with MCSE.
make_items <- function(J = 20L, seed = 777L) {
  # Draw fixed item parameters WITHOUT disturbing the caller's RNG stream
  # (earlier version called set.seed() globally, so every later simulation
  #  reused the seed-777 stream regardless of the caller's own set.seed()).
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  list(b_m = rnorm(J, -0.9, 0.25), b_d = rnorm(J, 0, 0.45), b_e = rnorm(J, -0.2, 0.25))
}

sim_irtree <- function(items, n, dif_items = 1:4, dif_node = "direction", shift = 0.6,
                       impact_theta = 0, impact_ers = 0, group_ratio = 0.5) {
  J <- length(items$b_d)
  g <- as.integer(seq_len(n) > round(n * group_ratio))          # 1 = focal
  theta <- rnorm(n, impact_theta * g); mrs <- rnorm(n); ers <- rnorm(n, impact_ers * g)
  Y <- matrix(NA_integer_, n, J)
  for (j in seq_len(J)) {
    s  <- if (j %in% dif_items) shift * g else 0
    M  <- rbinom(n, 1, plogis(items$b_m[j] + mrs   + if (dif_node == "middle")    s else 0))
    D  <- rbinom(n, 1, plogis(items$b_d[j] + theta + if (dif_node == "direction") s else 0))
    E  <- rbinom(n, 1, plogis(items$b_e[j] + ers   + if (dif_node == "extreme")   s else 0))
    Y[, j] <- ifelse(M == 1, 3L, ifelse(D == 1, ifelse(E == 1, 5L, 4L), ifelse(E == 1, 1L, 2L)))
  }
  list(Y = Y, group = g)
}

nominate_anchors <- function(J, dif_items, size = 10L, k_contam = 0L) {
  clean <- setdiff(seq_len(J), dif_items)
  sort(c(head(dif_items, k_contam), head(clean, size - k_contam)))
}

## ---------- Pseudo-item coding for each tree --------------------------
# Returns list of n x J matrices (NA = node not reached)
pseudo_items <- function(Y, tree = c("MDE", "EMD")) {
  tree <- match.arg(tree)
  if (tree == "MDE") {                    # correct tree
    mid <- (Y == 3) * 1
    dir <- ifelse(Y == 3, NA, (Y >= 4) * 1)
    ext <- ifelse(Y == 3, NA, (Y %in% c(1, 5)) * 1)
  } else {                                # misspecified tree: extremity first
    ext <- (Y == 1 | Y == 5) * 1
    mid <- ifelse(Y %in% c(1, 5), NA, (Y == 3) * 1)
    dir <- ifelse(Y == 3, NA, (Y >= 4) * 1)
  }
  dim(ext) <- dim(mid) <- dim(dir) <- dim(Y)
  list(middle = mid, direction = dir, extreme = ext)
}

## ---------- Fast logistic Wald test -------------------------------------
wald_glm <- function(z, x, g) {
  ok <- !is.na(z); z <- z[ok]; X <- cbind(1, x[ok], g[ok])
  if (length(unique(z)) < 2 || length(unique(X[, 3])) < 2) return(c(NA, NA, 0))
  f <- suppressWarnings(glm.fit(X, z, family = binomial()))
  bad <- as.numeric(!f$converged || f$boundary)
  V <- tryCatch(chol2inv(f$qr$qr[1:3, 1:3]), error = function(e) NULL)
  if (is.null(V)) return(c(NA, NA, 1))
  est <- f$coefficients[3]; se <- sqrt(V[3, 3])
  c(est, 2 * pnorm(-abs(est / se)), bad)
}

# Node-wise DIF for one anchor set: matching = logit of anchor proportion at same node
tree_dif_once <- function(P, g, anchors) {
  J <- ncol(P[[1]]); out <- list()
  for (nd in names(P)) {
    M <- P[[nd]]; A <- M[, anchors, drop = FALSE]
    sA <- rowSums(A, na.rm = TRUE); cA <- rowSums(!is.na(A))
    res <- t(vapply(seq_len(J), function(j) {
      s <- sA; c <- cA
      if (j %in% anchors) { v <- M[, j]; s <- s - ifelse(is.na(v), 0, v); c <- c - !is.na(v) }
      wald_glm(M[, j], qlogis((s + .5) / (c + 1)), g)
    }, numeric(3)))
    out[[nd]] <- res
  }
  est <- sapply(out, `[`, , 1); p <- sapply(out, `[`, , 2); bad <- sapply(out, `[`, , 3)
  padj <- matrix(p.adjust(p, "BH"), nrow(p), dimnames = dimnames(p))     # BH over all J x 3 tests
  list(est = est, p = padj, flag_node = !is.na(padj) & padj < .05,
       flag = rowSums(!is.na(padj) & padj < .05) > 0, nbad = sum(bad), nfit = length(bad))
}

ord_dif_once <- function(Y, g, anchors) {
  J <- ncol(Y); TA <- rowSums(Y[, anchors, drop = FALSE])
  res <- t(vapply(seq_len(J), function(j) {
    Tj <- if (j %in% anchors) TA - Y[, j] else TA
    f <- tryCatch(suppressWarnings(polr(factor(Y[, j]) ~ Tj + g, Hess = TRUE)), error = function(e) NULL)
    if (is.null(f)) return(c(NA, NA, 1))
    cf <- coef(summary(f))["g", ]
    c(cf[1], 2 * pnorm(-abs(cf[3])), as.numeric(f$convergence != 0))
  }, numeric(3)))
  padj <- p.adjust(res[, 2], "BH")
  list(est = res[, 1], p = padj, flag = !is.na(padj) & padj < .05, nbad = sum(res[, 3]), nfit = J)
}

# Purification (pre-specified): anchor_{t+1} = nominated anchors NOT flagged at round t.
# Purification may only REMOVE suspected items from the nominated set; it never adds
# untested items (an "all non-flagged items" rule let DIF items flow back into the
# anchor whenever a round flagged nothing, causing runaway contamination).
run_with_strategy <- function(fun, data_arg, g, nominated, strategy, max_iter = 5L) {
  r <- fun(data_arg, g, nominated)
  if (strategy == "nominated") return(c(r, iters = 1L))
  anc <- nominated; it <- 1L; nbad <- r$nbad; nfit <- r$nfit
  repeat {
    new_anc <- setdiff(nominated, which(r$flag))
    if (length(new_anc) < 3 || identical(new_anc, anc) || it >= max_iter) break
    anc <- new_anc; r <- fun(data_arg, g, anc); it <- it + 1L
    nbad <- nbad + r$nbad; nfit <- nfit + r$nfit
  }
  r$nbad <- nbad; r$nfit <- nfit; c(r, iters = it)
}

## ---------- One replication -> summary rows ----------------------------
one_rep <- function(cond, items, rep_id, dif_items = 1:4) {
  J <- length(items$b_d)
  shift <- if (cond$node == "none") 0 else cond$shift
  ratio <- if (!is.null(cond$ratio)) cond$ratio else 0.5
  dat <- sim_irtree(items, cond$N, dif_items, if (cond$node == "none") "direction" else cond$node,
                    shift, cond$imp_theta, cond$imp_ers, group_ratio = ratio)
  truth <- if (shift == 0) rep(FALSE, J) else seq_len(J) %in% dif_items
  anchors <- nominate_anchors(J, dif_items, 10L, cond$k)
  rows <- list()
  for (meth in c("tree_correct", "tree_wrong", "ordinal")) {
    for (strat in c("nominated", "purified")) {
      t0 <- proc.time()[["elapsed"]]
      r <- if (meth == "ordinal") run_with_strategy(ord_dif_once, dat$Y, dat$group, anchors, strat)
           else run_with_strategy(tree_dif_once, pseudo_items(dat$Y, if (meth == "tree_correct") "MDE" else "EMD"),
                                  dat$group, anchors, strat)
      secs <- proc.time()[["elapsed"]] - t0
      V <- sum(r$flag & !truth); S <- sum(r$flag & truth); R <- V + S
      loc <- est_true <- NA_real_
      if (meth != "ordinal" && any(truth)) {
        loc <- mean(r$flag_node[truth, cond$node])                       # flagged at node with the true label
        est_true <- mean(r$est[truth, cond$node], na.rm = TRUE)
      }
      if (meth == "ordinal" && any(truth)) est_true <- mean(r$est[truth], na.rm = TRUE)
      rows[[length(rows) + 1]] <- data.frame(
        cond_id = cond$cond_id, rep = rep_id, method = meth, strategy = strat,
        fp = V, n_null = sum(!truth), tp = S, n_dif = sum(truth),
        fdp = if (R > 0) V / R else 0, any_fp = as.integer(V > 0),
        loc = loc, est_true = est_true, nbad = r$nbad, nfit = r$nfit,
        iters = r$iters, secs = secs)
    }
  }
  do.call(rbind, rows)
}

## ---------- Condition grid ----------------------------------------------
build_grid <- function() {
  g1 <- expand.grid(N = c(500, 1000, 2000), k = 0:2, node = c("direction", "middle", "extreme"),
                    impact = c("none", "trait+ERS"), stringsAsFactors = FALSE)
  g0 <- expand.grid(N = c(500, 1000, 2000), k = 0L, node = "none",
                    impact = c("none", "trait+ERS"), stringsAsFactors = FALSE)
  g <- rbind(g0, g1)
  g$shift <- 0.6
  g$imp_theta <- ifelse(g$impact == "none", 0, 0.5)
  g$imp_ers   <- ifelse(g$impact == "none", 0, 0.5)
  g$contam_pct <- g$k * 10
  g$cond_id <- seq_len(nrow(g)); g
}

run_sim <- function(grid, R, items, cores = 2L, seed = 20261014L) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  tasks <- expand.grid(rep = seq_len(R), cond_id = grid$cond_id)
  res <- mclapply(seq_len(nrow(tasks)), function(i) {
    cond <- grid[grid$cond_id == tasks$cond_id[i], ]
    tryCatch(one_rep(cond, items, tasks$rep[i]), error = function(e) NULL)
  }, mc.cores = cores, mc.set.seed = TRUE)
  do.call(rbind, res)
}

summarise_sim <- function(raw, grid) {
  key <- c("cond_id", "method", "strategy")
  agg <- function(f) aggregate(raw[, f, drop = FALSE], raw[key], function(v) mean(v, na.rm = TRUE))
  a <- aggregate(cbind(fp, n_null, tp, n_dif, fdp, any_fp, nbad, nfit, iters, secs) ~ cond_id + method + strategy,
                 raw, sum)
  nrep <- aggregate(rep ~ cond_id + method + strategy, raw, length); names(nrep)[4] <- "R"
  a <- merge(merge(a, nrep), agg(c("loc", "est_true")), all.x = TRUE)
  a$type1 <- a$fp / a$n_null
  a$power <- ifelse(a$n_dif > 0, a$tp / a$n_dif, NA)
  a$FDR <- a$fdp / a$R; a$FWER <- a$any_fp / a$R
  a$mcse_type1 <- sqrt(a$type1 * (1 - a$type1) / a$n_null)
  a$mcse_power <- ifelse(a$n_dif > 0, sqrt(a$power * (1 - a$power) / a$n_dif), NA)
  a$nonconv_rate <- a$nbad / a$nfit; a$mean_iters <- a$iters / a$R; a$sec_per_rep <- a$secs / a$R
  a <- merge(grid, a, by = "cond_id")
  a[order(a$cond_id, a$method, a$strategy),
    c("cond_id", "N", "contam_pct", "node", "impact", "method", "strategy", "R", "type1", "mcse_type1",
      "power", "mcse_power", "FDR", "FWER", "loc", "est_true", "nonconv_rate", "mean_iters", "sec_per_rep")]
}
