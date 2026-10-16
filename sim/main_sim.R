# =====================================================================
# main_sim.R — main simulation driver (one process = one list of condition x replication tasks)
#   * Estimation/detection algorithms are read VERBATIM from 04g_confirmatory_paired.R (lines 23–66)
#     and 02b_latent_fast.R; nothing in them is changed.
#   * Methods per replication (all on the SAME generated data): OracleK (random clean K; benchmark),
#     RankK, RankK2, BH (sequential BH purification starting from the all-other screen).
#   * Exact reuse of identical fits (no change in results, verified against 4H):
#       - RankK2's first re-fit uses RankK's anchors -> identical to RankK's final fit (reused);
#       - if RankK2's re-selected anchors equal RankK's, its final fit is that same fit (reused);
#       - BH starts from the screen fit (as in 4H).
#   * Truth labels for node results are recomputed here from the condition (fit_all() passes
#     dif_node = 'extreme' to run_case_fast(), which only affects those label columns).
#   * One RDS file per condition x replication, written atomically (tmp + rename) => resumable:
#     tasks whose file exists are skipped; an interrupted task leaves no partial file.
# Environment variables:
#   IRTREE_GRID (conditions_56.csv) IRTREE_COND (comma list or 'all') IRTREE_REP_FROM IRTREE_REP_TO
#   IRTREE_SEED_SCHEME (main|legacy) IRTREE_ORACLE_DRAWS (1) IRTREE_OUT (results) IRTREE_REUSE (1)
#   IRTREE_WORKER / IRTREE_NWORKERS (split the task list across independent processes)
# =====================================================================
suppressMessages(source("02b_latent_fast.R", local = TRUE)); source("sim_core.R")
env <- function(n, d) Sys.getenv(n, d)
geti <- function(n, d) { x <- suppressWarnings(as.integer(env(n, d))); if (is.na(x)) stop("Invalid ", n); x }
GRID <- read_grid(env("IRTREE_GRID", "conditions_56.csv")); validate_grid(GRID)
CONDS <- env("IRTREE_COND", "all"); CONDS <- if (CONDS == "all") GRID$condition_id else strsplit(CONDS, ",")[[1]]
stopifnot(all(CONDS %in% GRID$condition_id))
REP_FROM <- geti("IRTREE_REP_FROM", "1"); REP_TO <- geti("IRTREE_REP_TO", "1")
SCHEME <- env("IRTREE_SEED_SCHEME", "main"); ORACLE_DRAWS <- geti("IRTREE_ORACLE_DRAWS", "1")
REUSE <- geti("IRTREE_REUSE", "1") == 1L; OUT <- env("IRTREE_OUT", "results")
WORKER <- geti("IRTREE_WORKER", "1"); NWORKERS <- geti("IRTREE_NWORKERS", "1")
Q <- 31L; MAXITER <- 5L; J <- 20L; ALL <- seq_len(J); ALPHA <- .05
params <- make_items(J, seed = 777L)
src4g <- readLines("04g_confirmatory_paired.R")
stopifnot(grepl("^fit_all <- function", src4g[23]), grepl("^metrics <- function", src4g[57]))
eval(parse(text = src4g[23:66]))            # fit_all, rank_k, bh_purify, metrics (K is read from the global K)

with_preserved_rng <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv()); old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()))
  set.seed(seed); force(expr)
}
link_means <- function(dat, tree, anchors) {
  # Linking estimates are recorded SEPARATELY from DIF testing: each node gets its own status
  # (OK / NONCONVERGED / ERROR), error text and warning count; nothing here affects DIF decisions.
  P <- pseudo_items(dat$Y, tree); g <- factor(ifelse(dat$group == 1, "focal", "reference"), levels = c("reference", "focal"))
  key <- c(middle = "M", direction = "D", extreme = "E")
  do.call(rbind, lapply(names(key), function(nd) {
    W <- P[[nd]]; storage.mode(W) <- "numeric"; colnames(W) <- paste0(key[nd], ALL)
    nw <- 0L; wt <- character()
    est <- tryCatch(withCallingHandlers({
      m <- mirt::multipleGroup(W, 1, group = g, itemtype = "Rasch", invariance = c(paste0(key[nd], anchors), "free_means", "free_var"),
                               quadpts = Q, verbose = FALSE, technical = list(NCYCLES = 1000))
      list(v = mirt::coef(m, simplify = TRUE)$focal$means[1], conv = isTRUE(mirt::extract.mirt(m, "converged")), err = "") },
      warning = function(w) { nw <<- nw + 1L; wt <<- c(wt, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) list(v = NA_real_, conv = FALSE, err = conditionMessage(e)))
    st <- if (nzchar(est$err)) "ERROR" else if (!est$conv || !is.finite(est$v)) "NONCONVERGED" else "OK"
    data.frame(node = nd, est_mean = est$v, link_converged = est$conv, link_status = st, link_error = est$err,
               link_n_warnings = nw, link_warning_text = paste(unique(wt), collapse = " | ")) }))
}
link_not_run <- function(reason) data.frame(node = c("middle", "direction", "extreme"), est_mean = NA_real_, link_converged = FALSE,
  link_status = "NOT_RUN", link_error = reason, link_n_warnings = 0L, link_warning_text = "")
# true focal mean per node: defined under MDE; under EMD only the direction pseudo-items are
# identical to MDE's (Y==3 -> NA, else Y>=4), so only that node has a known target there.
true_means <- function(tree, impact) {
  v <- c(middle = 0, direction = impact, extreme = impact); if (tree == "EMD") v[c("middle", "extreme")] <- NA; v }
timed <- function(expr, we) { t0 <- proc.time()[["elapsed"]]
  v <- withCallingHandlers(expr, warning = function(w) { we$w <- c(we$w, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(val = v, sec = proc.time()[["elapsed"]] - t0) }

run_one <- function(row, rep) {
  a <- cond_args(row); K <<- a$K; scenarios <<- list(X = list(shift = a$shift))
  seed <- data_seed(row, rep, SCHEME); set.seed(seed)
  dat <- sim_irtree(params, a$N, a$dif_items, a$dif_node, a$shift, a$impact_theta, a$impact_ers)
  tru <- a$truth; tree <- a$tree; we <- new.env(); we$w <- character(); T0 <- proc.time()[["elapsed"]]
  key <- c(middle = "M", direction = "D", extreme = "E")
  fix_truth <- function(z) { z$true_dif_item <- z$item %in% tru
    # node labels are comparable to the generating node only under MDE, or for the direction node under EMD
    # (identical pseudo-items); otherwise node-level truth is undefined -> NA
    lab_ok <- row$dif_node != "none" && (tree == "MDE" || row$dif_node == "direction")
    z$true_dif_node <- if (lab_ok) z$true_dif_item & z$node == key[a$dif_node] else if (length(tru)) NA else FALSE; z }
  rows <- list(); nodes <- list(); sel <- list(); links <- list(); issues <- list(); cache <- list()
  sc <- tryCatch(timed(fit_all(dat, tree, ALL, "X"), we), error = function(e) e)
  methods <- c(paste0("Oracle", seq_len(ORACLE_DRAWS)), "RankK", "RankK2", "BH")
  for (m in methods) {
    nw0 <- length(we$w)
    res <- tryCatch({
      if (inherits(sc, "error")) stop("screen failed: ", conditionMessage(sc))
      if (startsWith(m, "Oracle")) {
        d <- as.integer(sub("Oracle", "", m))
        anc <- with_preserved_rng(oracle_seed(seed, tree, d, SCHEME), sort(sample(setdiff(ALL, tru), a$K, replace = FALSE)))
        f <- timed(fit_all(dat, tree, anc, "X"), we); list(anchor = anc, fit = f$val, sec = f$sec, status = "ORACLE", hist = paste(anc, collapse = ";"))
      } else if (m == "RankK") {
        anc <- rank_k(sc$val, a$K); f <- timed(fit_all(dat, tree, anc, "X"), we); cache$A1 <- anc; cache$F1 <- f$val; cache$s1 <- f$sec
        list(anchor = anc, fit = f$val, sec = f$sec, status = "RANKED", hist = paste(anc, collapse = ";"))
      } else if (m == "RankK2") {
        A1 <- rank_k(sc$val, a$K)
        F1 <- if (REUSE && identical(A1, cache$A1)) list(val = cache$F1, sec = 0) else timed(fit_all(dat, tree, A1, "X"), we)
        A2 <- rank_k(F1$val, a$K)
        F2 <- if (REUSE && identical(A2, A1)) list(val = F1$val, sec = 0) else timed(fit_all(dat, tree, A2, "X"), we)
        list(anchor = A2, fit = F2$val, sec = F1$sec + F2$sec, status = if (identical(A1, A2)) "STABLE" else "CHANGED",
             hist = paste(paste(A1, collapse = ";"), paste(A2, collapse = ";"), sep = " -> "))
      } else {
        t0 <- proc.time()[["elapsed"]]
        b <- withCallingHandlers(bh_purify(dat, tree, "X", initial = sc$val),
               warning = function(w) { we$w <- c(we$w, conditionMessage(w)); invokeRestart("muffleWarning") })
        list(anchor = b$anchor, fit = b$fit, sec = proc.time()[["elapsed"]] - t0, status = b$status,
             hist = paste(b$history$anchors_after, collapse = " -> "))
      }
    }, error = function(e) list(error = conditionMessage(e)))
    nw <- length(we$w) - nw0
    if (!is.null(res$error)) {
      issues[[length(issues) + 1]] <- data.frame(method = m, type = "ERROR", n = 1L, text = res$error)
      rows[[length(rows) + 1]] <- data.frame(method = m, status = "ERROR", seconds = NA, n_warnings = nw)
      L <- link_not_run(paste("DIF method failed:", res$error)); L$method <- m; L$true_mean <- true_means(tree, row$impact)[L$node]
      links[[length(links) + 1]] <- L; next }
    if (nw) issues[[length(issues) + 1]] <- data.frame(method = m, type = "WARNING", n = nw, text = paste(unique(tail(we$w, nw)), collapse = " | "))
    z <- fix_truth(res$fit)
    node_tp <- if (length(tru)) sum(z$flag & z$true_dif_node) else NA
    mt <- metrics(res$fit, res$anchor, tru)
    lm <- timed(link_means(dat, tree, res$anchor), we)
    rows[[length(rows) + 1]] <- cbind(data.frame(method = m, status = res$status, seconds = res$sec, seconds_link = lm$sec, n_warnings = nw, link_ok_nodes = sum(lm$val$link_status == "OK")), mt,
      data.frame(node_power = if (length(tru)) node_tp / length(tru) else NA,
                 FP_node_tests = sum(z$flag & !z$true_dif_item), flags_on_dif_items_wrong_node = sum(z$flag & z$true_dif_item & !z$true_dif_node)))
    z$method <- m; nodes[[length(nodes) + 1]] <- z
    sel[[length(sel) + 1]] <- data.frame(method = m, anchors = paste(res$anchor, collapse = ";"), history = res$hist,
                                         contaminated = paste(intersect(res$anchor, tru), collapse = ";"))
    L <- lm$val; L$method <- m; L$true_mean <- true_means(tree, row$impact)[L$node]; links[[length(links) + 1]] <- L
  }
  meta <- data.frame(condition_id = row$condition_id, rep = rep, seed = seed, seed_scheme = SCHEME, N = a$N, shift = a$shift,
                     impact = row$impact, tree = tree, dif_node = row$dif_node, K = a$K,
                     screen_seconds = if (inherits(sc, "error")) NA else sc$sec, wall_seconds = proc.time()[["elapsed"]] - T0,
                     n_warnings_total = length(we$w), host = Sys.info()[["nodename"]], time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  add <- function(x) if (length(x)) { b <- rbind_fill(x); cbind(meta[rep(1, nrow(b)), c("condition_id", "rep", "seed", "tree")], b) }
  list(meta = meta, results = add(rows), nodes = add(nodes), anchors = add(sel), links = add(links), issues = add(issues))
}

# rbind for data frames with different columns (an ERROR row has fewer columns than a successful one)
rbind_fill <- function(x) { cols <- unique(unlist(lapply(x, names)))
  do.call(rbind, lapply(x, function(d) { for (k in setdiff(cols, names(d))) d[[k]] <- NA; d[, cols, drop = FALSE] })) }

# ---- Frozen-version guard: code hashes, R and package versions must match FROZEN_MANIFEST.csv ----
FROZEN_FILES <- c("main_sim.R", "sim_core.R", "02b_latent_fast.R", "irtree_dif_sim.R", "04g_confirmatory_paired.R")
PKGS <- c("mirt", "GPArotation", "Deriv", "dcurver", "Rcpp", "RcppArmadillo", "Matrix", "mgcv", "vegan", "SimDesign", "pbapply", "lattice")
current_manifest <- function() rbind(
  data.frame(kind = "file_md5", name = FROZEN_FILES, value = unname(tools::md5sum(FROZEN_FILES))),
  data.frame(kind = "R", name = "R.version", value = R.version.string),
  data.frame(kind = "RNG", name = "RNGkind", value = paste(RNGkind(), collapse = "/")),
  data.frame(kind = "package", name = PKGS, value = sapply(PKGS, function(p) as.character(packageVersion(p)))))
if (Sys.getenv("IRTREE_WRITE_MANIFEST") == "1") { write.csv(current_manifest(), "FROZEN_MANIFEST.csv", row.names = FALSE); message("manifest written"); quit(save = "no") }
if (Sys.getenv("IRTREE_REQUIRE_FROZEN", "1") == "1") {
  if (!file.exists("FROZEN_MANIFEST.csv")) stop("FROZEN_MANIFEST.csv missing")
  fm <- read.csv("FROZEN_MANIFEST.csv", stringsAsFactors = FALSE); cm <- current_manifest()
  mm <- merge(fm, cm, by = c("kind", "name"), all = TRUE); bad <- mm[is.na(mm$value.x) | is.na(mm$value.y) | mm$value.x != mm$value.y, ]
  if (nrow(bad)) { print(bad); stop("code/package versions differ from FROZEN_MANIFEST.csv") }
}
MANIFEST_ID <- substr(unname(tools::md5sum("FROZEN_MANIFEST.csv")), 1, 12)
if (Sys.getenv("IRTREE_NO_RUN") == "1") stop("IRTREE_NO_RUN set: functions loaded, no tasks run")

# Task order is replication-major (all conditions advance together); MDE/EMD of a cell are adjacent.
tasks <- expand.grid(condition_id = CONDS, rep = REP_FROM:REP_TO, stringsAsFactors = FALSE)
tasks <- tasks[seq_len(nrow(tasks)) %% NWORKERS == (WORKER - 1L) %% NWORKERS, , drop = FALSE]
BUDGET <- as.numeric(env("IRTREE_TIME_BUDGET_MIN", "Inf")) * 60; EXPECT <- as.numeric(env("IRTREE_EXPECT_REP_MIN", "11")) * 60
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); START <- proc.time()[["elapsed"]]
for (i in seq_len(nrow(tasks))) {
  id <- tasks$condition_id[i]; rp <- tasks$rep[i]; d <- file.path(OUT, id); dir.create(d, showWarnings = FALSE)
  f <- file.path(d, sprintf("rep_%04d.rds", rp))
  if (file.exists(f)) next
  if (proc.time()[["elapsed"]] - START + EXPECT > BUDGET) { message("time budget reached; stopping cleanly before ", id, " rep ", rp); break }
  res <- tryCatch(run_one(GRID[GRID$condition_id == id, ], rp), error = function(e) {
    row <- GRID[GRID$condition_id == id, ]
    list(meta = data.frame(condition_id = id, rep = rp, seed = data_seed(row, rp, SCHEME), seed_scheme = SCHEME, tree = row$tree,
                           fatal_error = conditionMessage(e)), results = NULL, nodes = NULL, anchors = NULL, links = NULL, issues = NULL) })
  res$meta$manifest_id <- MANIFEST_ID
  tmp <- paste0(f, ".tmp_", Sys.getpid()); saveRDS(res, tmp); file.rename(tmp, f)
  message(sprintf("%s rep %d done: %.0fs wall, method errors=%d", id, rp, if (is.null(res$meta$wall_seconds)) NA else res$meta$wall_seconds,
                  if (is.null(res$results)) NA else sum(res$results$status == "ERROR")))
}
writeLines(c(capture.output(sessionInfo()), paste("manifest:", MANIFEST_ID)), file.path(OUT, sprintf("session_info_worker%d.txt", WORKER)))
