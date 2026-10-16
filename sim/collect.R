# =====================================================================
# collect.R — merges per-replication RDS files and computes performance measures with MCSE
# Usage: Rscript collect.R <results_dir> <output_dir>
# Per condition x method: R valid / failed; power; node power; Type I error (FPR on clean items);
# FWER = P(>=1 clean item flagged); FDR = E[FDP]; anchor quality; linking bias and RMSE per node
# (OK link estimates only; link failures counted SEPARATELY from DIF failures); linking shift vs
# Oracle; paired contrasts between methods and between trees (same data) on several outcomes.
# MCSE: mean -> sd/sqrt(R); RMSE -> sqrt(var(e^2) / (4 R MSE))  (Morris, White, & Crowther, 2019).
# =====================================================================
args <- commandArgs(TRUE); IN <- if (length(args)) args[1] else "results"; OUT <- if (length(args) > 1) args[2] else "summary"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fs <- list.files(IN, pattern = "^rep_[0-9]+\\.rds$", recursive = TRUE, full.names = TRUE)
if (!length(fs)) stop("no results in ", IN)
rbind_fill <- function(x) { x <- Filter(Negate(is.null), x); if (!length(x)) return(NULL); cols <- unique(unlist(lapply(x, names)))
  do.call(rbind, lapply(x, function(d) { for (k in setdiff(cols, names(d))) d[[k]] <- NA; d[, cols, drop = FALSE] })) }
L <- lapply(fs, readRDS); grab <- function(k) rbind_fill(lapply(L, `[[`, k))
meta <- grab("meta"); stopifnot(!anyDuplicated(meta[, c("condition_id", "rep")]))
if (!"fatal_error" %in% names(meta)) meta$fatal_error <- NA
grid <- read.csv(file.path(dirname(normalizePath(IN)), "conditions_56.csv"), fileEncoding = "UTF-8-BOM")
names(grid) <- sub("^X\\.+", "", names(grid))
res <- merge(grab("results"), grid[, c("condition_id", "N", "shift", "impact", "dif_node")])
lk <- grab("links"); nodes <- grab("nodes"); anc <- grab("anchors"); iss <- grab("issues")
for (nm in c("meta", "res", "lk", "nodes", "anc")) write.csv(get(nm), file.path(OUT, paste0("all_", nm, ".csv")), row.names = FALSE)
if (!is.null(iss)) write.csv(iss, file.path(OUT, "all_issues.csv"), row.names = FALSE)
undef <- res$tree == "EMD" & res$dif_node != "direction"; res$node_power[undef] <- NA
mcse <- function(x) { x <- x[is.finite(x)]; if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ }
mm <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
r4 <- function(d) { for (k in names(d)) if (is.numeric(d[[k]])) d[[k]] <- round(d[[k]], 4); d }

# 1) Replication accounting
V <- do.call(rbind, lapply(split(meta, meta$condition_id), function(z) data.frame(condition_id = z$condition_id[1], tree = z$tree[1],
  reps_saved = nrow(z), reps_fatal = sum(!is.na(z$fatal_error)), manifest_ids = paste(unique(z$manifest_id), collapse = ";"),
  mean_wall_min = mm(z$wall_seconds) / 60)))
write.csv(r4(V), file.path(OUT, "replication_accounting.csv"), row.names = FALSE)

# 2) Detection and anchor quality
S <- do.call(rbind, lapply(split(res, list(res$condition_id, res$method), drop = TRUE), function(z) { o <- z[z$status != "ERROR", ]
  data.frame(condition_id = z$condition_id[1], tree = z$tree[1], N = z$N[1], shift = z$shift[1], impact = z$impact[1], dif_node = z$dif_node[1],
    method = z$method[1], R_valid = nrow(o), R_failed_DIF = sum(z$status == "ERROR"),
    power = mm(o$power), power_mcse = mcse(o$power), node_power = mm(o$node_power), node_power_mcse = mcse(o$node_power),
    typeI_FPR = mm(o$FPR), typeI_FPR_mcse = mcse(o$FPR), FWER = mm(o$any_FP), FWER_mcse = mcse(o$any_FP),
    FDR = mm(o$FDP), FDR_mcse = mcse(o$FDP), anchor_size = mm(o$anchor_size), anchor_contaminated = mm(o$anchor_contaminated),
    anchor_contaminated_mcse = mcse(o$anchor_contaminated), P_any_contamination = mm(as.numeric(o$anchor_contaminated > 0)),
    P_any_contamination_mcse = mcse(as.numeric(o$anchor_contaminated > 0)), anchor_precision = mm(o$anchor_precision),
    sec_method = mm(o$seconds), sec_link = mm(o$seconds_link)) }))
write.csv(r4(S), file.path(OUT, "summary_detection_mcse.csv"), row.names = FALSE)

# 3) Linking: status accounting (separate from DIF), bias/RMSE vs known target, shift vs Oracle
lk$grp <- ifelse(grepl("^Oracle", lk$method), "Oracle", lk$method); lk$ok <- lk$link_status == "OK"
lk$e <- ifelse(lk$ok, lk$est_mean - lk$true_mean, NA)
LB <- do.call(rbind, lapply(split(lk, list(lk$condition_id, lk$method, lk$node), drop = TRUE), function(z) { e <- z$e[is.finite(z$e)]; R <- length(e)
  mse <- if (R) mean(e^2) else NA
  data.frame(condition_id = z$condition_id[1], tree = z$tree[1], method = z$method[1], node = z$node[1],
    link_OK = sum(z$link_status == "OK"), link_NONCONVERGED = sum(z$link_status == "NONCONVERGED"), link_ERROR = sum(z$link_status == "ERROR"),
    link_NOT_RUN_DIF_failed = sum(z$link_status == "NOT_RUN"), mean_est = mm(z$est_mean[z$ok]), target = z$true_mean[1],
    bias = if (R) mean(e) else NA, bias_mcse = mcse(e), RMSE = sqrt(mse), RMSE_mcse = if (R > 1 && mse > 0) sqrt(var(e^2) / (4 * R * mse)) else NA) }))
write.csv(r4(LB), file.path(OUT, "summary_linking_bias_rmse.csv"), row.names = FALSE)
om <- aggregate(est_mean ~ condition_id + rep + node, lk[lk$grp == "Oracle" & lk$ok, ], mean); names(om)[4] <- "oracle_est"
LO <- merge(lk[lk$grp != "Oracle" & lk$ok, ], om); LO$d <- LO$est_mean - LO$oracle_est
LOS <- do.call(rbind, lapply(split(LO, list(LO$condition_id, LO$method, LO$node), drop = TRUE), function(z)
  data.frame(condition_id = z$condition_id[1], tree = z$tree[1], method = z$method[1], node = z$node[1], R_pairs = nrow(z),
             shift_vs_oracle = mm(z$d), mcse = mcse(z$d), RMS_vs_oracle = sqrt(mm(z$d^2)))))
if (!is.null(LOS)) write.csv(r4(LOS), file.path(OUT, "summary_linking_vs_oracle.csv"), row.names = FALSE)

# 4) Paired contrasts between methods (same data; complete pairs only)
ok <- res[res$status != "ERROR", ]; ok$grp <- ifelse(grepl("^Oracle", ok$method), "Oracle", ok$method)
outs <- c("power", "FPR", "any_FP", "FDP", "anchor_contaminated")
con <- list(`RankK2-RankK` = c("RankK2", "RankK"), `RankK2-Oracle` = c("RankK2", "Oracle"), `RankK-Oracle` = c("RankK", "Oracle"),
            `BH-Oracle` = c("BH", "Oracle"), `RankK2-BH` = c("RankK2", "BH"))
PC <- do.call(rbind, lapply(outs, function(v) { x <- ok[is.finite(ok[[v]]), c("condition_id", "rep", "grp", v)]; if (!nrow(x)) return(NULL)
  x <- aggregate(x[[v]], x[, c("condition_id", "rep", "grp")], mean)
  W <- reshape(x, idvar = c("condition_id", "rep"), timevar = "grp", direction = "wide"); names(W) <- sub("^x\\.", "", names(W))
  do.call(rbind, lapply(split(W, W$condition_id), function(z) do.call(rbind, lapply(names(con), function(k) {
    if (!all(con[[k]] %in% names(z))) return(NULL); d <- z[[con[[k]][1]]] - z[[con[[k]][2]]]; d <- d[is.finite(d)]
    data.frame(outcome = v, condition_id = z$condition_id[1], contrast = k, R_pairs = length(d), mean_diff = mm(d), mcse = mcse(d),
               wins = sum(d > 0), ties = sum(d == 0), losses = sum(d < 0)) })))) }))
write.csv(r4(PC), file.path(OUT, "paired_method_contrasts.csv"), row.names = FALSE)

# 5) Paired tree contrasts EMD - MDE (same data seed, same method)
TE <- merge(ok[ok$tree == "EMD", c("seed", "rep", "method", outs)], ok[ok$tree == "MDE", c("seed", "rep", "method", "condition_id", outs)],
            by = c("seed", "rep", "method"), suffixes = c(".EMD", ".MDE"))
TC <- do.call(rbind, lapply(outs, function(v) do.call(rbind, lapply(split(TE, list(TE$condition_id, TE$method), drop = TRUE), function(z) {
  d <- z[[paste0(v, ".EMD")]] - z[[paste0(v, ".MDE")]]; d <- d[is.finite(d)]; if (!length(d)) return(NULL)
  data.frame(outcome = v, MDE_condition = z$condition_id[1], method = z$method[1], R_pairs = length(d), EMD_minus_MDE = mean(d), mcse = mcse(d),
             EMD_higher = sum(d > 0), same = sum(d == 0), EMD_lower = sum(d < 0)) }))))
if (!is.null(TC)) write.csv(r4(TC), file.path(OUT, "paired_tree_contrasts.csv"), row.names = FALSE)
cat("files:", length(fs), " conditions:", length(unique(meta$condition_id)), " fatal reps:", sum(!is.na(meta$fatal_error)),
    " DIF method errors:", sum(res$status == "ERROR"), " link not OK:", sum(lk$link_status != "OK"), "\n")
