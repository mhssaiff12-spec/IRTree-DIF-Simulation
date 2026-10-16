# Verifies that data generation matches each condition (item 2 of the review).
# Part 1 (all 56 rows): mapping, sample size, 50/50 groups, seed uniqueness, MDE/EMD share data.
# Part 2 (14 distinct generating cells, large N = 40000): per node, Rasch multiple-group fit on the
#   MDE pseudo-items with the 16 clean items (5–20) as anchors; recovers (a) the intercept shift of
#   items 1–4 (expected = shift on dif_node, 0 elsewhere) and (b) the focal latent mean
#   (expected: middle 0, direction = impact, extreme = impact).
suppressMessages(source("02b_latent_fast.R")); source("sim_core.R")
g <- read_grid("conditions_56.csv"); validate_grid(g); p <- make_items(20, 777)
fp <- function(Y) sum(Y * outer(seq_len(nrow(Y)), seq_len(ncol(Y)), function(i, j) ((i * 31 + j * 17) %% 97)))
gen <- function(r, seed, N = NULL) { a <- cond_args(r); set.seed(seed)
  sim_irtree(p, if (is.null(N)) a$N else N, a$dif_items, a$dif_node, a$shift, a$impact_theta, a$impact_ers) }
P1 <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) { r <- g[i, ]; s <- data_seed(r, 1); d <- gen(r, s)
  data.frame(condition_id = r$condition_id, N = r$N, tree = r$tree, cell = cell_code(r), seed_rep1 = s,
             n_rows = nrow(d$Y), n_focal = sum(d$group == 1), categories_1to5 = all(d$Y %in% 1:5), fingerprint = fp(d$Y)) }))
cat("Part 1 rows:", nrow(P1), " N matches:", all(P1$n_rows == P1$N), " groups 50/50:", all(P1$n_focal == P1$N / 2),
    " categories ok:", all(P1$categories_1to5), "\n")
cat("distinct cells:", length(unique(P1$cell)), " (expected 28); each cell has exactly MDE+EMD:",
    all(table(P1$cell) == 2), "; MDE/EMD fingerprints equal within cell:",
    all(tapply(P1$fingerprint, P1$cell, function(x) length(unique(x)) == 1)), "\n")
seeds <- unlist(lapply(seq_len(nrow(g)), function(i) data_seed(g[i, ], 1:1000)))
cat("data seeds unique across cells for reps 1..1000 (MDE/EMD share):", !anyDuplicated(unique(seeds)) && length(unique(seeds)) == 28 * 1000, "\n")
k <- read_grid("conditions_core8.csv"); cat("core8 legacy cases:", paste(k$condition_id, sapply(seq_len(nrow(k)), function(i) legacy_case(k[i, ])), collapse = ", "), "\n")
write.csv(P1, "verify_part1_conditions.csv", row.names = FALSE)

# Part 2
cells <- g[!duplicated(g[, c("shift", "impact", "dif_node")]), ]; NN <- 40000L
key <- c(middle = "M", direction = "D", extreme = "E")
P2 <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) { r <- cells[i, ]; d <- gen(r, 4242L + i, NN)
  P <- pseudo_items(d$Y, "MDE"); grp <- factor(ifelse(d$group == 1, "focal", "reference"), levels = c("reference", "focal"))
  do.call(rbind, lapply(names(key), function(nd) {
    W <- P[[nd]]; storage.mode(W) <- "numeric"; colnames(W) <- paste0(key[nd], 1:20)
    m <- suppressWarnings(multipleGroup(W, 1, group = grp, itemtype = "Rasch", invariance = c(paste0(key[nd], 5:20), "free_means", "free_var"),
                                        quadpts = 31, verbose = FALSE, technical = list(NCYCLES = 2000)))
    cf <- coef(m, simplify = TRUE); dd <- cf$focal$items[1:4, "d"] - cf$reference$items[1:4, "d"]
    exp_shift <- if (r$dif_node == nd) r$shift else 0
    data.frame(shift = r$shift, impact = r$impact, dif_node = r$dif_node, node = nd,
               est_shift_items1to4 = round(mean(dd), 3), expected_shift = exp_shift,
               est_focal_mean = round(cf$focal$means[1], 3), expected_mean = c(middle = 0, direction = r$impact, extreme = r$impact)[[nd]],
               converged = extract.mirt(m, "converged")) })) }))
P2$shift_err <- round(P2$est_shift_items1to4 - P2$expected_shift, 3); P2$mean_err <- round(P2$est_focal_mean - P2$expected_mean, 3)
print(P2, row.names = FALSE)
cat("max |shift error|:", max(abs(P2$shift_err)), "  max |mean error|:", max(abs(P2$mean_err)), "  all converged:", all(P2$converged), "\n")
write.csv(P2, "verify_part2_dif_and_impact_recovery.csv", row.names = FALSE)
