# =====================================================================
# sim_core.R — shared definitions for the main simulation (no estimation code here)
#   read_grid()  : reads the condition grid (handles the UTF-8 BOM present in both CSVs)
#   validate_grid(): structural checks (completeness, duplicates, admissible values)
#   cond_args()  : maps one grid row to the arguments of sim_irtree()
#   cell_code()  : data-generation cell = N x shift x impact x dif_node (tree EXCLUDED,
#                  so MDE and EMD analyse the SAME data -> paired tree comparison)
#   data_seed()  : 'main' scheme (new, independent of all development stages) or
#                  'legacy' scheme (reproduces stages 4C–4J for the 4 cells A–D; verification only)
# =====================================================================
read_grid <- function(path) {
  g <- read.csv(path, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
  names(g) <- sub("^X\\.+", "", names(g)); g
}
validate_grid <- function(g) {
  need <- c("condition_id", "N", "shift", "impact", "tree", "dif_node", "K", "dif_items", "status")
  stopifnot(all(need %in% names(g)))
  stopifnot(!anyDuplicated(g$condition_id), !anyDuplicated(g[, c("N", "shift", "impact", "tree", "dif_node", "K", "dif_items")]))
  stopifnot(all(g$tree %in% c("MDE", "EMD")), all(g$dif_node %in% c("none", "middle", "direction", "extreme")))
  stopifnot(all((g$shift == 0) == (g$dif_node == "none")))       # null <=> no DIF node
  stopifnot(all(g$N %in% c(500, 1000, 2000)), all(g$shift %in% c(0, .3, .6)), all(g$impact %in% c(0, .5)))
  invisible(TRUE)
}
parse_items <- function(s) as.integer(strsplit(as.character(s), ";", fixed = TRUE)[[1]])
cond_args <- function(row) list(
  N = as.integer(row$N), shift = as.numeric(row$shift),
  impact_theta = as.numeric(row$impact), impact_ers = as.numeric(row$impact),  # 'impact' applies to theta AND ERS (as in stages 1–4J)
  dif_node = if (row$dif_node == "none") "extreme" else row$dif_node,       # irrelevant when shift = 0 (adds exactly 0)
  dif_items = parse_items(row$dif_items), tree = row$tree, K = as.integer(row$K),
  truth = if (row$shift > 0) parse_items(row$dif_items) else integer())
cell_code <- function(row) {
  match(row$N, c(500, 1000, 2000)) * 100L + match(row$shift, c(0, .3, .6)) * 10L - 10L +
    match(row$impact, c(0, .5)) * 4L - 4L + match(row$dif_node, c("none", "middle", "direction", "extreme")) - 1L
}
legacy_case <- function(row) {   # stages 4C–4J cases (N = 1000, extreme node, shift 0/.6, impact 0/.5)
  if (row$N != 1000 || !(row$dif_node %in% c("none", "extreme")) || row$shift == .3) return(NA_character_)
  c("A", "B", "C", "D")[1L + (row$shift > 0) + 2L * (row$impact > 0)]
}
data_seed <- function(row, rep, scheme = c("main", "legacy")) {
  scheme <- match.arg(scheme)
  if (scheme == "legacy") {
    cc <- legacy_case(row); if (is.na(cc)) stop("legacy seeds exist only for the 4 cells A–D")
    return(24000000L + match(cc, c("A", "B", "C", "D")) * 100000L + as.integer(rep))
  }
  stopifnot(rep >= 1, rep <= 9999)
  100000000L + cell_code(row) * 100000L + as.integer(rep)
}
# Oracle draw seed: legacy = 4H/4I formula (tree-dependent); main = tree-INDEPENDENT, so the
# MDE and EMD analyses of the same data use the same random clean anchor set (cleaner pairing).
oracle_seed <- function(dseed, tree, draw, scheme)
  if (scheme == "legacy") dseed + 10000L * match(tree, c("MDE", "EMD")) + draw else dseed + 10000L + draw
