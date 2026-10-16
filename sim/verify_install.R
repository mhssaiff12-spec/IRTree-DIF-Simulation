# Checks a new machine reproduces the reference replication exactly (C047, legacy seed, rep 1 = stage 4H case B rep 1).
# Run:  IRTREE_COND=C047 IRTREE_REP_FROM=1 IRTREE_REP_TO=1 IRTREE_SEED_SCHEME=legacy IRTREE_ORACLE_DRAWS=3 IRTREE_OUT=check Rscript main_sim.R
#       Rscript verify_install.R
a <- readRDS("reference_C047_legacy_rep1.rds"); b <- readRDS("check/C047/rep_0001.rds")
m <- merge(a$nodes[, c("method", "item", "node", "p")], b$nodes[, c("method", "item", "node", "p")], by = c("method", "item", "node"))
cat("node tests:", nrow(m), " max |p diff|:", max(abs(m$p.x - m$p.y)), "\n")
cat("anchors identical:", identical(a$anchors$anchors, b$anchors$anchors), "  metrics identical:",
    isTRUE(all.equal(a$results[, c("TP", "FP", "anchor_contaminated")], b$results[, c("TP", "FP", "anchor_contaminated")])), "\n")
cat(if (nrow(m) == 360 && max(abs(m$p.x - m$p.y)) < 1e-6) "PASS\n" else "CHECK VERSIONS (R, mirt)\n")
