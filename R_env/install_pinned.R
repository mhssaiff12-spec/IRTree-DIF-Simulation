# Installs the EXACT package versions used by the frozen simulation from the vendored source
# tarballs in R_env/vendor (no internet / CRAN needed). Dependency order is computed from DESCRIPTION.
# Usage: Rscript R_env/install_pinned.R <library_dir> [--with-recommended]
#   --with-recommended also installs lattice/Matrix/mgcv/MASS/nlme/cluster/codetools/class at the
#   R 4.3.3 versions (only needed if your R installation ships different versions).
args <- commandArgs(TRUE); lib <- if (length(args)) args[1] else Sys.getenv("R_LIBS_USER")
rec <- "--with-recommended" %in% args
dir.create(lib, recursive = TRUE, showWarnings = FALSE); .libPaths(c(lib, .libPaths()))
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
pins <- read.csv(file.path(here, "R_packages_pinned.csv"), stringsAsFactors = FALSE)
if (!rec) pins <- pins[pins$priority != "recommended", ]
tb <- file.path(here, "vendor", paste0(pins$package, "_", pins$version, ".tar.gz")); stopifnot(all(file.exists(tb)))
# read DESCRIPTION from each tarball
desc_of <- function(f, p) { td <- tempfile(); untar(f, files = file.path(p, "DESCRIPTION"), exdir = td); read.dcf(file.path(td, p, "DESCRIPTION")) }
need <- setNames(lapply(seq_len(nrow(pins)), function(i) { d <- desc_of(tb[i], pins$package[i])
  f <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(d)); x <- unlist(strsplit(paste(d[1, f], collapse = ","), ","))
  x <- trimws(sub("\\(.*", "", x)); intersect(x[nzchar(x)], pins$package) }), pins$package)
order <- character(); left <- pins$package
while (length(left)) { ready <- left[sapply(left, function(p) all(need[[p]] %in% order))]
  if (!length(ready)) stop("dependency cycle among: ", paste(left, collapse = ","))
  order <- c(order, ready); left <- setdiff(left, ready) }
for (p in order) {
  have <- tryCatch(as.character(packageVersion(p, lib.loc = lib)), error = function(e) NA)
  v <- pins$version[pins$package == p]
  if (!is.na(have) && have == as.character(package_version(v))) { message(p, " ", v, " already installed"); next }
  message("installing ", p, " ", v)
  install.packages(tb[pins$package == p], lib = lib, repos = NULL, type = "source", dependencies = FALSE, quiet = TRUE)
  if (!requireNamespace(p, lib.loc = lib, quietly = TRUE)) stop("failed to install ", p)
}
message("done: ", length(order), " packages in ", lib)
