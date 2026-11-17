# Diagnostic for one frozen replication only: C032 rep 92 RankK2.
# Compare A=original saved RDS, B=fresh unfiltered, C=fresh all-NA-removed.
# Never writes to RDS, main_sim.R, frozen source files, or original results branches.
suppressPackageStartupMessages(source("02b_latent_fast.R", local=TRUE))
source("sim_core.R")
suppressPackageStartupMessages(library(mirt))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L,args %in% c("session1","session2"))
session <- args[1]
q <- 31L; J <- 20L; K <- 10L
grid <- read_grid("conditions_56.csv")
row <- grid[grid$condition_id=="C032",,drop=FALSE]
stopifnot(nrow(row)==1L,row$tree=="EMD",row$N==1000L)
ca <- cond_args(row)
seed <- data_seed(row,92L,"main")
orig <- readRDS("../input/C032/rep_0092.rds")
stopifnot(orig$meta$manifest_id=="ab03dffdf580",
          orig$meta$rep==92L,orig$meta$seed==seed)
arow <- orig$anchors[orig$anchors$method=="RankK2",,drop=FALSE]
stopifnot(nrow(arow)==1L)
anchors <- as.integer(strsplit(as.character(arow$anchors[1]),";",fixed=TRUE)[[1]])
stopifnot(length(anchors)==10L,anyDuplicated(anchors)==0L)
dir.create("../diagnostic_out",showWarnings=FALSE)
write.csv(data.frame(condition="C032",rep=92,method="RankK2",session=session,
  seed=seed,manifest_id=as.character(orig$meta$manifest_id[1]),
  tree=as.character(row$tree[1]),anchors=paste(anchors,collapse=";"),
  original_anchors=paste(arow$anchors[1]),stringsAsFactors=FALSE),
  sprintf("../diagnostic_out/metadata_%s.csv",session),row.names=FALSE)

params <- make_items(J,seed=777L)
set.seed(seed)
dat <- sim_irtree(params,ca$N,ca$dif_items,ca$dif_node,ca$shift,ca$impact_theta,ca$impact_ers)
P <- pseudo_items(dat$Y,ca$tree)
group <- factor(ifelse(dat$group==1,"focal","reference"),levels=c("reference","focal"))
keys <- c(middle="M",direction="D",extreme="E")
rows <- list()
do_lr <- function(W,G,anc,cand) {
  alt <- fit_node(W,G,anc,q)
  nul <- fit_node(W,G,c(anc,cand),q)
  converged <- isTRUE(mirt::extract.mirt(alt,"converged")) && isTRUE(mirt::extract.mirt(nul,"converged"))
  LR <- max(0,2*(as.numeric(logLik(alt))-as.numeric(logLik(nul))))
  df <- mirt::extract.mirt(alt,"nest")-mirt::extract.mirt(nul,"nest")
  c(p=pchisq(LR,df,lower.tail=FALSE),LR=LR,df=df,converged=as.numeric(converged))
}
for(nd in names(keys)) {
  W<-P[[nd]]; storage.mode(W)<-"numeric";colnames(W)<-paste0(keys[[nd]],seq_len(J))
  keep<-rowSums(!is.na(W))>0L
  n_removed<-sum(!keep)
  for(j in seq_len(J)) {
    anc<-paste0(keys[[nd]],setdiff(anchors,j)); cand<-paste0(keys[[nd]],j)
    z<-tryCatch(suppressWarnings({
      b<-do_lr(W,group,anc,cand)
      if(session=="session1") {
        c<-do_lr(W[keep,,drop=FALSE],droplevels(group[keep]),anc,cand)
      } else c<-rep(NA_real_,4L)
      data.frame(item=j,node=keys[[nd]],p_B=unname(b["p"]),p_C=unname(c["p"]),
        LR_B=unname(b["LR"]),LR_C=unname(c["LR"]),
        converged_B=as.logical(b["converged"]),converged_C=as.logical(c["converged"]),
        removed=n_removed,error="",stringsAsFactors=FALSE)
    }),error=function(e)data.frame(item=j,node=keys[[nd]],p_B=NA_real_,p_C=NA_real_,
      LR_B=NA_real_,LR_C=NA_real_,converged_B=FALSE,converged_C=FALSE,
      removed=n_removed,error=conditionMessage(e),stringsAsFactors=FALSE))
    rows[[length(rows)+1L]]<-z
  }
  message(session," completed node ",nd)
}
out<-do.call(rbind,rows)
stopifnot(nrow(out)==60L,all(is.finite(out$p_B)),all(out$converged_B))
if(session=="session1") stopifnot(all(is.finite(out$p_C)),all(out$converged_C))
write.csv(out,sprintf("../diagnostic_out/%s_raw_p.csv",session),row.names=FALSE)
cat("anchors:",paste(anchors,collapse=","),"\n")
cat("session ",session," finished 60 tests\n")
