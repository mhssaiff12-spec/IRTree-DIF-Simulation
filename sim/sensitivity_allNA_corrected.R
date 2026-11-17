# Independent all-NA pseudo-item sensitivity. No original code/RDS modification.
# Run from sim/ with R_LIBS_USER pinned and one condition in IRTREE_COND.
suppressPackageStartupMessages(source("02b_latent_fast.R", local=TRUE))
source("sim_core.R")
Q <- 31L; J <- 20L; ALL <- seq_len(J)
cond <- Sys.getenv("IRTREE_COND")
stopifnot(grepl("^C0[2-5][0-9]$",cond))
meta <- read.csv("../input/replication_metadata.csv",stringsAsFactors=FALSE)
stopifnot(all(c("condition","replication","n_warnings_total") %in% names(meta)))
md <- meta[meta$condition==cond,]
stopifnot(nrow(md)==200L,!anyDuplicated(md$replication))
# reproducible selection, no overlaps: highest two, two random of remaining, one zero-warning
high <- md$replication[order(-md$n_warnings_total,md$replication)][1:2]
pool <- setdiff(md$replication,high)
set.seed(20261117L+as.integer(sub("C","",cond)))
random <- sample(pool,2L)
zero <- setdiff(md$replication[md$n_warnings_total==0],c(high,random))
stopifnot(length(zero)>0)
none <- sample(zero,1L)
selected <- data.frame(condition=cond,rep=c(high,random,none),
  selection=c("highest1","highest2","random1","random2","zero_warning"),
  n_warnings_total=md$n_warnings_total[match(c(high,random,none),md$replication)])
# Exact locked selection from the previously reviewed 140-row CSV.
locked <- read.csv("sensitivity_selected_140.csv", stringsAsFactors=FALSE)
stopifnot(nrow(locked)==140L, length(unique(paste(locked$condition,locked$rep)))==140L)
fixed <- locked[locked$condition==cond,,drop=FALSE]
stopifnot(nrow(fixed)==5L,identical(as.integer(fixed$rep),as.integer(selected$rep)),
          identical(as.character(fixed$selection),as.character(selected$selection)),
          identical(as.integer(fixed$n_warnings_total),as.integer(selected$n_warnings_total)))
selected <- fixed
pilot_rep <- Sys.getenv("IRTREE_PILOT_REP","")
if(nzchar(pilot_rep)) {
 stopifnot(cond=="C032",pilot_rep=="92",92L %in% selected$rep)
 selected <- selected[selected$rep==92L,,drop=FALSE]
}
dir.create("../sensitivity_out",showWarnings=FALSE)
write.csv(selected,sprintf("../sensitivity_out/%s_selection.csv",cond),row.names=FALSE)
if(Sys.getenv("IRTREE_SELECT_ONLY")=="1") quit(save="no")
grid <- read_grid("conditions_56.csv"); row<-grid[grid$condition_id==cond,,drop=FALSE]
stopifnot(nrow(row)==1L,as.integer(row$N)==1000L)
args <- cond_args(row)
params <- make_items(20L,777L)
key <- c(middle="M",direction="D",extreme="E")
fmt <- function(node) key[[node]]
get_original <- function(rep) {
 f<-sprintf("../input/%s/rep_%04d.rds",cond,rep)
 x<-readRDS(f)
 stopifnot(x$meta$manifest_id=="ab03dffdf580",x$meta$seed==data_seed(row,rep,"main"),
           nrow(x$results)==4L,nrow(x$nodes)==240L,nrow(x$links)==12L)
 x
}
ans<-list();decision_changes<-list();node_diagnostics<-list();errors<-list()
for(rep in selected$rep) {
 x<-get_original(rep)
 set.seed(data_seed(row,rep,"main"))
 dat<-sim_irtree(params,args$N,args$dif_items,args$dif_node,args$shift,
                 args$impact_theta,args$impact_ers)
 P<-pseudo_items(dat$Y,args$tree)
 groups<-factor(ifelse(dat$group==1L,"focal","reference"),levels=c("reference","focal"))
 for(method in c("Oracle1","RankK","RankK2","BH")) {
  arow<-x$anchors[x$anchors$method==method,,drop=FALSE]
  stopifnot(nrow(arow)==1L)
  anchors<-as.integer(strsplit(as.character(arow$anchors),";",fixed=TRUE)[[1]])
  orig<-x$nodes[x$nodes$method==method,,drop=FALSE]
  stopifnot(nrow(orig)==60L)
  out<-list();links<-list()
  for(nd in names(key)) {
   W<-P[[nd]];storage.mode(W)<-"numeric";colnames(W)<-paste0(key[[nd]],ALL)
   keep<-rowSums(!is.na(W))>0
   G<-droplevels(groups[keep]); WW<-W[keep,,drop=FALSE]
   cnt<-data.frame(condition=cond,rep=rep,method=method,tree=args$tree,
     dif_node=as.character(row$dif_node),shift=as.numeric(row$shift),impact=as.numeric(row$impact),
     node=key[[nd]],removed_reference=sum(!keep & groups=="reference"),
     removed_focal=sum(!keep & groups=="focal"),
     valid_reference=sum(keep & groups=="reference"),
     valid_focal=sum(keep & groups=="focal"))
   nw<-0L
   for(j in ALL) {
     anc<-paste0(key[[nd]],setdiff(anchors,j));cand<-paste0(key[[nd]],j)
     z<-tryCatch(withCallingHandlers({
       alt<-fit_node(WW,G,anc,Q);nul<-fit_node(WW,G,c(anc,cand),Q)
       conv<-isTRUE(mirt::extract.mirt(alt,"converged")) && isTRUE(mirt::extract.mirt(nul,"converged"))
       LR<-max(0,2*(as.numeric(logLik(alt))-as.numeric(logLik(nul))))
       df<-mirt::extract.mirt(alt,"nest")-mirt::extract.mirt(nul,"nest")
       data.frame(item=j,node=key[[nd]],p=pchisq(LR,df,lower.tail=FALSE),converged=conv,
                  error="",stringsAsFactors=FALSE)
     },warning=function(w){nw<<-nw+1L;invokeRestart("muffleWarning")}),
     error=function(e)data.frame(item=j,node=key[[nd]],p=NA_real_,
                          converged=FALSE,error=conditionMessage(e),stringsAsFactors=FALSE))
     out[[length(out)+1L]]<-z
   }
   link<-tryCatch(withCallingHandlers({
     m<-mirt::multipleGroup(WW,1,group=G,itemtype="Rasch",
       invariance=c(paste0(key[[nd]],anchors),"free_means","free_var"),
       quadpts=Q,verbose=FALSE,technical=list(NCYCLES=1000))
     fm <- as.numeric(unname(mirt::coef(m,simplify=TRUE)$focal$means))[1]
     cv <- isTRUE(mirt::extract.mirt(m,"converged"))
     list(value=fm,converged=cv,error="")
   },warning=function(w){nw<<-nw+1L;invokeRestart("muffleWarning")}),
     error=function(e)list(value=NA_real_,converged=FALSE,error=conditionMessage(e)))
   origlink<-x$links[x$links$method==method & x$links$node==nd,,drop=FALSE]
   stopifnot(nrow(origlink)==1L)
   cnt$original_link_mean<-as.numeric(origlink$est_mean[1])
   cnt$filtered_link_mean<-as.numeric(link$value)
   cnt$link_difference<-cnt$filtered_link_mean-cnt$original_link_mean
   cnt$link_converged<-isTRUE(link$converged)
   cnt$link_error<-link$error
   cnt$link_failed<-as.integer(!isTRUE(link$converged) || !is.finite(link$value))
   cnt$filtered_warning_count<-nw
   node_diagnostics[[length(node_diagnostics)+1L]]<-cnt
  }
  if(Sys.getenv("IRTREE_REFINED_CHECK")=="1" && cond=="C032" &&
      rep==92L && method=="RankK2") {
    fine<-list()
    for(nd in names(key)) {
      W<-P[[nd]]; storage.mode(W)<-"numeric"; colnames(W)<-paste0(key[[nd]],ALL)
      keep<-rowSums(!is.na(W))>0
      for(filtered in c(FALSE,TRUE)) {
        WD<-if(filtered) W[keep,,drop=FALSE] else W
        GD<-if(filtered) droplevels(groups[keep]) else groups
        for(j in ALL) {
          anc<-paste0(key[[nd]],setdiff(anchors,j)); cand<-paste0(key[[nd]],j)
          z<-tryCatch(suppressWarnings({
            alt<-mirt::multipleGroup(WD,1,group=GD,itemtype="Rasch",
              invariance=c(anc,"free_means","free_var"),method="EM",
              quadpts=Q,verbose=FALSE,TOL=1e-6,technical=list(NCYCLES=5000))
            nul<-mirt::multipleGroup(WD,1,group=GD,itemtype="Rasch",
              invariance=c(anc,cand,"free_means","free_var"),method="EM",
              quadpts=Q,verbose=FALSE,TOL=1e-6,technical=list(NCYCLES=5000))
            stopifnot(isTRUE(mirt::extract.mirt(alt,"converged")),
                      isTRUE(mirt::extract.mirt(nul,"converged")))
            LR<-max(0,2*(as.numeric(logLik(alt))-as.numeric(logLik(nul))))
            df<-mirt::extract.mirt(alt,"nest")-mirt::extract.mirt(nul,"nest")
            data.frame(item=j,node=key[[nd]],filtered=filtered,
                      p=pchisq(LR,df,lower.tail=FALSE),status="OK")
          }), error=function(e)data.frame(item=j,node=key[[nd]],filtered=filtered,
                                          p=NA_real_,status=conditionMessage(e)))
          fine[[length(fine)+1L]]<-z
        }
      }
    }
    fine<-do.call(rbind,fine)
    un<-fine[!fine$filtered,c("item","node","p","status")]
    fi<-fine[fine$filtered,c("item","node","p","status")]
    colnames(un)[3:4]<-c("p_unfiltered","status_unfiltered")
    colnames(fi)[3:4]<-c("p_filtered","status_filtered")
    pd<-merge(un,fi,by=c("item","node"),sort=FALSE)
    pd$p_abs_diff<-abs(pd$p_unfiltered-pd$p_filtered)
    stopifnot(nrow(pd)==60L)
    write.csv(pd,"../sensitivity_out/C032_rep92_RankK2_refined_p.csv",row.names=FALSE)
    cat("REFINED maximum |p_filtered-p_unfiltered| = ",
        max(pd$p_abs_diff,na.rm=TRUE),"\n")
  }
  result<-do.call(rbind,out)
  result$p_BH_all<-p.adjust(result$p,"BH")
  result$flag<-is.finite(result$p_BH_all) & result$p_BH_all < .05
  merged<-merge(orig[,c("item","node","p","flag")],result,by=c("item","node"),
               suffixes=c("_original","_filtered"),sort=FALSE)
  stopifnot(nrow(merged)==60L)
  merged$p_abs_diff<-abs(merged$p_original-merged$p_filtered)
  merged$changed_node<-merged$flag_original!=merged$flag_filtered
  origitem<-tapply(merged$flag_original,merged$item,any)
  newitem<-tapply(merged$flag_filtered,merged$item,any)
  item_same<-sum(origitem==newitem)
  v<-data.frame(condition=cond,rep=rep,method=method,tree=args$tree,
   dif_node=as.character(row$dif_node),shift=as.numeric(row$shift),
   max_p_abs_diff=if(any(is.finite(merged$p_abs_diff)))max(merged$p_abs_diff,na.rm=TRUE) else NA_real_,
   node_agreement=mean(!merged$changed_node),item_agreement=item_same/20,
   node_changes=sum(merged$changed_node),
   item_changes=sum(origitem!=newitem),
   nonconverged_filtered=sum(!result$converged),failed_tests=sum(!is.finite(result$p)))
  ans[[length(ans)+1L]]<-v
  if(Sys.getenv("IRTREE_REFINED_CHECK")=="1" && cond=="C032" && rep==92L && method=="RankK2") {
     write.csv(data.frame(condition=cond,rep=rep,method=method,
       max_p_abs_diff_original=v$max_p_abs_diff,
       max_p_abs_diff_refined=max(pd$p_abs_diff,na.rm=TRUE),
       nonfinite_refined=sum(!is.finite(pd$p_abs_diff))),
       "../sensitivity_out/C032_rep92_RankK2_refined_summary.csv",row.names=FALSE)
  }
  if(any(merged$changed_node) || item_same<20) {
    bad<-merged[merged$changed_node | merged$item %in% names(origitem)[origitem!=newitem],]
    bad$condition<-cond;bad$rep<-rep;bad$method<-method
    decision_changes[[length(decision_changes)+1L]]<-bad
  }
  if(any(!is.finite(result$p))||any(!result$converged))
    errors[[length(errors)+1L]]<-data.frame(condition=cond,rep=rep,method=method,
      note="nonfinite p or nonconvergence",n=sum(!is.finite(result$p)|!result$converged))
 }
 message(cond," rep=",rep," done")
}
summary<-do.call(rbind,ans);details<-do.call(rbind,node_diagnostics)
if(any(!is.finite(details$filtered_link_mean))) {
  stop("Filtered link means must be finite; inspect link_error before accepting results")
}
if(Sys.getenv("IRTREE_PILOT_REP")=="92") {
  stopifnot(nrow(details)==12L,all(details$link_converged),
            all(is.finite(details$filtered_link_mean)),all(is.finite(details$link_difference)))
}
write.csv(summary,sprintf("../sensitivity_out/%s_comparison.csv",cond),row.names=FALSE)
write.csv(details,sprintf("../sensitivity_out/%s_nodes.csv",cond),row.names=FALSE)
if(length(decision_changes)) {
 write.csv(do.call(rbind,decision_changes),
 sprintf("../sensitivity_out/%s_changes.csv",cond),row.names=FALSE)
} else {
 write.csv(data.frame(condition=character(),rep=integer(),method=character(),item=integer(),node=character()),
 sprintf("../sensitivity_out/%s_changes.csv",cond),row.names=FALSE)
}
if(length(errors)) write.csv(do.call(rbind,errors),sprintf("../sensitivity_out/%s_errors.csv",cond),row.names=FALSE)
cat(cond," completed 5 replications / 4 methods, node decision changes:",
 sum(summary$node_changes),"item decision changes:",sum(summary$item_changes),"\n")
