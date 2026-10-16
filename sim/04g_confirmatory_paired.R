# Stage 4G: Paired confirmatory pilot: RankK vs RankK2 + multiple randomized OracleK references
# Methods: randomized Oracle references, FixedK, BH purification, RankK and RankK2.
# Requires mirt. NOT validated in R in the authoring environment.
if(!requireNamespace('mirt',quietly=TRUE)) stop('Install mirt: install.packages("mirt")')
source('02b_latent_fast.R',local=TRUE)
getint <- function(n,default,minimum) {
 x <- suppressWarnings(as.integer(Sys.getenv(n,as.character(default))))
 if(is.na(x)||x<minimum) stop('Invalid ',n);x
}
REPS <- getint('IRTREE_REPS',1L,1L)
ORACLE_DRAWS <- getint('IRTREE_ORACLE_DRAWS',3L,1L)
N <- getint('IRTREE_N',1000L,200L)
Q <- getint('IRTREE_QUAD',31L,15L)
K <- getint('IRTREE_K',10L,6L)
MAXITER <- getint('IRTREE_ITER',5L,1L)
J <- 20L; ALL <- seq_len(J); ALPHA <- .05
if(K>=J) stop('IRTREE_K must be <20')
OUT <- 'stage4g_output';dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
params <- make_items(J,seed=777L)
scenarios <- list(A=list(shift=0,it=0,ie=0),B=list(shift=.6,it=0,ie=0),
                  C=list(shift=0,it=.5,ie=.5),D=list(shift=.6,it=.5,ie=.5))
truth <- function(cc) if(scenarios[[cc]]$shift>0) 1:4 else integer()
fit_all <- function(dat,tree,anchors,cc) {
 z <- run_case_fast(cc,dat,test_items=ALL,pure_anchors=anchors,
             dif_node='extreme',shift=scenarios[[cc]]$shift,tree=tree,quadpts=Q)
 if(nrow(z)!=3L*J||any(!is.finite(z$p))||any(z$df!=1L)||any(!z$converged))
   stop('Invalid/nonconverged node LR test')
 z$p_BH_all <- p.adjust(z$p,method='BH')
 z$flag <- z$p_BH_all<ALPHA
 z
}
# Ranking based on minimum unadjusted node p for each item, used ONLY to rank
# anchor candidates. High values are considered less suspicious. The
# final detection always uses BH across all 60 node tests.
rank_k <- function(fit, k) {
 s <- tapply(fit$p,fit$item,min)
 score <- as.numeric(s[as.character(ALL)])
 if(any(!is.finite(score))) stop('Non-finite ranking score')
 sort(ALL[order(-score,ALL)[seq_len(k)]])
}
# BH sequential purification comparator; no knowledge of actual DIF.
bh_purify <- function(dat,tree,cc,initial=NULL) {
 anc <- ALL; hist <- list(); status <- 'ITER_LIMIT'
 for(i in seq_len(MAXITER+1L)) {
   f <- if(i==1L && !is.null(initial)) initial else fit_all(dat,tree,anc,cc)
   bad <- intersect(anc,unique(f$item[f$flag]))
   new <- setdiff(anc,bad)
   st <- if(!length(bad)) 'STABLE' else if(length(new)<K) 'MIN_K' else if(i>MAXITER) 'ITER_LIMIT' else 'CONTINUE'
   hist[[length(hist)+1L]] <- data.frame(step=i,anchors_before=paste(anc,collapse=';'),
     anchors_after=paste(if(st=='CONTINUE')new else anc,collapse=';'),
     removed=if(st=='CONTINUE')paste(bad,collapse=';') else '',status=st)
   if(st!='CONTINUE') {status<-st;break}
   anc<-new
 }
 list(anchor=anc,fit=f,history=do.call(rbind,hist),status=status)
}
metrics <- function(fit,anc,tru) {
 called <- unique(fit$item[fit$flag]);tp<-length(intersect(called,tru));fp<-length(setdiff(called,tru));ncall<-length(called)
 data.frame(TP=tp,FP=fp,FN=length(tru)-tp,
  FPR=fp/(J-length(tru)),power=if(length(tru))tp/length(tru) else NA_real_,
  FDP=if(ncall)fp/ncall else 0,any_FP=as.integer(fp>0),
  anchor_size=length(anc),anchor_contaminated=length(intersect(anc,tru)),
  anchor_precision=length(setdiff(anc,tru))/length(anc),
  clean_items_removed=length(setdiff(setdiff(ALL,tru),anc)),
  contaminated_items_excluded=length(setdiff(tru,anc)))
}
rows<-list();nodes<-list();selection<-list();failures<-list()
append_select<-function(cc,rr,tree,method,step,anc,detail) {
 selection[[length(selection)+1L]] <<- data.frame(case=cc,rep=rr,tree=tree,method=method,
 step=step,anchor_size=length(anc),anchors=paste(anc,collapse=';'),detail=detail)
}
for(cc in names(scenarios)) for(rr in seq_len(REPS)) {
 seed<-24000000L+match(cc,names(scenarios))*100000L+rr
 set.seed(seed);s<-scenarios[[cc]]
 dat<-sim_irtree(params,N,dif_items=1:4,dif_node='extreme',shift=s$shift,
                 impact_theta=s$it,impact_ers=s$ie)
 tru<-truth(cc)
 for(tree in c('MDE','EMD')) {
   # Cache all-other screening: no duplicate computation across RankK variants.
   screen<-tryCatch(fit_all(dat,tree,ALL,cc),error=function(e)e)
   for(method in c(paste0('Oracle',seq_len(ORACLE_DRAWS)),'FixedK','RankK','RankK2','BH')) {
      t0<-proc.time()[['elapsed']];warnings<-character()
      result<-tryCatch(withCallingHandlers({
        if(inherits(screen,'error')) stop(conditionMessage(screen))
        if(startsWith(method,'Oracle')) {
          draw <- as.integer(sub('Oracle','',method))
          set.seed(seed + 10000L * match(tree,c('MDE','EMD')) + draw)
          anc<-sort(sample(setdiff(ALL,tru),K,replace=FALSE))
          append_select(cc,rr,tree,method,0L,anc,'Oracle only; sampled from genuinely DIF-free items')
          list(anchor=anc,fit=fit_all(dat,tree,anc,cc),status='ORACLE')
        } else if(method=='FixedK') {
          anc<-tail(ALL,K)
          append_select(cc,rr,tree,method,0L,anc,'Prespecified, no truth used')
          list(anchor=anc,fit=fit_all(dat,tree,anc,cc),status='FIXED')
        } else if(method=='RankK') {
          anc<-rank_k(screen,K)
          append_select(cc,rr,tree,method,1L,anc,'Ranked once by minimum node p')
          list(anchor=anc,fit=fit_all(dat,tree,anc,cc),status='RANKED')
        } else if(method=='RankK2') {
          anc1<-rank_k(screen,K)
          append_select(cc,rr,tree,method,1L,anc1,'First screening')
          re<-fit_all(dat,tree,anc1,cc)
          anc2<-rank_k(re,K)
          append_select(cc,rr,tree,method,2L,anc2,'Second screening; fixed K')
          list(anchor=anc2,fit=fit_all(dat,tree,anc2,cc),status=if(identical(anc1,anc2))'STABLE' else 'CHANGED')
        } else {
          b<-bh_purify(dat,tree,cc,initial=screen)
          for(ii in seq_len(nrow(b$history))) {
            bh<-b$history[ii,]
            append_select(cc,rr,tree,method,bh$step,
              as.integer(strsplit(if(bh$status=='CONTINUE')bh$anchors_after else bh$anchors_before,';',fixed=TRUE)[[1]]),bh$status)
          }
          list(anchor=b$anchor,fit=b$fit,status=b$status)
        }
      },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
      error=function(e)list(error=conditionMessage(e)))
      info<-data.frame(case=cc,rep=rr,seed=seed,tree=tree,method=method,
           seconds=proc.time()[['elapsed']]-t0,warnings=length(warnings),
           status=if(is.null(result$error))result$status else 'ERROR',
           error=if(is.null(result$error))'' else result$error)
      if(!is.null(result$error)) {
        failures[[length(failures)+1L]]<-info
        rows[[length(rows)+1L]]<-cbind(info,as.data.frame(as.list(setNames(rep(NA_real_,11),c('TP','FP','FN','FPR','power','FDP','any_FP','anchor_size','anchor_contaminated','anchor_precision','clean_items_removed')))),contaminated_items_excluded=NA_real_)
      } else {
        rows[[length(rows)+1L]]<-cbind(info,metrics(result$fit,result$anchor,tru))
        z<-result$fit;z$rep<-rr;z$seed<-seed;z$method<-method
        z$node_truth_interpretable<-tree=='MDE'
        nodes[[length(nodes)+1L]]<-z
      }
      message(cc,' rep=',rr,' ',tree,' ',method,' ',info$status)
   }
 }
 write.csv(do.call(rbind,rows),file.path(OUT,'replication_results.csv'),row.names=FALSE)
 if(length(nodes))write.csv(do.call(rbind,nodes),file.path(OUT,'node_results.csv'),row.names=FALSE)
 if(length(selection))write.csv(do.call(rbind,selection),file.path(OUT,'anchor_selection_history.csv'),row.names=FALSE)
 if(length(failures))write.csv(do.call(rbind,failures),file.path(OUT,'failures.csv'),row.names=FALSE)
}
r<-do.call(rbind,rows)
fields<-c('FPR','power','FDP','any_FP','anchor_size','anchor_contaminated','anchor_precision','clean_items_removed','contaminated_items_excluded')
sm<-do.call(rbind,lapply(split(r,interaction(r$case,r$tree,r$method,drop=TRUE)),function(z) {
 do.call(rbind,lapply(fields,function(v) {
 x<-if(v %in% names(z))z[[v]] else numeric()
 x<-x[is.finite(x)];m<-length(x)
 data.frame(case=z$case[1],tree=z$tree[1],method=z$method[1],metric=v,
 attempted=nrow(z),valid=m,mean=if(m)mean(x) else NA_real_,
 MCSE=if(m>1)sd(x)/sqrt(m) else NA_real_)
 }))
}))
write.csv(sm,file.path(OUT,'summary_mcse.csv'),row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),
 paste('REPS',REPS,'N',N,'Q',Q,'K',K,'MAXITER',MAXITER),
 'Truth is passed only to randomized Oracle references and evaluations.',
 'RankK and RankK2 use min raw p to screen; no claim of calibrated post-selection inference.',
 'Final BH is across 60 node LR tests, independent of selection method.',
 'FixedK is clean only in current generator, not a realistic unknown-clean reference.',
 paste('ORACLE_DRAWS',ORACLE_DRAWS),
 'Equal anchor length does not ensure identical information/quality.',
 'EMD node labels not equivalent to MDE node labels; item-level outcomes only.',
 'Latent dimensions assumed independent; this is a feasibility analysis.',
 'Software not executed in R in preparation environment.'),file.path(OUT,'session_info.txt'))
print(sm)

# Per-replication PAIRED method differences. A replication is the independent unit.
paired <- list()
for (cc in names(scenarios)) for (tt in c('MDE','EMD')) {
  for (rr in seq_len(REPS)) {
    z <- r[r$case==cc & r$tree==tt & r$rep==rr,]
    a <- z[z$method=='RankK',]; b <- z[z$method=='RankK2',]
    if(nrow(a)!=1L||nrow(b)!=1L) next
    for(v in c('power','FPR','FDP','any_FP','anchor_contaminated','clean_items_removed')) {
      paired[[length(paired)+1L]] <- data.frame(case=cc,tree=tt,rep=rr,metric=v,
        RankK=a[[v]],RankK2=b[[v]],difference=b[[v]]-a[[v]],
        valid=is.finite(a[[v]]) && is.finite(b[[v]]))
    }
  }
}
if(length(paired)) {
  pd<-do.call(rbind,paired)
  write.csv(pd,file.path(OUT,'paired_differences.csv'),row.names=FALSE)
  ps<-do.call(rbind,lapply(split(pd,interaction(pd$case,pd$tree,pd$metric,drop=TRUE)),function(d) {
    x<-d$difference[d$valid];nn<-length(x)
    data.frame(case=d$case[1],tree=d$tree[1],metric=d$metric[1],n=nn,
      mean_difference=if(nn)mean(x) else NA_real_,
      MCSE=if(nn>1)sd(x)/sqrt(nn) else NA_real_)
  }))
  write.csv(ps,file.path(OUT,'paired_summary.csv'),row.names=FALSE)
}
# Within-data Oracle variability due exclusively to its different valid anchor subsets.
ora<-r[grepl('^Oracle[0-9]+$',r$method),]
if(nrow(ora)) {
  oo<-do.call(rbind,lapply(split(ora,interaction(ora$case,ora$tree,ora$rep,drop=TRUE)),function(z){
    data.frame(case=z$case[1],tree=z$tree[1],rep=z$rep[1],
       valid_oracles=sum(is.finite(z$FPR)),
       range_power=if(sum(is.finite(z$power))>1)diff(range(z$power,na.rm=TRUE)) else NA_real_,
       range_FPR=if(sum(is.finite(z$FPR))>1)diff(range(z$FPR,na.rm=TRUE)) else NA_real_)
  }))
  write.csv(oo,file.path(OUT,'oracle_variability.csv'),row.names=FALSE)
}
cat('Stage 4G complete. Outputs in ',OUT,'\n')
