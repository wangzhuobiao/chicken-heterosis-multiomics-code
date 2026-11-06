args <- commandArgs(trailingOnly=TRUE)
out <- normalizePath(args[1], winslash='/', mustWork=TRUE)
suppressPackageStartupMessages(library(DESeq2));library(jsonlite)
gate<-fromJSON(file.path(Sys.getenv('A1_INPUT_ROOT'),'09_output/Article1_R2_Stage4_20261029/RNA_internal_audit.json'));stopifnot(gate$status=='PASS_INTERNAL_ANALYSIS_INPUTS')
root<-Sys.getenv('A1_INPUT_ROOT');d<-file.path(root,'Article1_core_scripts_20260812','output_duodenum_transcriptome_core')
dest<-file.path(out,'F1_trial100');dir.create(dest,showWarnings=FALSE)
raw<-read.delim(Sys.getenv('A1_RNA_COUNTS'),check.names=FALSE,row.names=1,fileEncoding='UTF-8')
meta<-read.delim(file.path(root,'05_duodenum_RNAseq','sample_metadata.txt'),stringsAsFactors=FALSE)
cts<-as.matrix(raw[,meta$sample_id]);storage.mode(cts)<-'integer';rownames(meta)<-meta$sample_id;meta$group<-factor(meta$group,levels=c('XX','YY','YX'))
v<-read.csv(file.path(d,'Table3_VST_expression_matrix.csv'),check.names=FALSE,row.names=1);universe<-rownames(v)
cs<-list(YY_vs_XX=c('group','YY','XX'),YX_vs_YY=c('group','YX','YY'),YX_vs_XX=c('group','YX','XX'))
patterns<-c('not_tested','YX_high_transgressive','YX_low_transgressive','YY_dominant','XX_dominant','YX_intermediate','parental_similar_or_unresolved','unclassified')
tol<-c(0,.02,.05,.10)
classify<-function(means,res,t){
 a<-res[[1]]$pvalue<.05;b<-res[[2]]$pvalue<.05;c<-res[[3]]$pvalue<.05;a[is.na(a)]<-FALSE;b[is.na(b)]<-FALSE;c[is.na(c)]<-FALSE
 y<-means[,'YX'];lo<-pmin(means[,'YY'],means[,'XX']);hi<-pmax(means[,'YY'],means[,'XX'])
 p<-rep(7L,length(y));p[!a]<-6L;p[a&y>=lo-t&y<=hi+t&b&c]<-5L;p[a&b&!c]<-4L;p[a&!b&c]<-3L;p[y<lo-t&b&c]<-2L;p[y>hi+t&b&c]<-1L;p
}
baseline_res<-lapply(names(cs),function(n){r<-read.csv(file.path(d,paste0('DEG_',n,'.csv')),check.names=FALSE);r[match(universe,r$gene_id),]})
baseline_means<-sapply(c('YY','YX','XX'),function(g)rowMeans(v[,meta$sample_id[meta$group==g],drop=FALSE]))
base<-sapply(tol,function(t)classify(baseline_means,baseline_res,t));rownames(base)<-universe
state<-array(NA_integer_,c(length(universe),100,4),dimnames=list(universe,NULL,as.character(tol)))
set.seed(20261029);rng<-RNGkind();log<-data.frame();counts<-data.frame();trans<-data.frame();total<-proc.time()[3]
writeLines(c(capture.output(sessionInfo()),paste('seed',20261029),'B=100 attempts; tolerance=0,0.02,0.05,0.10; original nominal P<0.05; count>=5 in>=3; DESeq2 design~group; vst blind=FALSE'),file.path(dest,'session_and_settings.txt'))
for(i in seq_len(100)){
 start<-proc.time()[3];warnings<-character();sampled<-unlist(lapply(c('YY','YX','XX'),function(g)sample(which(meta$group==g),8,replace=TRUE)))
 ans<-tryCatch(withCallingHandlers({
  x<-cts[,sampled,drop=FALSE];colnames(x)<-paste0('b',i,'_',seq_len(24));mm<-meta[sampled,,drop=FALSE];rownames(mm)<-colnames(x)
  x<-x[rowSums(x>=5)>=3,,drop=FALSE];dds<-DESeqDataSetFromMatrix(x,mm,~group);dds<-DESeq(dds,quiet=TRUE)
  res<-lapply(cs,function(ct)as.data.frame(results(dds,contrast=ct,alpha=.05)))
  vm<-assay(vst(dds,blind=FALSE));means<-sapply(c('YY','YX','XX'),function(g)rowMeans(vm[,mm$group==g,drop=FALSE]))
  cls<-sapply(tol,function(t)classify(means,res,t));all<-matrix(0L,length(universe),4);j<-match(universe,rownames(means));ok<-!is.na(j);all[ok,]<-cls[j[ok],]
  saveRDS(list(attempt=i,sampled_original_ids=meta$sample_id[sampled],results=res,mean_VST=means,classes=all,warnings=warnings),file.path(dest,sprintf('replicate_%03d.rds',i)))
  all
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 elapsed<-proc.time()[3]-start
 if(inherits(ans,'error')){err<-conditionMessage(ans);saveRDS(list(attempt=i,sampled_ids=meta$sample_id[sampled],error=err,warnings=warnings),file.path(dest,sprintf('FAILED_%03d.rds',i)));status<-'FAILED'}else{
  err<-'';status<-'SUCCESS';state[,i,]<-ans
  for(k in seq_along(tol)){
   tab<-table(factor(ans[,k],levels=0:7));counts<-rbind(counts,data.frame(attempt=i,tolerance=tol[k],pattern=patterns,count=as.integer(tab)))
   tt<-as.data.frame(table(baseline=factor(base[,k],levels=0:7),bootstrap=factor(ans[,k],levels=0:7)));tt$attempt<-i;tt$tolerance<-tol[k];trans<-rbind(trans,tt)
  }
 }
 log<-rbind(log,data.frame(attempt=i,status=status,seconds=elapsed,error=err,warnings=paste(warnings,collapse=' | ')))
 write.table(log,file.path(dest,'attempt_log.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
 cat(sprintf('attempt %d/100 %s %.2fs\n',i,status,elapsed));flush.console()
}
success<-which(log$status=='SUCCESS');agreements<-data.frame()
for(k in seq_along(tol)){
 sm<-matrix(state[,success,k],nrow=length(universe));eval<-rowSums(sm!=0);agree<-rowSums(sm==base[,k]);originalagree<-rowSums(sm==base[,2])
 agreements<-rbind(agreements,data.frame(gene_id=universe,tolerance=tol[k],baseline_pattern=patterns[base[,k]+1],successful_replicates=length(success),evaluable_replicates=eval,
 agreement_all_successful=if(length(success))agree/length(success) else NA_real_,agreement_given_evaluable=ifelse(eval>0,agree/eval,NA_real_),agreement_with_original_tolerance002=if(length(success))originalagree/length(success) else NA_real_))
}
write.table(agreements,file.path(dest,'gene_classification_stability.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(counts,file.path(dest,'class_count_distribution.tsv'),sep='\t',row.names=FALSE,quote=FALSE);write.table(trans,file.path(dest,'classification_transitions.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
saveRDS(list(states=state,baseline=base,patterns=patterns,RNGkind=rng,seed=20261029),file.path(dest,'all_classification_states.rds'))
write_json(list(attempts=100,successes=length(success),failures=sum(log$status=='FAILED'),failure_rate=mean(log$status=='FAILED'),wall_seconds=unname(proc.time()[3]-total),automatic_500=FALSE),file.path(dest,'summary.json'),auto_unbox=TRUE,pretty=TRUE)
