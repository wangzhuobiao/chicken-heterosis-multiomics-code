# Audit saved results only: no DESeq(), vst(), results(), or new statistical fit.
args <- commandArgs(trailingOnly=TRUE)
out <- normalizePath(args[1],winslash='/',mustWork=TRUE)
root <- Sys.getenv('A1_INPUT_ROOT'); old <- file.path(root,'09_output/Article1_R2_Stage4_20261029')
src <- Sys.getenv('A1_F1_RUN'); core <- file.path(root,'Article1_core_scripts_20260812/output_duodenum_transcriptome_core')
options(digits=17,scipen=30)
write_tsv <- function(x,name)write.table(x,file.path(out,name),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
obj <- readRDS(file.path(src,'all_classification_states.rds'))
st <- obj$states; bas <- obj$baseline; pats <- obj$patterns; tol <- as.numeric(dimnames(st)[[3]])
genes <- dimnames(st)[[1]]; N <- length(genes); B <- dim(st)[2]
stopifnot(N==19429,B==100,identical(tol,c(0,.02,.05,.1)),!anyNA(st),!anyDuplicated(genes))
prior <- read.csv(file.path(core,'F1_inheritance_standalone2/Table_F1_expression_pattern_all_genes.csv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
prior <- prior[match(genes,prior$gene_id),]; stopifnot(!anyNA(prior$gene_id),all(pats[bas[,2]+1]==prior$pattern))
sel <- bas[,2] %in% 1:5;stopifnot(sum(sel)==1534)
subsets <- list(all_19429=rep(TRUE,N),original_1534=sel)
meta <- read.delim(file.path(root,'05_duodenum_RNAseq/sample_metadata.txt'),stringsAsFactors=FALSE)
raw_path <- Sys.getenv('A1_RNA_COUNTS')
raw <- read.delim(raw_path,check.names=FALSE,row.names=1,fileEncoding='UTF-8')
cts <- as.matrix(raw[,meta$sample_id]); storage.mode(cts)<-'integer'
v <- read.csv(file.path(core,'Table3_VST_expression_matrix.csv'),check.names=FALSE,row.names=1)
v <- v[genes,meta$sample_id]
cts_names <- c('YY_vs_XX','YX_vs_YY','YX_vs_XX')
base_res <- lapply(cts_names,function(n){r<-read.csv(file.path(core,paste0('DEG_',n,'.csv')),row.names=1,check.names=FALSE);r[genes,]})
base_means <- sapply(c('YY','YX','XX'),function(g)rowMeans(v[,meta$sample_id[meta$group==g],drop=FALSE]))
# Independent explicit first-match classifier preserves original rule priority.
classify <- function(means,res,t){
 a<-res[[1]]$pvalue<.05;b<-res[[2]]$pvalue<.05;c<-res[[3]]$pvalue<.05
 a[is.na(a)]<-FALSE;b[is.na(b)]<-FALSE;c[is.na(c)]<-FALSE
 y<-means[,'YX'];lo<-pmin(means[,'YY'],means[,'XX']);hi<-pmax(means[,'YY'],means[,'XX'])
 conditions<-cbind(y>hi+t&b&c,y<lo-t&b&c,a&!b&c,a&b&!c,a&y>=lo-t&y<=hi+t&b&c,!a,TRUE)
 apply(conditions,1,function(r)which(r)[1])
}
baseline_recomputed <- sapply(tol,function(t)classify(base_means,base_res,t))
stopifnot(all(baseline_recomputed==bas))
write_tsv(data.frame(gene_id=genes,symbol=prior$symbol,original_pattern=prior$pattern,original_1534=sel,baseline_t0=pats[bas[,1]+1],baseline_t002=pats[bas[,2]+1],baseline_t005=pats[bas[,3]+1],baseline_t010=pats[bas[,4]+1]),'baseline_gene_classes.tsv')
bc<-do.call(rbind,lapply(seq_along(tol),function(k)data.frame(tolerance=tol[k],pattern=pats,count=as.integer(table(factor(bas[,k],levels=0:7))))))
write_tsv(bc,'baseline_class_counts.tsv')
checks<-list();dups<-list();counts<-list();changes<-list();change_genes<-list();evidence<-list();pdiag<-list();sampling<-list()
pairs<-combn(1:4,2);set.seed(obj$seed)
for(i in 1:B){
 r<-readRDS(file.path(src,sprintf('replicate_%03d.rds',i)))
 ids<-r$sampled_original_ids;replay<-unlist(lapply(c('YY','YX','XX'),function(g)sample(which(meta$group==g),8,replace=TRUE)))
 mm<-meta[match(ids,meta$sample_id),];xx<-cts[,ids,drop=FALSE];colnames(xx)<-paste0('b',i,'_',seq_len(24))
 keep<-rowSums(xx>=5)>=3;expected<-rownames(xx)[keep];rid<-rownames(r$mean_VST)
 j<-match(genes,rid);ok<-!is.na(j);rec<-matrix(0L,N,4)
 local<-sapply(tol,function(t)classify(r$mean_VST,r$results,t));rec[ok,]<-local[j[ok],]
 checks[[i]]<-data.frame(attempt=i,draws=length(ids),unique_individuals=length(unique(ids)),unique_instance_names=length(unique(colnames(xx))),YY_draws=sum(mm$group=='YY'),YX_draws=sum(mm$group=='YX'),XX_draws=sum(mm$group=='XX'),rng_replay_identical=identical(ids,meta$sample_id[replay]),filtered_genes=length(expected),filter_matches_saved_ID_order=identical(expected,rid),all_contrast_gene_order_matches=all(vapply(r$results,function(a)identical(rownames(a),rid),logical(1))),recomputed_matches_replicate=all(rec==r$classes),replicate_matches_aggregate=all(r$classes==st[,i,]),warnings=paste(r$warnings,collapse=' | '))
 for(g in c('YY','YX','XX'))dups[[length(dups)+1]]<-data.frame(attempt=i,group=g,draws=sum(mm$group==g),unique_individuals=length(unique(ids[mm$group==g])),duplicate_draws=sum(mm$group==g)-length(unique(ids[mm$group==g])))
 sampling[[i]]<-data.frame(attempt=i,instance_name=colnames(xx),original_id=ids,group=mm$group)
 for(q in 1:3){a<-r$results[[q]];pdiag[[length(pdiag)+1]]<-data.frame(attempt=i,contrast=names(r$results)[q],N_genes=nrow(a),NA_pvalue=sum(is.na(a$pvalue)),NA_padj=sum(is.na(a$padj)),NA_log2FC=sum(is.na(a$log2FoldChange)),NA_stat=sum(is.na(a$stat)),NA_lfcSE=sum(is.na(a$lfcSE)))}
 for(s in names(subsets)){
   use<-subsets[[s]]
   for(k in 1:4)counts[[length(counts)+1]]<-data.frame(attempt=i,tolerance=tol[k],subset=s,pattern=pats,count=as.integer(table(factor(st[use,i,k],levels=0:7))))
   for(p in 1:ncol(pairs)){
     a<-pairs[1,p];b<-pairs[2,p];changed<-which(use&st[,i,a]!=st[,i,b])
     changes[[length(changes)+1]]<-data.frame(attempt=i,subset=s,N_genes=sum(use),tolerance_a=tol[a],tolerance_b=tol[b],changed_genes=length(changed))
     if(length(changed))change_genes[[length(change_genes)+1]]<-data.frame(attempt=i,subset=s,gene_id=genes[changed],symbol=prior$symbol[changed],original_pattern=prior$pattern[changed],tolerance_a=tol[a],tolerance_b=tol[b],class_a=pats[st[changed,i,a]+1],class_b=pats[st[changed,i,b]+1])
   }
 }
 # Gene-level evidence for t=0 vs .10 differences and boundary eligibility.
 changed<-which(rec[,1]!=rec[,4]); a<-r$results[[1]]$pvalue;b<-r$results[[2]]$pvalue;c<-r$results[[3]]$pvalue
 if(length(changed)){
  jj<-j[changed];m<-r$mean_VST[jj,,drop=FALSE]
  evidence[[length(evidence)+1]]<-data.frame(attempt=i,gene_id=genes[changed],original_1534=sel[changed],original_pattern=prior$pattern[changed],mean_YY=m[,'YY'],mean_YX=m[,'YX'],mean_XX=m[,'XX'],p_YY_XX=a[jj],p_YX_YY=b[jj],p_YX_XX=c[jj],upper_margin=m[,'YX']-pmax(m[,'YY'],m[,'XX']),lower_margin=pmin(m[,'YY'],m[,'XX'])-m[,'YX'],class_t0=pats[rec[changed,1]+1],class_t010=pats[rec[changed,4]+1])
 }
 if(i%%10==0){cat('audited',i,'/100\n');flush.console()}
}
checks<-do.call(rbind,checks);write_tsv(checks,'per_replicate_method_checks.tsv')
stopifnot(all(checks$recomputed_matches_replicate),all(checks$replicate_matches_aggregate),all(checks$filter_matches_saved_ID_order),all(checks$all_contrast_gene_order_matches),all(checks$draws==24),all(checks$unique_instance_names==24))
write_tsv(do.call(rbind,sampling),'sample_draw_instances.tsv');write_tsv(do.call(rbind,dups),'duplicate_sampling_checks.tsv');write_tsv(do.call(rbind,pdiag),'contrast_result_NA_diagnostics.tsv')
counts<-do.call(rbind,counts);changes<-do.call(rbind,changes)
write_tsv(counts,'per_replicate_class_counts.tsv');write_tsv(changes,'tolerance_changes_per_replicate.tsv')
write_tsv(do.call(rbind,change_genes),'tolerance_changed_gene_instances.tsv');write_tsv(do.call(rbind,evidence),'tolerance_boundary_gene_evidence.tsv')
change_summary<-do.call(rbind,lapply(split(changes,interaction(changes$subset,changes$tolerance_a,changes$tolerance_b,drop=TRUE)),function(x)data.frame(subset=x$subset[1],tolerance_a=x$tolerance_a[1],tolerance_b=x$tolerance_b[1],N_genes=x$N_genes[1],attempts=nrow(x),total_gene_attempt_changes=sum(x$changed_genes),replicates_with_changes=sum(x$changed_genes>0),min=min(x$changed_genes),mean=mean(x$changed_genes),median=median(x$changed_genes),max=max(x$changed_genes))))
write_tsv(change_summary,'02_tolerance_changes_summary.tsv')
summ<-list();detail<-list();matrices<-list();overall<-list()
for(k in 1:4){
 sm<-matrix(st[,,k],nrow=N);agree<-rowMeans(sm==bas[,2]);same<-rowMeans(sm==bas[,k]);evaluable<-rowSums(sm!=0)
 detail[[k]]<-data.frame(gene_id=genes,symbol=prior$symbol,tolerance=tol[k],original_pattern=prior$pattern,baseline_at_tolerance=pats[bas[,k]+1],original_1534=sel,retention_original=agree,agreement_same_tolerance=same,evaluable_replicates=evaluable)
 for(code in 0:7)detail[[k]][[paste0('to_',pats[code+1])]]<-rowSums(sm==code)
 for(s in names(subsets)){
  use<-subsets[[s]]
  overall[[length(overall)+1]]<-data.frame(tolerance=tol[k],subset=s,N_genes=sum(use),mean_original_agreement=mean(agree[use]),median_original_agreement=median(agree[use]),mean_same_tolerance_agreement=mean(same[use]),median_same_tolerance_agreement=median(same[use]))
  tt<-table(baseline=factor(rep(bas[use,2],B),levels=0:7),bootstrap=factor(as.vector(sm[use,,drop=FALSE]),levels=0:7));long<-as.data.frame(tt);long$baseline<-pats[as.integer(as.character(long$baseline))+1];long$bootstrap<-pats[as.integer(as.character(long$bootstrap))+1];long$tolerance<-tol[k];long$subset<-s;long$reference<-'frozen baseline tolerance=0.02'
  denominators<-table(factor(bas[use,2],levels=0:7))*B;long$row_probability<-ifelse(denominators[match(long$baseline,pats)]>0,long$Freq/as.numeric(denominators[match(long$baseline,pats)]),NA_real_);matrices[[length(matrices)+1]]<-long
  for(code in 0:7){
   ix<-which(use&bas[,2]==code);cc<-counts[counts$tolerance==tol[k]&counts$subset==s&counts$pattern==pats[code+1],'count']
   summ[[length(summ)+1]]<-data.frame(tolerance=tol[k],subset=s,original_pattern=pats[code+1],original_N=length(ix),mean_retention=if(length(ix))mean(agree[ix]) else NA_real_,median_retention=if(length(ix))median(agree[ix]) else NA_real_,Q25_retention=if(length(ix))unname(quantile(agree[ix],.25)) else NA_real_,Q75_retention=if(length(ix))unname(quantile(agree[ix],.75)) else NA_real_,bootstrap_count_min=min(cc),bootstrap_count_Q025=unname(quantile(cc,.025)),bootstrap_count_mean=mean(cc),bootstrap_count_median=median(cc),bootstrap_count_Q975=unname(quantile(cc,.975)),bootstrap_count_max=max(cc))
  }
 }
}
summ<-do.call(rbind,summ);detail<-do.call(rbind,detail);matrices<-do.call(rbind,matrices)
write_tsv(summ,'01_all_patterns_stability_summary.tsv');write_tsv(summ[summ$original_pattern%in%pats[2:6],],'01_five_patterns_stability_summary.tsv')
write_tsv(do.call(rbind,overall),'overall_stability_high_precision.tsv')
write_tsv(matrices,'03_class_transition_matrices_long.tsv');write_tsv(detail[detail$original_1534,],'04_original1534_gene_stability_details.tsv');write_tsv(detail,'all19429_gene_stability_details.tsv')
for(k in 1:4)for(s in names(subsets)){
 x<-matrices[matrices$tolerance==tol[k]&matrices$subset==s,];m<-matrix(x$Freq,8,8,dimnames=list(pats,pats));write_tsv(data.frame(original_pattern=rownames(m),m,check.names=FALSE),paste0('matrix_t',tol[k],'_',s,'.tsv'))
}
saveRDS(list(seed=obj$seed,RNGkind=obj$RNGkind,tolerances=tol,independent_classifier=classify),file.path(out,'audit_settings.rds'))
writeLines(c(capture.output(sessionInfo()),paste('source seed',obj$seed),paste('source RNG',paste(obj$RNGkind,collapse=', ')),'No DESeq2 fitting, VST recomputation, permutation or expanded bootstrap performed.'),file.path(out,'audit_session.txt'))
cat('COMPLETE. RNG replay matches:',sum(checks$rng_replay_identical),'/100\n')
