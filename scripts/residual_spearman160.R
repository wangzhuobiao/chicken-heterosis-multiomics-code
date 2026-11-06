O<-normalizePath(commandArgs(trailingOnly=TRUE)[1],winslash='/',mustWork=TRUE);R<-Sys.getenv('A1_INPUT_ROOT');options(stringsAsFactors=FALSE,digits=17,scipen=30)
rd<-function(p)read.csv(p,check.names=FALSE,fileEncoding='UTF-8-BOM');ts<-function(x,n)write.table(x,file.path(O,n),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
C<-file.path(R,'Article1_multiomics/Article1_multiomics_scripts/output_multiomics_phenotype_anchored_article1_core');S<-Sys.getenv('A1_ASSOCIATIONS')
script<-file.path(Sys.getenv('A1_CODE'),'frozen_definitions.R')
expr<-parse(script,encoding='UTF-8');for(e in expr)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])[1]%in%c('safe_numeric','residualize','safe_cor_test','outcome_specs'))eval(e)
MIN_COMPLETE_N<-12
ph<-rd(file.path(C,'Table_core_phenotypes_for_integration.csv'));ph$group<-factor(ph$group,levels=c('YY','YX','XX'));old<-rd(file.path(C,'Table_integrated_module_scores_long.csv'));new<-read.delim(file.path(S,'C_modules/all20_updated_scores.tsv'),check.names=FALSE,fileEncoding='UTF-8-BOM');ref<-rd(file.path(C,'Table_module_phenotype_residual_spearman.csv'))
stopifnot(nrow(ph)==24,nrow(old)==480,nrow(new)==480,nrow(ref)==160,length(outcome_specs)==8,!anyDuplicated(ph$sample_id))
key<-function(d)paste(d$module_id,d$outcome,sep='|')
calculate<-function(scores){out<-list();i<-1
 for(mid in unique(old$module_id)){
  a<-scores[scores$module_id==mid,];stopifnot(nrow(a)==24,!anyDuplicated(a$sample_id),setequal(a$sample_id,ph$sample_id));a<-a[match(ph$sample_id,a$sample_id),];dat<-ph;dat$score<-a$score
  for(v in names(outcome_specs)){
   cov<-outcome_specs[[v]];mr<-residualize(dat$score,dat,c('group'));yr<-residualize(dat[[v]],dat,cov);ct<-safe_cor_test(mr,yr,method='spearman')
   out[[i]]<-data.frame(omics=a$omics[1],module=a$module[1],module_id=mid,outcome=v,module_adjustment='group',outcome_adjustment=paste(cov,collapse=' + '),n=unname(ct['n']),spearman_rho=unname(ct['estimate']),pvalue=unname(ct['pvalue']));i<-i+1
  }
 }
 z<-do.call(rbind,out);z$padj_within_outcome<-NA_real_;for(v in unique(z$outcome)){ix<-z$outcome==v;z$padj_within_outcome[ix]<-p.adjust(z$pvalue[ix],'BH')};z$padj_global<-p.adjust(z$pvalue,'BH');z
}
baseline<-calculate(old);j<-match(key(baseline),key(ref));stopifnot(!anyNA(j));checks<-baseline[,c('module_id','outcome','module_adjustment','outcome_adjustment','n')]
for(f in c('spearman_rho','pvalue','padj_within_outcome','padj_global'))checks[[paste0('abs_error_',f)]]<-abs(baseline[[f]]-ref[[f]][j])
stopifnot(all(baseline$n==ref$n[j]),all(baseline$module_adjustment==ref$module_adjustment[j]),all(baseline$outcome_adjustment==ref$outcome_adjustment[j]))
checks$PASS<-apply(checks[,grep('abs_error_',names(checks)),drop=FALSE],1,function(x)all(is.finite(x)&x<1e-10));ts(checks,'历史160复现逐项审核.tsv');ts(baseline,'历史160复现结果.tsv');if(!all(checks$PASS))stop('OLD_REPLICATION_FAILED; no new sensitivity analysis permitted')
z<-calculate(new);stopifnot(nrow(z)==160,all(z$n==24));ts(z,'01_160项残差Spearman完整结果.tsv')
l<-read.delim(file.path(S,'D_models/core160_final.tsv'),check.names=FALSE,fileEncoding='UTF-8-BOM');j<-match(key(z),key(l));o<-match(key(z),key(ref));stopifnot(!anyNA(j))
d<-z;d$linear_beta<-l$standardized_beta[j];d$linear_P<-l$pvalue[j];d$linear_within_FDR<-l$padj_within_outcome[j];d$linear_global_FDR<-l$padj_global[j];d$direction_relation<-ifelse(sign(d$spearman_rho)==sign(d$linear_beta),'same_direction','opposite_direction');d$old_sensitivity_rho<-ref$spearman_rho[o];d$old_sensitivity_P<-ref$pvalue[o];d$old_sensitivity_within_FDR<-ref$padj_within_outcome[o];d$old_sensitivity_global_FDR<-ref$padj_global[o]
for(t in c(.05,.1)){tag<-if(t==.05)'005'else '010';for(f in c('within','global')){a<-if(f=='within')z$padj_within_outcome else z$padj_global;b<-if(f=='within')d$linear_within_FDR else d$linear_global_FDR;oldq<-if(f=='within')d$old_sensitivity_within_FDR else d$old_sensitivity_global_FDR;d[[paste0(f,'_FDR',tag,'_diff_linear')]]<-(a<t)!=(b<t);d[[paste0(f,'_FDR',tag,'_changed_from_old_sensitivity')]]<-(a<t)!=(oldq<t)}}
d$nominal_diff_linear<-(z$pvalue<.05)!=(d$linear_P<.05);d$nominal_changed_from_old_sensitivity<-(z$pvalue<.05)!=(d$old_sensitivity_P<.05);d$direction_changed_from_old_sensitivity<-sign(z$spearman_rho)!=sign(d$old_sensitivity_rho);ts(d,'02_线性模型与敏感性结果对照.tsv')
ts(d[d$direction_relation=='opposite_direction'|d$nominal_diff_linear|d$within_FDR010_diff_linear|d$global_FDR010_diff_linear,],'方向或统计状态差异条目.tsv')
ts(data.frame(family=c('old_sensitivity160','new_sensitivity160','Stage6_linear160'),nominal005=c(sum(ref$pvalue<.05),sum(z$pvalue<.05),sum(l$pvalue<.05)),within010=c(sum(ref$padj_within_outcome<.1),sum(z$padj_within_outcome<.1),sum(l$padj_within_outcome<.1)),global010=c(sum(ref$padj_global<.1),sum(z$padj_global<.1),sum(l$padj_global<.1)),within005=c(sum(ref$padj_within_outcome<.05),sum(z$padj_within_outcome<.05),sum(l$padj_within_outcome<.05)),global005=c(sum(ref$padj_global<.05),sum(z$padj_global<.05),sum(l$padj_global<.05))),'统计状态汇总.tsv')
writeLines(c(capture.output(sessionInfo()),'No RNG calls; seed not applicable. Only module-phenotype160 residual Spearman computed; no cross-omics132 rerun.'),file.path(O,'sessionInfo.txt'));cat('OLD_REPLICATION_PASS; NEW160_COMPLETE\n')
