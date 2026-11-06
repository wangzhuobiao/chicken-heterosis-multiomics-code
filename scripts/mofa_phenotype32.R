options(digits=17)
O<-Sys.getenv('A1_OUTPUT');S<-Sys.getenv('A1_MOFA_REFERENCE')
tryCatch({
 # Use exactly the three Stage7 function/specification definitions, not its executable workflow.
 nodes<-parse(file.path(Sys.getenv('A1_CODE'),'frozen_definitions.R'))
 for(e in nodes) if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]]) %in% c('spec','fit_one'))eval(e)
 safe_numeric <- function(x) {x<-trimws(as.character(x));x[x %in% c('','NA','NaN','NULL','null','-','--')]<-NA_character_;suppressWarnings(as.numeric(gsub(',','',x,fixed=TRUE)))}
 stopifnot(exists('spec'),exists('safe_numeric'),exists('fit_one'))
 pheno<-read.csv(file.path(Sys.getenv('A1_INPUT_ROOT'),'Article1_multiomics/Article1_multiomics_scripts/output_multiomics_phenotype_anchored_article1_core/Table_core_phenotypes_for_integration.csv'),check.names=FALSE)
 pheno$group<-factor(pheno$group,levels=c('YY','YX','XX'))
 z<-read.delim(file.path(O,'results/stage7_aligned_factor_scores.tsv'),check.names=FALSE)
 stopifnot(nrow(z)==24,!anyDuplicated(z$sample_id),setequal(z$sample_id,pheno$sample_id))
 z<-z[match(pheno$sample_id,z$sample_id),];stopifnot(all(as.character(pheno$group)==z$group))
 df<-cbind(pheno,z[,grep('^Factor',names(z)),drop=FALSE]);factors<-grep('^Factor',names(df),value=TRUE);rows<-list()
 for(p in names(spec))for(f in factors)rows[[length(rows)+1L]]<-fit_one(df,p,f,spec[[p]])
 res<-do.call(rbind,rows);stopifnot(nrow(res)==32,all(res$n==24))
 res$padj_BH_within_phenotype<-ave(res$pvalue,res$phenotype,FUN=function(x)p.adjust(x,'BH'));res$padj_BH_global<-p.adjust(res$pvalue,'BH')
 res$within_family_n<-4L;res$global_family_n<-32L
 res$nominal_P_lt_0.05<-res$pvalue<.05;res$within_FDR_lt_0.10<-res$padj_BH_within_phenotype<.1;res$within_FDR_lt_0.05<-res$padj_BH_within_phenotype<.05;res$global_FDR_lt_0.10<-res$padj_BH_global<.1;res$global_FDR_lt_0.05<-res$padj_BH_global<.05
 write.table(res,file.path(O,'results/factor_phenotype_complete.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
 old<-read.delim(file.path(S,'results/new_factor_phenotype_complete.tsv'),check.names=FALSE)
 cmp<-merge(old,res,by=c('phenotype','factor'),suffixes=c('_stage7','_stage7_1'));stopifnot(nrow(cmp)==32)
 # Tolerances read from the immutable pre-training file, never changed after result inspection.
 tol<-jsonlite::fromJSON(file.path(dirname(Sys.getenv('A1_CODE')),'config','tolerances_frozen.json'))$phenotype_statistics
 measures<-c('standardized_beta','SE','t_value','pvalue','model_R2','adjusted_R2','padj_BH_within_phenotype','padj_BH_global')
 for(k in measures){a<-cmp[[paste0(k,'_stage7')]];b<-cmp[[paste0(k,'_stage7_1')]];cmp[[paste0(k,'_abs_delta')]]<-abs(a-b);cmp[[paste0(k,'_PASS')]]<-abs(a-b)<=tol$atol+tol$rtol*abs(a)}
 flags<-c('nominal_P_lt_0.05','within_FDR_lt_0.10','within_FDR_lt_0.05','global_FDR_lt_0.10','global_FDR_lt_0.05')
 for(k in flags)cmp[[paste0(k,'_same')]]<-cmp[[paste0(k,'_stage7')]]==cmp[[paste0(k,'_stage7_1')]]
 cmp$n_same<-cmp$n_stage7==cmp$n_stage7_1;cmp$model_same<-cmp$model_stage7==cmp$model_stage7_1
 cmp$all_PASS<-apply(cmp[,c(paste0(measures,'_PASS'),paste0(flags,'_same'),'n_same','model_same')],1,all)
 write.table(cmp,file.path(O,'results/factor_phenotype_comparison.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 stopifnot(all(cmp$all_PASS))
 dir.create(file.path(O,'logs'),showWarnings=FALSE)
 capture.output(sessionInfo(),file=file.path(O,'logs/R_factor_phenotype_sessionInfo.txt'))
 writeLines('FACTOR_PHENOTYPE_NUMERICAL_PASS',file.path(O,'results/factor_phenotype_PASS.txt'))
 print(c(N=nrow(res),nominal=sum(res$nominal_P_lt_0.05),within10=sum(res$within_FDR_lt_0.10),global10=sum(res$global_FDR_lt_0.10),max_P_delta=max(cmp$pvalue_abs_delta)))
},error=function(e){writeLines(paste('PHENOTYPE_VALIDATION_FAIL',conditionMessage(e)),file.path(O,'FAILURE_factor_phenotype.txt'));stop(e)})
