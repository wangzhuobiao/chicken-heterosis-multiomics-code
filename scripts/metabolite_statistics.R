args<-commandArgs(trailingOnly=TRUE);O<-normalizePath(args[1],winslash='/',mustWork=TRUE);R<-Sys.getenv('A1_INPUT_ROOT')
options(stringsAsFactors=FALSE,digits=17,scipen=30)
for(d in c('A_DE','B_F1','C_modules','D_models','E_MOFA','logs'))dir.create(file.path(O,d),showWarnings=FALSE)
rd<-function(p)read.csv(p,check.names=FALSE,fileEncoding='UTF-8-BOM')
ts<-function(x,p)write.table(x,file.path(O,p),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
repair<-Sys.getenv('A1_REPAIR');core<-file.path(R,'Article1_core_scripts_20260812');M<-file.path(core,'output_plasma_metabolome_article1_core');I<-file.path(R,'Article1_multiomics/Article1_multiomics_scripts');C<-file.path(I,'output_multiomics_phenotype_anchored_article1_core')
gate<-jsonlite::fromJSON(file.path(repair,'validation/independent_output_QA.json'));stopifnot(gate$status=='INDEPENDENT_SERIALIZED_MATRIX_VERIFICATION_PASS')
map<-read.delim(file.path(repair,'matrix/feature_key_mapping_5272.tsv'),check.names=FALSE,fileEncoding='UTF-8-BOM');meta<-read.delim(file.path(repair,'validation/sample_order_24.tsv'),fileEncoding='UTF-8-BOM');s<-meta$sample_id
map$is_ms2<-as.logical(map$is_ms2);map$expected_numeric_repair<-as.logical(map$expected_numeric_repair)
x<-as.matrix(rd(file.path(repair,'matrix/processed_log2_5272x24.csv')));rownames(x)<-map$legacy_feature_id;stopifnot(all(colnames(x)==s),all(dim(x)==c(5272,24)),!anyDuplicated(map$new_feature_key),all(table(meta$group)==8),all(is.finite(x)))
ann<-read.delim(file.path(repair,'matrix/feature_annotations_5272.tsv'),check.names=FALSE,fileEncoding='UTF-8-BOM');stopifnot(all(ann$legacy_feature_id==rownames(x)));ms<-map$is_ms2;stopifnot(sum(ms)==383)
oldx<-rd(file.path(M,'tables/Table2_log2_processed_abundance_matrix.csv'));ox<-as.matrix(oldx[,s]);rownames(ox)<-oldx$feature_id;stopifnot(all(rownames(ox)==rownames(x)))
contr<-list(YY_vs_XX=c('YY','XX'),YX_vs_YY=c('YX','YY'),YX_vs_XX=c('YX','XX'));res<-list();logs<-list();start<-proc.time()[3]
for(ct in names(contr)){
 g<-contr[[ct]];a<-x[,meta$group==g[1],drop=FALSE];b<-x[,meta$group==g[2],drop=FALSE]
 p<-W<-rep(NA_real_,nrow(x));status<-rep('OK',nrow(x));warn<-rep('',nrow(x))
 for(i in seq_len(nrow(x))){z<-tryCatch(withCallingHandlers(wilcox.test(a[i,],b[i,],alternative='two.sided',exact=FALSE,correct=TRUE),warning=function(w){warn[i]<<-conditionMessage(w);invokeRestart('muffleWarning')}),error=function(e){status[i]<<-conditionMessage(e);NULL});if(!is.null(z)){p[i]<-z$p.value;W[i]<-unname(z$statistic)};if(!is.finite(p[i])&&status[i]=='OK')status[i]<-'NOT_TESTABLE'}
 delta<-rowMeans(a)-rowMeans(b);qsub<-rep(NA_real_,nrow(x));qsub[ms]<-p.adjust(p[ms],method='BH',n=383)
 z<-data.frame(legacy_feature_id=rownames(x),new_feature_key=map$new_feature_key,contrast=ct,numerator=g[1],denominator=g[2],N1=8,N2=8,mean_log2processed_1=rowMeans(a),mean_log2processed_2=rowMeans(b),delta_mean_log2processed=delta,direction=ifelse(delta>0,'Up',ifelse(delta<0,'Down','No_change')),W=W,rank_biserial=2*W/64-1,pvalue=p,q_all5272_primary=p.adjust(p,'BH',n=5272),q_MS2_383_secondary=qsub,test_status=status,warnings=warn)
 z<-cbind(z,ann[,setdiff(names(ann),c('legacy_feature_id','new_feature_key')),drop=FALSE]);res[[ct]]<-z
 cat(ct,'done\n');flush.console()
}
all<-do.call(rbind,res);ts(all,'A_DE/all5272_three_contrasts.tsv');ts(all[all$is_ms2,],'A_DE/MS2_383_three_contrasts.tsv')
old<-rd(file.path(M,'tables/Table5_all_pairwise_differential_metabolites.csv'));j<-match(paste(all$legacy_feature_id,all$contrast),paste(old$feature_id,old$contrast));stopifnot(!anyNA(j));comp<-all[,1:16];comp$old_P<-old$pvalue[j];comp$old_full5272_BH<-old$padj[j];comp$P_change<-all$pvalue-comp$old_P;comp$primary_q_change<-all$q_all5272_primary-comp$old_full5272_BH;comp$old_mean_difference<-old$log2FC[j];comp$abundance_repaired<-map$expected_numeric_repair[match(comp$legacy_feature_id,map$legacy_feature_id)];ts(comp,'A_DE/new_old_feature_statistics.tsv')
counts<-do.call(rbind,lapply(names(res),function(ct)do.call(rbind,lapply(c('all5272','MS2_383'),function(sub){z<-res[[ct]];if(sub=='MS2_383')z<-z[z$is_ms2,];data.frame(contrast=ct,subset=sub,N=nrow(z),nominal005=sum(z$pvalue<.05,na.rm=TRUE),primary_FDR005=sum(z$q_all5272_primary<.05,na.rm=TRUE),primary_FDR010=sum(z$q_all5272_primary<.1,na.rm=TRUE),secondary_FDR005=if(sub=='MS2_383')sum(z$q_MS2_383_secondary<.05,na.rm=TRUE)else NA,secondary_FDR010=if(sub=='MS2_383')sum(z$q_MS2_383_secondary<.1,na.rm=TRUE)else NA,not_testable=sum(!is.finite(z$pvalue)))}))))
ts(counts,'A_DE/counts_FDR_summary.tsv')
# Extract the exact frozen classifier function without executing the historical script.
code<-parse(file.path(Sys.getenv('A1_CODE'),'frozen_definitions.R'),encoding='UTF-8');for(e in code)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('classify_f1_pattern')))eval(e)
PARENTAL_DIFF_MIN<-.30;POSITION_TOLERANCE<-.15;MEAN_TOLERANCE_MIN<-.15
means<-function(m)sapply(c('YY','YX','XX'),function(g)rowMeans(m[,meta$group==g,drop=FALSE]))
bm<-means(ox);nm<-means(x);oldpat<-rd(file.path(M,'tables/Table8_F1_metabolite_inheritance_patterns.csv'));jj<-match(rownames(x)[ms],oldpat$feature_id);base<-apply(bm[ms,,drop=FALSE],1,function(v)classify_f1_pattern(v[1],v[2],v[3]));stopifnot(all(base==oldpat$pattern[jj]))
nc<-apply(nm[ms,,drop=FALSE],1,function(v)classify_f1_pattern(v[1],v[2],v[3]));fp<-data.frame(legacy_feature_id=rownames(x)[ms],new_feature_key=map$new_feature_key[ms],name=map$annotation_name[ms],old_class=base,new_class=nc,mean_YY=nm[ms,1],mean_YX=nm[ms,2],mean_XX=nm[ms,3],parental_range=abs(nm[ms,1]-nm[ms,3]),tolerance=pmax(abs(nm[ms,1]-nm[ms,3])*.15,.15));fp$old_classifiable<-fp$old_class!='unresolved';fp$new_classifiable<-fp$new_class!='unresolved';ts(fp,'B_F1/all383_new_old_classification.tsv');ts(fp[fp$legacy_feature_id=='M397T467.1',],'B_F1/Mitragynine_classification.tsv')
cc<-table(factor(nc,levels=c('YX_intermediate','YY_dominant','XX_dominant','YX_high_transgressive','YX_low_transgressive','unresolved')));fsummary<-data.frame(pattern=names(cc),N=as.integer(cc));fsummary$percent_all383<-fsummary$N/383*100;fsummary$percent_classified<-ifelse(fsummary$pattern=='unresolved',NA,fsummary$N/sum(nc!='unresolved')*100);ts(fsummary,'B_F1/META9_corrected_new_counts_percentages.tsv')
