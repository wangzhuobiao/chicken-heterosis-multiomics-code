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
# Module scores must reproduce frozen legacy scores before any update.
scores<-rd(file.path(C,'Table_integrated_module_scores_long.csv'));ph<-rd(file.path(C,'Table_core_phenotypes_for_integration.csv'));ph$group<-factor(ph$group,levels=c('YY','YX','XX'));samples<-ph$sample_id
members<-rd(file.path(C,'Table_metabolite_feature_to_module_mapping.csv'));affectedmods<-unique(members$module[members$feature_id%in%map$legacy_feature_id[map$expected_numeric_repair]])
zs<-function(v){d<-sd(v);if(!is.finite(d)||d==0)rep(0,length(v))else as.numeric((v-mean(v))/d)}
scoremod<-function(m,mod){f<-unique(members$feature_id[members$module==mod]);m<-m[f,samples,drop=FALSE];m<-m[apply(m,1,sd)>0,,drop=FALSE];colMeans(t(scale(t(m))))}
updates<-list();nscores<-scores
for(mod in unique(members$module)){
 ids<-which(scores$omics=='Plasma_metabolome'&scores$module==mod);bs<-scoremod(ox,mod);newscore<-scoremod(x,mod);expected<-bs[match(scores$sample_id[ids],names(bs))];stopifnot(max(abs(expected-scores$score[ids]))<1e-10)
 if(!mod%in%affectedmods)stopifnot(max(abs(newscore-bs))<1e-10)
 nscores$score[ids]<-newscore[match(scores$sample_id[ids],names(newscore))]
 updates[[mod]]<-data.frame(module=mod,sample_id=scores$sample_id[ids],old_score=scores$score[ids],new_score=nscores$score[ids],delta=nscores$score[ids]-scores$score[ids],affected=mod%in%affectedmods,baseline_max_error=max(abs(expected-scores$score[ids])))
}
ts(do.call(rbind,updates),'C_modules/new_old_plasma_scores.tsv');ts(nscores,'C_modules/all20_updated_scores.tsv');members$new_feature_key<-map$new_feature_key[match(members$feature_id,map$legacy_feature_id)];stopifnot(!anyNA(members$new_feature_key));ts(members,'C_modules/frozen_members_with_composite_keys.tsv')
# Existing model arithmetic verification: residual projection, no substitute model.
oldmodel<-rd(file.path(C,'Table_module_phenotype_linear_models.csv'));newmodel<-oldmodel;checks<-list()
for(i in seq_len(nrow(oldmodel))){r<-oldmodel[i,];a<-scores[scores$module_id==r$module_id,];a<-a[match(samples,a$sample_id),];b<-nscores[nscores$module_id==r$module_id,];b<-b[match(samples,b$sample_id),];dat<-ph;dat$z_module<-zs(a$score);dat$z_outcome<-zs(dat[[r$outcome]]);cov<-if(r$outcome=='AF_wt')c('group','Live_wt')else 'group';Q<-qr(model.matrix(reformulate(cov),dat));rx<-qr.resid(Q,dat$z_module);ry<-qr.resid(Q,dat$z_outcome);beta<-sum(rx*ry)/sum(rx^2);df<-24-Q$rank-1;se<-sqrt(sum((ry-beta*rx)^2)/df/sum(rx^2));pv<-2*pt(-abs(beta/se),df);err<-max(abs(c(beta-r$standardized_beta,se-r$standard_error,pv-r$pvalue)));stopifnot(err<1e-10)
 affected<-r$omics=='Plasma_metabolome'&&r$module%in%affectedmods
 if(affected){dat$z_module<-zs(b$score);fit<-lm(reformulate(c('z_module',cov),'z_outcome'),dat);co<-coef(summary(fit))['z_module',];newmodel$standardized_beta[i]<-co[1];newmodel$standard_error[i]<-co[2];newmodel$pvalue[i]<-co[4]}
 checks[[i]]<-data.frame(module_id=r$module_id,outcome=r$outcome,old_arithmetic_error=err,refitted=affected)
}
newmodel$padj_global<-p.adjust(newmodel$pvalue,'BH',n=160);for(o in unique(newmodel$outcome)){ix<-newmodel$outcome==o;newmodel$padj_within_outcome[ix]<-p.adjust(newmodel$pvalue[ix],'BH',n=20)}
ts(newmodel,'D_models/core160_final.tsv');ts(do.call(rbind,checks),'D_models/model_source_verification.tsv');mdiff<-merge(oldmodel,newmodel,by=c('module_id','outcome'),suffixes=c('_old','_new'));ts(mdiff,'D_models/core160_old_new.tsv')
oldcor<-rd(file.path(C,'Table_cross_omics_group_adjusted_correlations.csv'));newcor<-oldcor
residscore<-function(tbl,mid){a<-tbl[tbl$module_id==mid,];v<-a$score[match(samples,a$sample_id)];qr.resid(qr(model.matrix(~group,ph)),v)}
corchecks<-list()
for(i in seq_len(nrow(oldcor))){r<-oldcor[i,];a<-residscore(scores,r$module_id1);b<-residscore(scores,r$module_id2);z<-suppressWarnings(cor.test(a,b,method='spearman',exact=FALSE));err<-max(abs(c(unname(z$estimate)-r$spearman_rho,z$p.value-r$pvalue)));stopifnot(err<1e-10)
 affected<-(r$omics1=='Plasma_metabolome'&&r$module1%in%affectedmods)||(r$omics2=='Plasma_metabolome'&&r$module2%in%affectedmods)
 if(affected){z<-suppressWarnings(cor.test(residscore(nscores,r$module_id1),residscore(nscores,r$module_id2),method='spearman',exact=FALSE));newcor$spearman_rho[i]<-unname(z$estimate);newcor$pvalue[i]<-z$p.value}
 corchecks[[i]]<-data.frame(module_id1=r$module_id1,module_id2=r$module_id2,baseline_error=err,recomputed=affected)
}
newcor$padj_global<-p.adjust(newcor$pvalue,'BH',n=132);for(g in unique(newcor$pair_type)){ix<-newcor$pair_type==g;newcor$padj_within_pair_type[ix]<-p.adjust(newcor$pvalue[ix],'BH',n=sum(ix))}
newcor$nominal_edge<-abs(newcor$spearman_rho)>=.45&newcor$pvalue<.05;newcor$fdr_edge<-abs(newcor$spearman_rho)>=.45&newcor$padj_within_pair_type<.1
ts(newcor,'D_models/cross132_final.tsv');ts(do.call(rbind,corchecks),'D_models/correlation_source_verification.tsv');ts(merge(oldcor,newcor,by=c('module_id1','module_id2'),suffixes=c('_old','_new')),'D_models/cross132_old_new.tsv')
ts(data.frame(family=c('core160','cross132'),N=c(160,132),nominal005=c(sum(newmodel$pvalue<.05),sum(newcor$pvalue<.05)),withinFDR010=c(sum(newmodel$padj_within_outcome<.1),sum(newcor$padj_within_pair_type<.1)),globalFDR010=c(sum(newmodel$padj_global<.1),sum(newcor$padj_global<.1)),withinFDR005=c(sum(newmodel$padj_within_outcome<.05),sum(newcor$padj_within_pair_type<.05)),globalFDR005=c(sum(newmodel$padj_global<.05),sum(newcor$padj_global<.05))),'D_models/count_summary.tsv')
writeLines(capture.output(sessionInfo()),file.path(O,'logs/sessionInfo.txt'))
