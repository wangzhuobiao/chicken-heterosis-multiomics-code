O<-normalizePath(commandArgs(trailingOnly=TRUE)[1],winslash='/',mustWork=TRUE);R<-Sys.getenv('A1_INPUT_ROOT');dir.create(file.path(O,'E_MOFA'),showWarnings=FALSE);I<-file.path(R,'Article1_multiomics/Article1_multiomics_scripts');repair<-Sys.getenv('A1_REPAIR');h5<-file.path(I,'output_MOFA2_sensitivity/model/Article1_MOFA2_model.hdf5');s<-rhdf5::h5read(h5,'/samples/group1')
rd<-function(p)read.csv(p,check.names=FALSE,fileEncoding='UTF-8-BOM');ts<-function(x,n)write.table(x,file.path(O,'E_MOFA',n),sep='\t',row.names=FALSE,quote=TRUE)
sel<-function(m){v<-apply(m,1,var);ok<-is.finite(v)&v>0&rowSums(is.finite(m))==ncol(m);m[which(ok)[order(v[ok],decreasing=TRUE)[1:300]],,drop=FALSE]}
z<-function(m){a<-t(scale(t(m)));a[!is.finite(a)]<-0;a}
records<-list()
for(view in c('Duodenum_RNA','Cecal_metagenome','Plasma_metabolome')){
 if(view=='Duodenum_RNA'){d<-rd(file.path(I,'Table3_VST_expression_matrix.csv'));id<-d[[1]]}
 if(view=='Cecal_metagenome'){d<-read.delim(file.path(I,'KEGG_Pathway_TPM_matrix.tsv'),check.names=FALSE);id<-d[[1]]}
 if(view=='Plasma_metabolome'){d<-rd(file.path(I,'Table2_log2_processed_abundance_matrix.csv'));id<-d$feature_id}
 m<-as.matrix(d[,s,drop=FALSE]);storage.mode(m)<-'numeric';rownames(m)<-id
 if(anyDuplicated(id)){m<-rowsum(m,id,reorder=TRUE)/as.numeric(table(id)[rownames(rowsum(m,id,reorder=TRUE))])}
 if(view=='Cecal_metagenome')m<-log2(m+1)
 old<-sel(m);hf<-as.character(rhdf5::h5read(h5,paste0('/features/',view)));hd<-rhdf5::h5read(h5,paste0('/data/',view,'/group1'))
 if(!all(dim(hd)==c(300,24)))hd<-t(hd)
 zz<-z(m[hf,s,drop=FALSE]);err<-max(abs(zz-hd));quantized<-matrix(readBin(writeBin(as.double(zz),raw(),size=4),double(),n=length(zz),size=4),nrow=300);qe<-max(abs(quantized-hd));records[[view]]<-data.frame(view=view,source_selected_set_matches=setequal(rownames(old),hf),source_selected_order_matches=identical(rownames(old),hf),stored_input_max_abs_difference=err,float32_roundtrip_error=qe,reusable_fixed_input=(view!='Plasma_metabolome'&&qe<1e-12))
 if(view=='Plasma_metabolome'){
  map<-read.delim(file.path(repair,'matrix/feature_key_mapping_5272.tsv'),fileEncoding='UTF-8-BOM');nm<-as.matrix(rd(file.path(repair,'matrix/processed_log2_5272x24.csv')));rownames(nm)<-map$legacy_feature_id;nm<-nm[,s];nn<-sel(nm)
  ts(data.frame(new_order=1:300,legacy_feature_id=rownames(nn),new_feature_key=map$new_feature_key[match(rownames(nn),map$legacy_feature_id)],variance=apply(nn,1,var),previous_model_member=rownames(nn)%in%hf),'repaired_top300.tsv')
  ts(data.frame(legacy_feature_id=union(hf,rownames(nn)),old_selected=union(hf,rownames(nn))%in%hf,new_selected=union(hf,rownames(nn))%in%rownames(nn)),'old_new_selection_union.tsv')
  jsonlite::write_json(list(old_selected_set_reproduced=setequal(rownames(old),hf),old_feature_order_reproduced=identical(rownames(old),hf),new_old_overlap=length(intersect(rownames(nn),hf)),entered=setdiff(rownames(nn),hf),left=setdiff(hf,rownames(nn)),retrain_required_for_corrected_results=TRUE,training_run=FALSE),file.path(O,'E_MOFA/input_summary.json'),pretty=TRUE,auto_unbox=TRUE)
 }
}
ts(do.call(rbind,records),'three_views_raw_source_reuse_checks.tsv');writeLines(c(capture.output(sessionInfo()),'HDF5 read only. No training.'),file.path(O,'E_MOFA/input_check_session.txt'))
attrs<-rhdf5::h5readAttributes(h5,'/training_opts');capture.output(str(attrs),file=file.path(O,'E_MOFA/training_opts_attributes.txt'));saveRDS(attrs,file.path(O,'E_MOFA/training_opts_attributes.rds'))
tv<-as.numeric(rhdf5::h5read(h5,'/training_opts'));ts(data.frame(parameter=as.character(attrs$names),stored_value=tv),'stored_training_options.tsv')
