# Only the authorized CW correction and direct dependencies; source files are never sourced/executed.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr)})
O<-Sys.getenv('A1_OUTPUT');C<-file.path(Sys.getenv('A1_INPUT_ROOT'),'Article1_core_scripts_20260812'); code_dir<-Sys.getenv('A1_CODE'); record<-read.delim(Sys.getenv('A1_CW_RECORD_FILE'),colClasses=c('character','character','numeric'))
stopifnot(identical(names(record),c('wing_id','group','expected_cw_g')),nrow(record)==1L,
          all(complete.cases(record)),nzchar(record$wing_id),record$group=='YY',
          is.finite(record$expected_cw_g),record$expected_cw_g>0)
target<-record$wing_id;
 dir.create(file.path(O,'results'),recursive=TRUE); dir.create(file.path(O,'audit')); dir.create(file.path(O,'logs'));PH<-file.path(C,'output_story_phenotype_article1_core');AF<-file.path(C,'output_abdominal_fat_article1_core')
wr<-function(d,n)write.table(d,file.path(O,'results',paste0(n,'.tsv')),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
old<-read.csv(file.path(PH,'Table0_clean_phenotype_with_derived_traits.csv'));stopifnot(sum(old$wing_id==target,na.rm=TRUE)==1, old$Carcass_wt[match(target,old$wing_id)]==record$expected_cw_g, old$group[match(target,old$wing_id)]==record$group, !old$is_any_omics[match(target,old$wing_id)]);new<-old;new$Carcass_wt[new$wing_id==target]<-NA_real_;new$Carcass_live_pct<-new$Carcass_wt/new$Live_wt*100
old$group<-factor(old$group,levels=c('YY','YX','XX'));new$group<-factor(new$group,levels=c('YY','YX','XX'))
for(e in parse(file.path(code_dir,'phenotype_population.R'))){if(is.call(e)&&is.symbol(e[[1]])&&as.character(e[[1]])%in%c('<-','=')&&((is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')))||identical(e[[2]],as.name('trait_labels'))))eval(e)}
GROUP_LEVELS<-c('YY','YX','XX');LOWER_IS_BETTER_TRAITS<-c('AF_wt','AF_live_pct');BOOT_N<-2000
# Numerically identical bootstrap sampling and percentile arithmetic, avoiding construction of 2,000 small tibbles.
bs<-function(d,tr,kind='pheno',unit='%'){
 v<-lapply(GROUP_LEVELS,function(g){x<-d[[tr]][d$group==g];x[!is.na(x)]});names(v)<-GROUP_LEVELS
 metric<-function(yy,yx,xx){mp<-(yy+xx)/2;hi<-max(yy,xx);lo<-min(yy,xx)
  if(kind=='pheno')return(c(MPH_percent=(yx-mp)/mp*100,BPH_high_parent_percent=(yx-hi)/hi*100,Directional_BPH_percent=if(tr%in%LOWER_IS_BETTER_TRAITS)(lo-yx)/lo*100 else(yx-hi)/hi*100,F1_position=(yx-lo)/(hi-lo)))
  if(unit=='%')c(Favorable_MPH=(mp-yx)/mp*100,Favorable_BPH=(lo-yx)/lo*100,F1_position=(yx-lo)/(hi-lo))else c(Favorable_MPH=mp-yx,Favorable_BPH=lo-yx,F1_position=(yx-lo)/(hi-lo))
 }
 obs<-do.call(metric,unname(lapply(v,mean)));boot<-replicate(2000,do.call(metric,unname(lapply(v,function(x)mean(sample(x,length(x),replace=TRUE),na.rm=TRUE)))))
 out<-c(N_YY=length(v[[1]]),N_YX=length(v[[2]]),N_XX=length(v[[3]]),Mean_YY=mean(v[[1]]),Mean_YX=mean(v[[2]]),Mean_XX=mean(v[[3]]),obs)
 for(j in seq_along(obs)){nm<-names(obs)[j];prefix<-if(nm=='BPH_high_parent_percent')'BPH_high'else sub('_percent$','',nm);out[paste0(prefix,'_CI_low')]<-quantile(boot[j,],.025,na.rm=TRUE);out[paste0(prefix,'_CI_high')]<-quantile(boot[j,],.975,na.rm=TRUE)}
 as.data.frame(as.list(out))
}
het<-read.csv(file.path(PH,'Table5_heterosis_and_F1_relative_position.csv'));states<-list();rep<-list();updates<-list();set.seed(20260610)
for(tr in het$Trait){states[[tr]]<-.Random.seed;z<-bs(old,tr);h<-het[het$Trait==tr,];cs<-intersect(names(z),names(h));err<-max(abs(unlist(z[cs])-unlist(h[cs])),na.rm=TRUE);rep[[tr]]<-data.frame(branch='phenotype',trait=tr,max_absolute_error=err)
 if(tr%in%c('Carcass_wt','Carcass_live_pct')){afterold<-.Random.seed;.Random.seed<-states[[tr]];zz<-bs(new,tr);updates[[tr]]<-cbind(Trait=tr,zz);.Random.seed<-afterold}
}
wr(bind_rows(rep),'Historical_bootstrap_reproduction_pheno');stopifnot(max(sapply(rep,function(x)x$max_absolute_error))<1e-8)
hetnew<-het
for(tr in names(updates)){i<-match(tr,hetnew$Trait);for(n in intersect(names(updates[[tr]]),names(hetnew)))hetnew[i,n]<-updates[[tr]][[n]]}
for(tr in names(updates)){i<-match(tr,hetnew$Trait);hetnew$F1_pattern[i]<-classify_position(hetnew$Mean_YY[i],hetnew$Mean_YX[i],hetnew$Mean_XX[i])}
wr(hetnew,'Phenotype_heterosis_complete_updated');wr(bind_rows(updates),'CW_and_CW_LW_heterosis_updated')
gt<-read.csv(file.path(PH,'Table3_global_group_tests.csv'));gtn<-gt
for(tr in c('Carcass_wt','Carcass_live_pct')){a<-global_trait_test(old,tr);h<-gt[gt$Trait==tr,];stopifnot(max(abs(a$ANOVA_P-h$ANOVA_P),abs(a$Kruskal_P-h$Kruskal_P))<1e-9);z<-global_trait_test(new,tr);i<-match(tr,gtn$Trait);for(n in intersect(names(z),names(gtn)))gtn[i,n]<-z[[n]]}
gtn$ANOVA_FDR<-p.adjust(gtn$ANOVA_P,'BH');gtn$Kruskal_FDR<-p.adjust(gtn$Kruskal_P,'BH');wr(gtn,'Full_phenotype_global_tests_BH_updated')
su<-read.csv(file.path(PH,'Table2_trait_summary_by_group.csv'));sun<-su
for(tr in c('Carcass_wt','Carcass_live_pct'))for(g in GROUP_LEVELS){x<-new[[tr]][new$group==g];x<-x[!is.na(x)];i<-which(sun$Trait==tr & sun$group==g);if(!length(i))i<-which(sun$Trait==tr & sun$Group==g);for(n in c('N','Mean','SD','SE','Median','Min','Max'))sun[i,n]<-switch(n,N=length(x),Mean=mean(x),SD=sd(x),SE=sd(x)/sqrt(length(x)),Median=median(x),Min=min(x),Max=max(x))}
wr(sun,'Phenotype_group_summary_updated');wr(bind_rows(lapply(c('Carcass_wt','Carcass_live_pct'),function(t)pairwise_trait_test(new,t))),'CW_pairwise_BH_updated')

# Abdominal-fat helper definitions only; no original top-level fits or plots are executed.
for(e in parse(file.path(code_dir,'phenotype_abdominal_fat.R'))){if(is.call(e)&&is.symbol(e[[1]])&&as.character(e[[1]])%in%c('<-','=')&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')))eval(e)}
aold<-read.csv(file.path(AF,'Table0_abdominal_fat_clean_data.csv'));stopifnot(sum(aold$wing_id==target,na.rm=TRUE)==1,aold$Carcass_wt[match(target,aold$wing_id)]==record$expected_cw_g,aold$group[match(target,aold$wing_id)]==record$group);anew<-aold;anew$Carcass_wt[anew$wing_id==target]<-NA_real_;anew$AF_carcass_pct<-anew$AF_wt/anew$Carcass_wt*100
aold$group<-factor(aold$group,levels=GROUP_LEVELS);anew$group<-factor(anew$group,levels=GROUP_LEVELS)
fitold<-lm(AF_wt~Carcass_wt,data=aold,na.action=na.exclude);fitnew<-lm(AF_wt~Carcass_wt,data=anew,na.action=na.exclude)
stopifnot(max(abs(residuals(fitold)-aold$AF_resid_carcasswt),na.rm=TRUE)<1e-8)
anew$AF_resid_carcasswt<-residuals(fitnew);anew$AF_expected_by_Carcass_wt<-fitted(fitnew)
wr(anew,'Abdominal_fat_clean_CW_missing');wr(data.frame(term=names(coef(fitnew)),old_beta=coef(fitold),new_beta=coef(fitnew)),'CW_residual_model_coefficients')
sen<-read.csv(file.path(AF,'Table6_zero_AF_sensitivity_heterosis.csv'));snew<-sen;set.seed(20260610);checks<-list();scenarios<-unique(sen$Scenario)
repl<-c(0,NA,.5,1,min(aold$AF_wt[aold$AF_wt>0],na.rm=TRUE)/2)
for(k in seq_along(scenarios)){sc<-scenarios[k];aa<-make_sensitivity_dataset(aold,repl[k],k==2);bb<-make_sensitivity_dataset(anew,repl[k],k==2)
 for(tr in c('AF_wt_sens','AF_live_pct_sens','AF_carcass_pct_sens','AF_evis_pct_sens')){st<-.Random.seed;z<-bs(aa,tr,'af');h<-sen[sen$Scenario==sc&sen$Trait==tr,];cs<-intersect(names(z),names(h));err<-max(abs(unlist(z[cs])-unlist(h[cs])),na.rm=TRUE);checks[[length(checks)+1]]<-data.frame(branch=sc,trait=tr,max_absolute_error=err)
  if(tr=='AF_carcass_pct_sens'){after<-.Random.seed;.Random.seed<-st;zz<-bs(bb,tr,'af');.Random.seed<-after;i<-which(snew$Scenario==sc&snew$Trait==tr);for(n in intersect(names(zz),names(snew)))snew[i,n]<-zz[[n]];snew$Kruskal_P[i]<-kruskal.test(bb[[tr]]~bb$group)$p.value}
 }}
main<-read.csv(file.path(AF,'Table5_favorable_heterosis_abdominal_fat_traits.csv'));mnew<-main;afposition<-list()
for(i in seq_len(nrow(main))){tr<-main$Trait[i];st<-.Random.seed;z<-bs(aold,tr,'af',main$Unit[i]);cs<-intersect(names(z),names(main));err<-max(abs(unlist(z[cs])-unlist(main[i,cs])),na.rm=TRUE);checks[[length(checks)+1]]<-data.frame(branch='AF5_main',trait=tr,max_absolute_error=err)
 if(tr%in%c('AF_carcass_pct','AF_resid_carcasswt')){after<-.Random.seed;.Random.seed<-st;zz<-bs(anew,tr,'af',main$Unit[i]);.Random.seed<-after;afposition[[tr]]<-cbind(Trait=tr,zz);for(n in intersect(names(zz),names(mnew)))mnew[i,n]<-zz[[n]]}}
wr(bind_rows(checks),'Historical_AF_bootstrap_reproduction');stopifnot(max(sapply(checks,function(x)x$max_absolute_error))<1e-8)
for(tr in c('AF_carcass_pct','AF_resid_carcasswt')){i<-match(tr,mnew$Trait);mnew$Favorable_parent[i]<-if(mnew$Mean_YY[i]<=mnew$Mean_XX[i])'YY'else'XX'}
wr(mnew,'AF5_heterosis_updated');wr(snew,'AF6_zero_sensitivity_updated');wr(bind_rows(afposition),'AF_CW_and_residual_position_updated')
# Only affected group summaries/tests are updated; the full original BH family is retained.
ags<-read.csv(file.path(AF,'Table2_abdominal_fat_trait_summary_by_group.csv'));agsnew<-ags
for(tr in c('Carcass_wt','AF_carcass_pct','AF_resid_carcasswt'))for(g in GROUP_LEVELS){x<-anew[[tr]][anew$group==g];x<-x[!is.na(x)];i<-which(agsnew$Trait==tr & agsnew$group==g);for(n in intersect(c('N','Mean','SD','SE','Median','Min','Max'),names(agsnew)))agsnew[i,n]<-switch(n,N=length(x),Mean=mean(x),SD=sd(x),SE=sd(x)/sqrt(length(x)),Median=median(x),Min=min(x),Max=max(x))}
wr(agsnew,'AF_group_summary_updated')
extras<-bind_rows(lapply(c('AF_carcass_pct','AF_evis_pct'),function(tr){z<-global_trait_test(anew,tr);z$Trait_label<-if(tr=='AF_carcass_pct')'AFW/CW'else'AFW/EW';z}))
g36<-bind_rows(gtn,extras);stopifnot(nrow(g36)==36);g36$ANOVA_FDR<-p.adjust(g36$ANOVA_P,'BH');g36$Kruskal_FDR<-p.adjust(g36$Kruskal_P,'BH');wr(g36,'Full_phenotype_global_tests_BH_updated')
for(i in seq_len(nrow(hetnew))){j<-match(hetnew$Trait[i],g36$Trait);for(k in c('ANOVA_P','ANOVA_FDR','Kruskal_P','Kruskal_FDR'))if(k%in%names(hetnew))hetnew[i,k]<-g36[j,k]};wr(hetnew,'Phenotype_heterosis_complete_updated')
zz<-new[new$group=='YY' & is.finite(new$Carcass_wt),];sel<-tolower(as.character(zz$is_any_omics))=='true';wr(data.frame(group='YY',trait='Carcass_wt',Wilcoxon_P=wilcox.test(zz$Carcass_wt[sel],zz$Carcass_wt[!sel],exact=FALSE)$p.value),'S5A_CW_pairwise_P_updated')
wr(bind_rows(lapply(c('AF_carcass_pct','AF_resid_carcasswt'),function(tr){z<-anew[anew$AF_wt>0 & is.finite(anew[[tr]]),];data.frame(Trait=tr,group=GROUP_LEVELS,N=sapply(GROUP_LEVELS,function(g)sum(z$group==g)),Mean=sapply(GROUP_LEVELS,function(g)mean(z[z$group==g,tr])),SD=sapply(GROUP_LEVELS,function(g)sd(z[z$group==g,tr])),Median=sapply(GROUP_LEVELS,function(g)median(z[z$group==g,tr])),Min=sapply(GROUP_LEVELS,function(g)min(z[z$group==g,tr])),Max=sapply(GROUP_LEVELS,function(g)max(z[z$group==g,tr])))})),'AF7_positive_CW_summary_updated')
affected<-c('Carcass_wt','AF_carcass_pct','AF_resid_carcasswt')
aftests<-bind_rows(lapply(affected,function(tr){z<-anew[is.finite(anew[[tr]]),];data.frame(Trait=tr,ANOVA_P=summary(aov(z[[tr]]~z$group))[[1]][['Pr(>F)']][1],Kruskal_P=kruskal.test(z[[tr]]~z$group)$p.value)}));wr(aftests,'Affected_AF_global_tests')
wr(bind_rows(lapply(affected,function(tr)pairwise_trait_test(anew,tr))),'Affected_AF_pairwise_BH')
for(i in seq_len(nrow(agsnew))){tr<-agsnew$Trait[i];if(tr%in%affected && 'Kruskal_P'%in%names(agsnew))agsnew$Kruskal_P[i]<-aftests$Kruskal_P[match(tr,aftests$Trait)]};wr(agsnew,'AF_group_summary_updated')
cat(capture.output(summary(fitnew)),sep='\n')
cat('CW_CORE_RECALCULATION_COMPLETE; no omics models refitted\n');capture.output(sessionInfo(),file=file.path(O,'logs/R_sessionInfo_R3.txt'));saveRDS(list(pheno_states=states,seed=20260610,bootstrap=2000),file.path(O,'audit/bootstrap_protocol_states_R3.rds'))
