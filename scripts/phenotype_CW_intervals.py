from pathlib import Path
import os
import json,math,copy
import pandas as pd,numpy as np
R=Path(os.environ["A1_INPUT_ROOT"]); O=Path(os.environ["A1_OUTPUT"]); S3=R/"09_output/Article1_R2_Stage3_20261029"; C=R/"Article1_core_scripts_20260812"
def rd(p):return pd.read_csv(p,sep='\t' if str(p).endswith('.tsv')else',')
def out(n,d):d.to_csv(O/'results'/n,sep='\t',index=False,float_format='%.17g')
old=rd(C/'output_abdominal_fat_article1_core/Table0_abdominal_fat_clean_data.csv');new=rd(O/'results/Abdominal_fat_clean_CW_missing.tsv')
groups=['YY','YX','XX'];trs=['AF_wt','AF_live_pct','AF_carcass_pct','AF_evis_pct'];rng=np.random.default_rng(20261029)
base=rd(S3/'01_AF/AF_positive_group_means_CI.tsv');pm=base.copy();validation=[]
for g in groups:
    a=old[(old.group==g)&(old.AF_wt>0)][trs].to_numpy();st=copy.deepcopy(rng.bit_generator.state);boots=a[rng.integers(0,len(a),size=(2000,len(a)))].mean(axis=1)
    for j,t in enumerate(trs):
        x=base[(base.group==g)&(base.trait==t)].iloc[0];ci=np.quantile(boots[:,j],[.025,.975]);validation.append(abs(ci[0]-x.group_mean_CI_low));validation.append(abs(ci[1]-x.group_mean_CI_high))
    if g=='YY':
        v=new[(new.group==g)&(new.AF_wt>0)].AF_carcass_pct.dropna().to_numpy();r2=np.random.default_rng();r2.bit_generator.state=st;bs=v[r2.integers(0,len(v),size=(2000,len(v)))].mean(axis=1);ci=np.quantile(bs,[.025,.975]);idx=pm.index[(pm.group==g)&(pm.trait=='AF_carcass_pct')][0]
        for k,vv in dict(N=len(v),mean=v.mean(),SD=v.std(ddof=1),SE=v.std(ddof=1)/np.sqrt(len(v)),group_mean_CI_low=ci[0],group_mean_CI_high=ci[1]).items():pm.loc[idx,k]=vv
assert max(validation)<1e-10;out('Positive_group_means_CI_updated.tsv',pm)
cw_record=rd(Path(os.environ['A1_CW_RECORD_FILE']))
assert list(cw_record)==['wing_id','group','expected_cw_g'] and len(cw_record)==1
# Preserve the original PCG stream position per trait; only authorized changed cells are regenerated.
clean=rd(C/'output_story_phenotype_article1_core/Table0_clean_phenotype_with_derived_traits.csv');cn=clean.copy();cn.loc[cn.wing_id==cw_record.iloc[0]['wing_id'],'Carcass_wt']=np.nan
for d in [clean,cn]:d['AF_carcass_pct']=d.AF_wt/d.Carcass_wt*100;d['AF_evis_pct']=d.AF_wt/d.Evis_wt*100
sub=rd(S3/'02_subset/selected_vs_unselected_effects_CI.tsv');sn=sub.copy();rng=np.random.default_rng(20261029)
trb=['BW_wk15','ADG_0_15','Live_wt','Carcass_wt','Evis_wt']+trs
for g in groups:
    for t in trb:
        st=copy.deepcopy(rng.bit_generator.state);a=clean[clean.group==g];x=a[a.is_any_omics][t].dropna().to_numpy();y=a[~a.is_any_omics][t].dropna().to_numpy()
        xb=x[rng.integers(0,len(x),size=(2000,len(x)))].mean(axis=1);yb=y[rng.integers(0,len(y),size=(2000,len(y)))].mean(axis=1)
        if g=='YY'and t in ['Carcass_wt','AF_carcass_pct']:
            a=cn[cn.group==g];xx=a[a.is_any_omics][t].dropna().to_numpy();yy=a[~a.is_any_omics][t].dropna().to_numpy();r2=np.random.default_rng();r2.bit_generator.state=st
            xbb=xx[r2.integers(0,len(xx),size=(2000,len(xx)))].mean(axis=1);ybb=yy[r2.integers(0,len(yy),size=(2000,len(yy)))].mean(axis=1)
            idx=sn.index[(sn.group==g)&(sn.trait==t)][0]
            vals=dict(N_selected=len(xx),N_unselected=len(yy),N_missing_unselected=int(len(a[~a.is_any_omics])-len(yy)),mean_selected=xx.mean(),mean_unselected=yy.mean(),effect_mean_difference_selected_minus_unselected=xx.mean()-yy.mean())
            for prefix,v in [('selected_mean',xbb),('unselected_mean',ybb),('effect',xbb-ybb)]:lo,hi=np.quantile(v,[.025,.975]);vals[prefix+'_CI_low']=lo;vals[prefix+'_CI_high']=hi
            for k,v in vals.items():sn.loc[idx,k]=v
sn=sn.drop(columns=['physical_sample_pairing','random_selection_validated'],errors='ignore')
sn['R3_CW_value_override']=(sn.group=='YY')&sn.trait.isin(['Carcass_wt','AF_carcass_pct'])
out('S5B_profiled_comparison_updated.tsv',sn)

# Exact Stage3 definitions and R3 positive-only AFW/CW pairwise branch.
def wilcox(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);x=x[np.isfinite(x)];y=y[np.isfinite(y)]
    n,m=len(x),len(y);assert n and m
    vals=np.concatenate([x,y]);order=np.argsort(vals,kind='stable');ranks=np.empty(len(vals),float);ties=[]
    i=0
    while i<len(order):
        j=i+1
        while j<len(order) and vals[order[j]]==vals[order[i]]:j+=1
        ranks[order[i:j]]=(i+1+j)/2;ties.append(j-i);i=j
    u=float(ranks[:n].sum()-n*(n+1)/2);N=n+m
    variance=n*m/12*((N+1)-sum(t**3-t for t in ties)/(N*(N-1)))
    delta=u-n*m/2
    if variance<=0:return dict(W_U=u,z='',P=1.0,variance=variance,ties=max(ties)>1)
    corrected=delta-(.5 if delta>0 else -.5 if delta<0 else 0)
    z=corrected/math.sqrt(variance);p=math.erfc(abs(z)/math.sqrt(2))
    return dict(W_U=u,z=z,P=min(1.0,p),variance=variance,ties=max(ties)>1)

def bh(ps):
    ps=np.asarray(ps,float);order=np.argsort(ps,kind='stable');q=np.empty(len(ps));s=ps[order]*len(ps)/np.arange(1,len(ps)+1)
    q[order]=np.minimum(1,np.minimum.accumulate(s[::-1])[::-1]);return q

an=rd(O/'results/Abdominal_fat_clean_CW_missing.tsv');pairs=[]
for a,b in [('YY','YX'),('YY','XX'),('YX','XX')]:
    x=an[(an.group==a)&(an.AF_wt>0)].AF_carcass_pct.dropna().values;y=an[(an.group==b)&(an.AF_wt>0)].AF_carcass_pct.dropna().values;z=wilcox(x,y);pairs.append(dict(trait='AF_carcass_pct',group1=a,group2=b,N1=len(x),N2=len(y),W_U_group1=z['W_U'],z_continuity_corrected=z['z'],P_two_sided=z['P'],ties_present=z['ties']))
for r,q in zip(pairs,bh([r['P_two_sided']for r in pairs])):r['P_BH_within_trait']=q
pw=pd.DataFrame(pairs);pw.to_csv(O/'results/Positive_AF_CW_pairwise_updated.tsv',sep='\t',index=False,float_format='%.17g')
