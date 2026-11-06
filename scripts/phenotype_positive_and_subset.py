from pathlib import Path
import os
import csv,json,hashlib,sys,platform,math,time,itertools,shutil,datetime
from collections import Counter
import numpy as np
import openpyxl

sys.stdout.reconfigure(encoding='utf-8')
O=Path(os.environ["A1_OUTPUT"]); ROOT=Path(os.environ["A1_INPUT_ROOT"])
C=ROOT/'Article1_core_scripts_20260812'; AF=C/'output_abdominal_fat_article1_core'; PH=C/'output_story_phenotype_article1_core'; RNA=C/'output_duodenum_transcriptome_core'
for name in ['01_AF','02_subset','03_F1','04_candidates','04_candidates/legacy120','04_candidates/core160','05_reproducibility']:(O/name).mkdir(exist_ok=True,parents=True)
SEED=20261029; B=2000; GROUPS=['YY','YX','XX']; PAIRS=[('YY','YX'),('YY','XX'),('YX','XX')]
AF_TRAITS=['AF_wt','AF_live_pct','AF_carcass_pct','AF_evis_pct']
TRAITS=['BW_wk15','ADG_0_15','Live_wt','Carcass_wt','Evis_wt']+AF_TRAITS
LABEL={'AF_wt':'AFW','AF_live_pct':'AFW/LW','AF_carcass_pct':'AFW/CW','AF_evis_pct':'AFW/EW'}
UNITS={t:('g/day' if t=='ADG_0_15' else '%' if t.endswith('_pct') else 'g') for t in TRAITS}
INPUTS={}; LOG=[]
def say(s):
    LOG.append(s); print(s,flush=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,rows,fields=None):
    rows=list(rows)
    assert rows or fields, name
    with (O/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def md(name,s):(O/name).write_text(s.strip()+'\n',encoding='utf-8-sig')
def jout(name,s):(O/name).write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf-8')
def reg(p,role='source input'):
    p=Path(p)
    if str(p) not in INPUTS:INPUTS[str(p)]={'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'role':role}
    return p
def table(p):
    p=reg(p)
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t' if p.suffix in ['.tsv','.txt'] else ','))
def number(s):
    if s is None or str(s).strip() in ['','NA','nan','None']:return float('nan')
    return float(s)
def present(s):return s is not None and str(s).strip() not in ['','NA','nan','None']
def same(a,b):
    a,b=number(a),number(b)
    return (not math.isfinite(a) and not math.isfinite(b)) or abs(a-b)<1e-9

settings={
 'scope':'new AF positive-only means/CIs and within-trait pairwise Wilcoxon; record-based selected-vs-not-selected effect estimates; correction candidates only',
 'seed':SEED,'replicates':B,'rng':'numpy.random.Generator(PCG64); separate generators for Task A and Task B, each initialized once with 20261029',
 'CI':'pointwise 95% percentile; np.quantile([0.025,0.975],method=linear) (R quantile type 7)',
 'A':{'population':'frozen cleaned cohort, AF_wt>0; no additional exclusions',
      'order':GROUPS,'traits':AF_TRAITS,'resampling':'within genetic group, resample individuals n times with replacement; same sampled rows across all four AF traits; 2000 replicates',
      'pairwise':'two-sided independent Wilcoxon rank-sum/Mann-Whitney; average ties; normal approximation; tie-adjusted variance; continuity correction 0.5; equivalent to R exact=FALSE, correct=TRUE',
      'BH_family':'three pairwise comparisons YY-YX, YY-XX, YX-XX separately per AF trait; no pooled correction across traits',
      'historical_MPH_BPH':'copy all AF5/AF6/AF7 source bytes; transfer existing CIs with explicit branch labels; no recomputation',
      'Fisher':'copy historical 3x2 overall P=0.04114; not recomputed; not pairwise'},
 'B':{'traits':TRAITS,'membership':'exact frozen 24-bird key by (group,wing_id); eight per group; no ID edits',
      'estimand':'mean(selected) - mean(not selected), in original trait units; descriptive association, not causal',
      'order':'YY,YX,XX; traits in listed order; selected resampled first, then unselected; finite values only within each trait/status cell',
      'bootstrap':'independent resampling within each genetic-group/selection-status/trait cell, preserving observed finite N; difference CIs from paired independent bootstrap means',
      'zero_handling':'all original AF zero records retained; missing only if already missing in frozen tables',
      'P_values':'none requested or computed; no multiple-testing family of new P-values; CIs are pointwise, exploratory, not simultaneous',
      'historical_comparison':'copy existing subset-vs-full-cohort representativeness table separately; never use overlapping subset/full as independent samples'},
 'C':{'identity_gate':'NOT_APPROVED/UNVERIFIED; no new F1 analysis, including fallback VST bootstrap',
      'planned_tolerances':[0,.02,.05,.10],'planned_full_refit_replicates':100,'automatic_extension_to_500':False},
 'limitations':['individual-level bootstrap assumes exchangeable independent birds within cells; cage identifiers/cluster effects not available or modeled',
                'no automatic outlier removal, threshold revision, ID correction or new missing-value imputation',
                'sample physical identity and random selection mechanism remain unverified']}
jout('05_reproducibility/analysis_settings_PREDEFINED.json',settings)

clean_p=PH/'Table0_clean_phenotype_with_derived_traits.csv'; af_p=AF/'Table0_abdominal_fat_clean_data.csv'; key_p=PH/'Table0_omics_sample_key_for_integration.csv'
clean=table(clean_p);af=table(af_p);keys=table(key_p)
def key(r):return r['group'],r['wing_id']
assert len({key(r) for r in clean})==len(clean)==234
assert len({key(r) for r in af})==len(af)==234
assert Counter(r['group'] for r in keys)=={'YY':8,'YX':8,'XX':8}
selected={key(r) for r in keys};assert len(selected)==24
cm={key(r):r for r in clean};am={key(r):r for r in af}
assert cm.keys()==am.keys() and selected<=cm.keys()
for k,r in cm.items():
    for t in ['BW_wk15','Live_wt','Carcass_wt','Evis_wt','AF_wt','AF_live_pct']:assert same(r[t],am[k][t]),f'Source conflict {k} {t}'
    assert (r['is_any_omics']=='TRUE')==(k in selected)
    for t in ['AF_carcass_pct','AF_evis_pct']:r[t]=am[k][t]
    if k in selected:
        kr=next(x for x in keys if key(x)==k)
        for f in ['cecum_meta_id','duodenum_rna_id']:assert kr[f]==r[f]
        for t in ['BW_wk15','ADG_0_15','Live_wt','Carcass_wt','Evis_wt','AF_wt','AF_live_pct']:assert same(kr[t],r[t])

# Check AF records against the original supplied workbook without repairing it.
raw_p=reg(Path(os.environ["A1_PHENOTYPE_XLSX"]),'frozen original phenotype workbook, read only')
wb=openpyxl.load_workbook(raw_p,read_only=True,data_only=True);rawrows=list(wb['Sheet1'].iter_rows(values_only=True));wb.close()
rawmap={};heads=list(rawrows[0])
for i,v in enumerate(rawrows[1:],2):
    if v[0] in GROUPS and v[1] is not None:rawmap.setdefault((str(v[0]),str(v[1])),[]).append((i,dict(zip(heads,v))))
afcheck=[];membership=[]
for k,r in cm.items():
    assert len(rawmap.get(k,[]))==1,f'Ambiguous raw key {k}'
    excelrow,rv=rawmap[k][0]
    if math.isfinite(number(r['AF_wt'])):
        assert same(rv['腹脂重'],r['AF_wt']),f'AF original conflict {k}'
        assert number(r['AF_wt'])>=0
        for field,rawfield in [('Live_wt','活重'),('Carcass_wt','屠体重'),('Evis_wt','全净堂重')]:assert same(r[field],rv[rawfield])
        afcheck.append({'group':k[0],'wing_id':k[1],'source_original':str(raw_p),'source_cell':f'Sheet1!AT{excelrow}',
                        'original_AFW':rv['腹脂重'],'frozen_AFW':r['AF_wt'],'is_zero':number(r['AF_wt'])==0,'selected':k in selected,'check':'PASS'})
    membership.append({'group':k[0],'wing_id':k[1],'cecum_meta_id':r['cecum_meta_id'],'duodenum_D_id':r['duodenum_rna_id'],
        'selected_frozen_key':k in selected,'source_key':str(key_p),'original_row':excelrow,'original_path':str(raw_p),
        'AFW':r['AF_wt'],'is_slaughter_as_stored':r['is_slaughter'],'physical_omics_identity':'UNVERIFIED'})
write('05_reproducibility/AF_raw_vs_clean_validation.tsv',afcheck)
membership.sort(key=lambda r:(GROUPS.index(r['group']),r['wing_id']))
write('02_subset/selected_vs_unselected_membership.tsv',membership)
assert len(afcheck)==59
say('PREFLIGHT: all 59 nonmissing AFW and body-size records match original workbook; 24 selected IDs unchanged.')

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

# Verify the implementation against actual prior R output before new pairwise tests.
historical_p=table(AF/'Table2B_pairwise_wilcoxon_AF_traits.csv');validation=[]
for t in AF_TRAITS:
    tests=[wilcox([number(r[t]) for r in af if r['group']==g1],[number(r[t]) for r in af if r['group']==g2]) for g1,g2 in PAIRS]
    qs=bh([v['P'] for v in tests])
    for (g1,g2),v,q in zip(PAIRS,tests,qs):
        h=next(r for r in historical_p if r['Trait']==t and {r['Group1'],r['Group2']}=={g1,g2})
        delta=abs(float(q)-float(h['Wilcoxon_BH_P']));assert delta<1e-12,(t,g1,g2,q,h)
        validation.append({'trait':t,'group1':g1,'group2':g2,'historical_R_BH':h['Wilcoxon_BH_P'],'reproduced_BH':q,'absolute_difference':delta,
            'pass_tolerance_1e12':True,'source_path':str(AF/'Table2B_pairwise_wilcoxon_AF_traits.csv'),'source_fields':'Trait;Group1;Group2;Wilcoxon_BH_P'})
write('05_reproducibility/Wilcoxon_12_historical_R_validation.tsv',validation)
say('VALIDATION: all 12 retained-zero historical R pairwise BH P-values match within 1e-12.')

start=time.perf_counter();rngA=np.random.default_rng(SEED);meansA=[];countsA=[];bootA={}
for g in GROUPS:
    positive=[r for r in af if r['group']==g and number(r['AF_wt'])>0]
    mat=np.array([[number(r[t]) for t in AF_TRAITS] for r in positive]);assert np.isfinite(mat).all()
    inds=rngA.integers(0,len(mat),size=(B,len(mat)))
    boots=mat[inds].mean(axis=1)
    bootA[g]=boots
    for j,t in enumerate(AF_TRAITS):
        lo,hi=np.quantile(boots[:,j],[.025,.975],method='linear')
        meansA.append({'trait':t,'label':LABEL[t],'unit':UNITS[t],'scenario':'Exclude zeros (AF_wt>0)','group':g,'N':len(mat),
            'mean':mat[:,j].mean(),'SD':mat[:,j].std(ddof=1),'SE':mat[:,j].std(ddof=1)/math.sqrt(len(mat)),
            'group_mean_CI_low':lo,'group_mean_CI_high':hi,'CI_target':'GROUP MEAN, not MPH/BPH','CI_method':'individual within-group percentile bootstrap; quantile linear/type7','B':B,'seed':SEED,
            'source_path':str(af_p),'source_fields':f'group;AF_wt>0;{t}','source_SHA256':INPUTS[str(af_p)]['sha256']})
        countsA.append({'task':'A','group':g,'trait':t,'status':'positive AF_wt','N':len(mat),'zero_count':0})
write('01_AF/AF_positive_group_means_CI.tsv',meansA)
np.savez_compressed(O/'01_AF/TaskA_bootstrap_group_mean_replicates.npz',**bootA)
pairA=[]
for t in AF_TRAITS:
    tests=[]
    for g1,g2 in PAIRS:
        x=[number(r[t]) for r in af if r['group']==g1 and number(r['AF_wt'])>0]
        y=[number(r[t]) for r in af if r['group']==g2 and number(r['AF_wt'])>0]
        test=wilcox(x,y);tests.append({'trait':t,'label':LABEL[t],'unit':UNITS[t],'scenario':'Exclude zeros','group1':g1,'group2':g2,'N1':len(x),'N2':len(y),
            'W_U_group1':test['W_U'],'z_continuity_corrected':test['z'],'P_two_sided':test['P'],'P_BH_within_trait':'','ties_present':test['ties'],
            'method':'Wilcoxon rank-sum normal approximation, average ties, tie variance correction, continuity=TRUE','family':f'{t}: three YY-YX/YY-XX/YX-XX comparisons','family_size':3,
            'source_path':str(af_p),'source_fields':f'group;AF_wt>0;{t}','source_SHA256':INPUTS[str(af_p)]['sha256']})
    for r,q in zip(tests,bh([r['P_two_sided'] for r in tests])):r['P_BH_within_trait']=float(q)
    pairA+=tests
write('01_AF/AF_positive_pairwise_Wilcoxon.tsv',pairA)

# Historical branches copied as bytes and also merged in a clearly labelled long-form supplement.
historical=[]
for fname,branch in [('Table5_favorable_heterosis_abdominal_fat_traits.csv','AF5_main_historical'),('Table6_zero_AF_sensitivity_heterosis.csv','AF6_sensitivity_historical'),('Table7_positive_only_AF_trait_summary.csv','AF7_positive_descriptive_historical')]:
    p=reg(AF/fname);dest=O/'01_AF'/('HISTORICAL_UNCHANGED_'+fname);shutil.copyfile(p,dest);assert sha(p)==sha(dest)
    for i,r in enumerate(table(p),1):
        t=r['Trait'].removesuffix('_sens')
        if t in AF_TRAITS:historical.append({'branch':branch,'source_path':str(p),'source_record':i,'source_SHA256':sha(p),'original_record':r})
jout('01_AF/historical_AF_branches.json',historical)
fpath=reg(AF/'Table3B_measurable_abdominal_fat_fisher_test.txt');assert 'p-value = 0.04114' in fpath.read_text(encoding='utf-8-sig')
write('01_AF/Fisher_overall_existing.tsv',[{'test':'historical Fisher exact; overall 3x2; not pairwise','YY_measurable':20,'YY_nonmeasurable':0,'YX_measurable':15,'YX_nonmeasurable':5,'XX_measurable':17,'XX_nonmeasurable':2,
    'P_as_printed':'0.04114','newly_computed':False,'source_path':str(fpath),'source_locator':'text L8-12;L19-23','source_SHA256':sha(fpath)}])

# Full Reviewer #2.2 supplement: disjoint statistic columns prevent branch/CI confusion.
cols=['record_type','trait','scenario','branch','group','group1','group2','N','N_YY','N_YX','N_XX','mean','mean_YY','mean_YX','mean_XX',
    'group_mean_CI_low','group_mean_CI_high','MPH','MPH_CI_low','MPH_CI_high','BPH','BPH_CI_low','BPH_CI_high','overall_Kruskal_P','overall_Fisher_P','pairwise_raw_P','pairwise_BH_P',
    'unit','unit_as_stored','CI_target','method','source_path','source_locator','source_SHA256','status']
full=[]
def addfull(**kw):
    r={k:'' for k in cols};r.update(kw);full.append(r)
for r in meansA:addfull(record_type='new_positive_group_mean',trait=r['trait'],scenario=r['scenario'],branch='Stage3_A',group=r['group'],N=r['N'],mean=r['mean'],
    group_mean_CI_low=r['group_mean_CI_low'],group_mean_CI_high=r['group_mean_CI_high'],unit=r['unit'],CI_target='group mean',method=r['CI_method']+'; B=2000 seed=20261029',
    source_path=r['source_path'],source_locator=r['source_fields'],source_SHA256=r['source_SHA256'],status='NEW_APPROVED_STATISTIC')
for r in pairA:addfull(record_type='new_positive_pairwise',trait=r['trait'],scenario=r['scenario'],branch='Stage3_A',group1=r['group1'],group2=r['group2'],
    pairwise_raw_P=r['P_two_sided'],pairwise_BH_P=r['P_BH_within_trait'],unit=r['unit'],method=r['method']+'; '+r['family'],source_path=r['source_path'],source_locator=r['source_fields'],source_SHA256=r['source_SHA256'],status='NEW_APPROVED_STATISTIC')
for v in historical:
    r=v['original_record'];t=r['Trait'].removesuffix('_sens')
    if 'Favorable_MPH' in r:
        addfull(record_type='existing_heterosis',trait=t,scenario=r.get('Scenario','main retained-zero AF5'),branch=v['branch'],N_YY=r['N_YY'],N_YX=r['N_YX'],N_XX=r['N_XX'],
            mean_YY=r['Mean_YY'],mean_YX=r['Mean_YX'],mean_XX=r['Mean_XX'],MPH=r['Favorable_MPH'],MPH_CI_low=r['Favorable_MPH_CI_low'],MPH_CI_high=r['Favorable_MPH_CI_high'],
            BPH=r['Favorable_BPH'],BPH_CI_low=r['Favorable_BPH_CI_low'],BPH_CI_high=r['Favorable_BPH_CI_high'],overall_Kruskal_P=r.get('Kruskal_P',''),
            unit=UNITS[t],unit_as_stored=r['Unit'],CI_target='existing MPH/BPH bootstrap intervals; not group mean CIs',method='historical only; not recomputed',source_path=v['source_path'],source_locator=f'CSV logical record {v["source_record"]}; Scenario/Trait',source_SHA256=v['source_SHA256'],status='HISTORICAL_BRANCH_PRESERVED')
    else:addfull(record_type='existing_positive_mean',trait=t,scenario='Exclude zeros',branch=v['branch'],group=r['group'],N=r['N'],mean=r['Mean'],unit=UNITS[t],method='historical AF7 descriptive mean; no historical mean CI',source_path=v['source_path'],source_locator=f'CSV logical record {v["source_record"]}',source_SHA256=v['source_SHA256'],status='HISTORICAL_BRANCH_PRESERVED')
addfull(record_type='existing_Fisher_overall',trait='AFW>0 status',scenario='all original zero/nonzero records',branch='historical_FISH',N_YY=20,N_YX=20,N_XX=19,
    overall_Fisher_P='0.04114',method='overall 3x2 Fisher; not pairwise; not recomputed',source_path=str(fpath),source_locator='L8-12;L19-23',source_SHA256=sha(fpath),status='HISTORICAL_OUTPUT_REUSED')
write('01_AF/Reviewer2_2_complete_supplement.tsv',full)
say('TASK A COMPLETE: 12 positive-only group means/CIs; 12 pairwise tests; historical AF5/6/7/FISH branches preserved.')

rngB=np.random.default_rng(SEED);subset=[];bootB={}
for g in GROUPS:
    for t in TRAITS:
        sx=np.array([number(r[t]) for r in clean if r['group']==g and key(r) in selected]);ux=np.array([number(r[t]) for r in clean if r['group']==g and key(r) not in selected])
        x=sx[np.isfinite(sx)];y=ux[np.isfinite(ux)];assert len(x)>=2 and len(y)>=2
        xb=x[rngB.integers(0,len(x),size=(B,len(x)))].mean(axis=1);yb=y[rngB.integers(0,len(y),size=(B,len(y)))].mean(axis=1)
        bootB[g+'__'+t+'__selected_mean']=xb;bootB[g+'__'+t+'__unselected_mean']=yb
        xlo,xhi=np.quantile(xb,[.025,.975],method='linear');ylo,yhi=np.quantile(yb,[.025,.975],method='linear');dlo,dhi=np.quantile(xb-yb,[.025,.975],method='linear')
        subset.append({'group':g,'trait':t,'label':LABEL.get(t,t),'unit':UNITS[t],'N_selected':len(x),'N_unselected':len(y),'N_missing_selected':len(sx)-len(x),'N_missing_unselected':len(ux)-len(y),
            'mean_selected':x.mean(),'selected_mean_CI_low':xlo,'selected_mean_CI_high':xhi,'mean_unselected':y.mean(),'unselected_mean_CI_low':ylo,'unselected_mean_CI_high':yhi,
            'effect_mean_difference_selected_minus_unselected':x.mean()-y.mean(),'effect_CI_low':dlo,'effect_CI_high':dhi,'effect_type':'unstandardized difference of means',
            'effect_unit':'percentage points' if t.endswith('_pct') else UNITS[t],
            'N_zero_selected':int(np.sum(x==0)),'N_zero_unselected':int(np.sum(y==0)),
            'CI_method':'2000 independent within-status/trait bootstrap means; percentile/type7; pointwise not simultaneous','seed':SEED,'B':B,
            'source_phenotype':str(clean_p),'source_AF_extra_ratios':str(af_p) if t in ['AF_carcass_pct','AF_evis_pct'] else '',
            'source_phenotype_SHA256':INPUTS[str(clean_p)]['sha256'],'source_AF_SHA256':INPUTS[str(af_p)]['sha256'],
            'membership_source':str(key_p),'membership_SHA256':INPUTS[str(key_p)]['sha256'],'physical_sample_pairing':'UNVERIFIED','random_selection_validated':False,
            'interpretation':'recorded selection association; comparison is selected vs NOT selected; not subset vs overlapping full cohort'})
write('02_subset/selected_vs_unselected_effects_CI.tsv',subset)
np.savez_compressed(O/'02_subset/TaskB_bootstrap_group_means_replicates.npz',**bootB)
assert all(r['N_selected']==8 and r['N_unselected']==12 for r in subset if r['group']=='YX' and r['trait'] in AF_TRAITS)
oldrep=reg(PH/'Table6_omics_sample_representativeness.csv');shutil.copyfile(oldrep,O/'02_subset/HISTORICAL_subset_vs_full_cohort.csv')
say('TASK B COMPLETE: 27 within-group selected-vs-unselected trait comparisons; no new P-values; AF zeros retained.')

