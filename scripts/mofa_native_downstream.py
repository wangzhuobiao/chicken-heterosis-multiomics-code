"""Independent native HDF5 validation and frozen-tolerance comparisons. NO training/writes to HDF5."""
import sys
import os
sys.dont_write_bytecode=True
from pathlib import Path
import csv,json,re,hashlib,traceback
import numpy as np,h5py
from scipy.optimize import linear_sum_assignment
O=Path(os.environ['A1_OUTPUT']); S=Path(os.environ['A1_MOFA_REFERENCE']); A=Path(os.environ['A1_MOFA_ARCHIVE']); C=Path(__file__).resolve().parents[1]/'config'; inputs=Path(os.environ.get('A1_MOFA_INPUTS') or A/'inputs')
for folder in ('results','validation'):(O/folder).mkdir(exist_ok=True)
T=json.loads((C/'tolerances_frozen.json').read_text());cfg=json.loads((C/'training_config.json').read_text());samples=json.loads((inputs/'samples.json').read_text());views=cfg['views']
def rd(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def wr(p,rs):
 with (O/p).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
def dec(x):return [v.decode() if isinstance(v,bytes) else str(v) for v in np.atleast_1d(x)]
def close(a,b,rule):return bool(np.allclose(a,b,**T[rule],equal_nan=False))
archive=[];num=[]
def ck(rows,field,passed,detail):rows.append(dict(check=field,PASS=bool(passed),detail=str(detail)))
def parse_log(p):
 s=p.read_text(encoding='utf-8-sig');out=[]
 for m in re.finditer(r'Iteration (\d+): time=([\d.]+),(?: ELBO=([-\d.]+), deltaELBO=([-\d.]+) \(([-\d.]+)%\),)? Factors=(\d+)',s):out.append(dict(iteration=int(m[1]),n_factors=int(m[6]),ELBO_log_rounded=m[3] or '',deltaELBO_log_rounded=m[4] or '',relative_delta_percent=m[5] or ''))
 return s,out
try:
 model=A/'model/Article1_MOFA2_Stage7_1_native.hdf5';old=S/'model/Article1_MOFA2_Stage7_repaired.hdf5'
 log,trace=parse_log(A/'logs/training_once.log');oldlog,oldtrace=parse_log(S/'logs/training_once.log')
 ck(archive,'training_converged','Converged!' in log,'native stdout convergence marker')
 ck(num,'iteration_count_and_factor_trajectory',[(r['iteration'],r['n_factors']) for r in trace]==[(r['iteration'],r['n_factors']) for r in oldtrace],f'old/new {len(oldtrace)}/{len(trace)} logged iterations')
 wr('results/convergence_from_log.tsv',trace)
 elcomp=[]
 for a,b in zip(oldtrace,trace):
  for field,rule in [('ELBO_log_rounded','ELBO_vs_printed_log'),('deltaELBO_log_rounded','ELBO_delta_vs_printed_log'),('relative_delta_percent','printed_relative_delta_percent')]:
   if a[field]!='':
    ok=b[field]!='' and close(float(a[field]),float(b[field]),rule)
    elcomp.append(dict(iteration=a['iteration'],field=field,Stage7_printed=a[field],Stage7_1_printed=b[field],abs_delta=abs(float(a[field])-float(b[field])) if b[field] else '',PASS=ok))
 wr('results/ELBO_logged_comparison.tsv',elcomp);ck(num,'all_printed_ELBO_points',all(r['PASS'] for r in elcomp),'tolerances fixed to original printed precision')
 with h5py.File(model,'r') as h,h5py.File(old,'r') as a:
  roots=['data','expectations','features','samples','views','groups','intercepts','model_options','training_opts','training_stats','variance_explained']
  for x in roots:ck(archive,'root_'+x,x in h,'native group/dataset exists')
  assert all(r['PASS'] for r in archive),'Native root validation failed'
  inventory=[]
  def visit(n,x):
   if isinstance(x,h5py.Dataset):inventory.append(dict(dataset=n,shape=str(x.shape),dtype=str(x.dtype),attributes=json.dumps({k:dec(v) for k,v in x.attrs.items()}),nan_count=int(np.isnan(x[...]).sum()) if np.issubdtype(x.dtype,np.floating) else 0))
  h.visititems(visit);wr('validation/HDF5_inventory.tsv',inventory)
  ck(archive,'view_order',dec(h['views/views'][...])==views,'ordered view labels')
  ck(archive,'sample_order',dec(h['samples/group1'][...])==samples==dec(a['samples/group1'][...]),'24 frozen samples in exact order')
  ck(archive,'group_names',dec(h['groups/groups'][...])==['group1'],'MOFA single group, genetic groups used downstream')
  Z=h['expectations/Z/group1'][...].T;Z0=a['expectations/Z/group1'][...].T
  ck(archive,'Z_shape_finite',Z.shape==(24,4) and np.isfinite(Z).all(),str(Z.shape));ck(num,'retained_factors',Z.shape[1]==Z0.shape[1],f'{Z0.shape[1]} / {Z.shape[1]}')
  assert all(r['PASS'] for r in archive) and all(r['PASS'] for r in num),'Initial numeric/archive validation failed'
  C=np.corrcoef(Z0.T,Z.T)[:4,4:];ii,jj=linear_sum_assignment(-abs(C));sg=np.sign(C[ii,jj]);aligned=Z[:,jj]*sg
  matches=[dict(Stage7_factor=f'Factor{i+1}',native_factor=f'Factor{j+1}',Pearson_r=float(C[i,j]),sign_to_Stage7=int(s)) for i,j,s in zip(ii,jj,sg)]
  wr('results/factor_match.tsv',matches)
  wr('results/all_factor_correlations.tsv',[dict(Stage7_factor=f'Factor{i+1}',native_factor=f'Factor{j+1}',Pearson_r=float(C[i,j])) for i in range(4) for j in range(4)])
  wr('results/stage7_aligned_factor_scores.tsv',[dict(sample_id=s,group=s[:2],**{f'Factor{k+1}':float(v) for k,v in enumerate(row)}) for s,row in zip(samples,aligned)])
  wr('validation/native_factor_scores.tsv',[dict(sample_id=s,group=s[:2],**{f'Factor{k+1}':float(v) for k,v in enumerate(row)}) for s,row in zip(samples,Z)])
  wz=[]
  for n,s in enumerate(samples):
   for i in range(4):wz.append(dict(kind='Z',view='ALL',element_id=s,Stage7_factor=f'Factor{i+1}',native_factor=f'Factor{jj[i]+1}',Stage7_value=float(Z0[n,i]),aligned_new_value=float(aligned[n,i]),abs_delta=abs(float(Z0[n,i]-aligned[n,i])),PASS=close(Z0[n,i],aligned[n,i],'W_Z')))
  R2=h['variance_explained/r2_per_factor/group1'][...];RT=h['variance_explained/r2_total/group1'][...]
  ck(archive,'R2_shape_finite',R2.shape==(3,4) and RT.shape==(3,) and np.isfinite(R2).all() and np.isfinite(RT).all(),f'{R2.shape}, {RT.shape}')
  stats=h['training_stats'];fac=stats['number_factors'][...];times=stats['time'][...];elb=stats['elbo'][...]
  ck(archive,'native_training_stats_lengths',len(fac)==len(times)==len(elb)==len(trace),'native [:i] excludes convergence iteration; inspect source, no manual extension')
  ck(archive,'native_training_stats_valid',np.isfinite(fac).all() and (fac>0).all() and np.isfinite(times).all() and (times>=0).all() and np.isfinite(elb).any() and not np.isinf(elb).any(),'NaN ELBO expected at unassessed iterations')
  expected_indices=[0]+list(range(cfg['training']['startELBO'],len(elb),cfg['training']['freqELBO']))
  ck(archive,'ELBO_evaluation_positions',list(np.flatnonzero(np.isfinite(elb)))==expected_indices,'initial ELBO plus startELBO/freqELBO positions')
  statrows=[dict(iteration=i,n_factors=int(fac[i]),time_seconds=float(times[i]),ELBO_native_full_precision=float(elb[i]) if np.isfinite(elb[i]) else '',ELBO_evaluated=bool(np.isfinite(elb[i]))) for i in range(len(elb))]
  wr('results/native_training_stats.tsv',statrows)
  for r in trace:
   i=r['iteration']
   if i<len(fac):ck(num,f'native_factor_count_iteration_{i}',int(fac[i])==r['n_factors'],'native vs log')
   if i<len(elb) and r['ELBO_log_rounded']!='':ck(archive,f'native_ELBO_vs_log_{i}',close(elb[i],float(r['ELBO_log_rounded']),'ELBO_vs_printed_log'),'full native value vs logged rounding')
  names=dec(h['training_opts'].attrs['names']);opts=dict(zip(names,h['training_opts'][...].tolist()));expected_opts=dict(maxiter=cfg['training']['iter'],freqELBO=cfg['training']['freqELBO'],start_elbo=cfg['training']['startELBO'],gpu_mode=cfg['training']['gpu_mode'],stochastic=False,seed=cfg['training']['seed'])
  ck(archive,'training_opts_values_and_names',opts==expected_opts,json.dumps(opts))
  for k in ['likelihoods','spikeslab_factors','spikeslab_weights','ard_factors','ard_weights']:
   expected=cfg['likelihoods'] if k=='likelihoods' else [str(cfg['model'][k])]
   ck(archive,'model_option_'+k,dec(h['model_options/'+k][...])==expected,str(dec(h['model_options/'+k][...])))
  varcmp=[];weightlong=[];oldvar=rd(S/'results/variance_explained_comparison.tsv')
  for vi,v in enumerate(views):
   f=json.loads((inputs/f'{v}_features.json').read_text());ck(archive,'feature_order_'+v,dec(h['features/'+v][...])==f==dec(a['features/'+v][...]),'exact300 ordered IDs')
   Y=h['data/'+v+'/group1'][...];ex=np.load(A/'validation'/f'{v}_expected_backend_data.npy');inter=h['intercepts/'+v+'/group1'][...];exi=np.nanmean(np.load(inputs/f'{v}.npy').astype(np.float64),axis=0) # preserve frozen NPY layout as actual set_data_matrix does
   ck(archive,'saved_input_'+v,Y.dtype==ex.dtype and np.array_equal(Y,ex),'float64 promotion then group centering; bitwise equality')
   ck(archive,'intercept_'+v,np.array_equal(inter,exi),'exact pre-centering column means')
   W=h['expectations/W/'+v][...];W0=a['expectations/W/'+v][...];wa=W[jj]*sg[:,None]
   ck(archive,'W_shape_finite_'+v,W.shape==(4,300) and np.isfinite(W).all(),str(W.shape))
   ck(num,'W_allclose_'+v,close(W0,wa,'W_Z'),f'max abs delta={np.max(abs(W0-wa))}')
   for i in range(4):
    for k,feat in enumerate(f):wz.append(dict(kind='W',view=v,element_id=feat,Stage7_factor=f'Factor{i+1}',native_factor=f'Factor{jj[i]+1}',Stage7_value=float(W0[i,k]),aligned_new_value=float(wa[i,k]),abs_delta=abs(float(W0[i,k]-wa[i,k])),PASS=close(W0[i,k],wa[i,k],'W_Z')))
   native_total=100*(1-np.square(Y-Z@W).sum()/np.square(Y).sum())
   ck(archive,'native_total_R2_formula_'+v,close(RT[vi],native_total,'variance_percent_internal'),f'native={RT[vi]}, recomputed={native_total}')
   for i in range(4):
    j=jj[i];r=next(r for r in oldvar if r['view']==v and r['new_factor']==f'Factor{i+1}');ov=float(r['new_R2_percent']);nv=float(R2[vi,j]);rv=100*(1-np.square(Y-np.outer(Z[:,j],W[j])).sum()/np.square(Y).sum())
    ck(archive,f'native_factor_R2_formula_{v}_{j+1}',close(nv,rv,'variance_percent_internal'),f'{nv} vs {rv}')
    varcmp.append(dict(view=v,Stage7_factor=f'Factor{i+1}',native_factor=f'Factor{j+1}',Stage7_helper_R2_percent=ov,native_R2_percent=nv,abs_delta=abs(ov-nv),PASS=close(ov,nv,'variance_percent_vs_Stage7_helper')))
   ov=float(next(r for r in oldvar if r['view']==v and r['new_factor']=='ALL')['new_R2_percent']);nv=float(RT[vi]);varcmp.append(dict(view=v,Stage7_factor='ALL',native_factor='ALL',Stage7_helper_R2_percent=ov,native_R2_percent=nv,abs_delta=abs(ov-nv),PASS=close(ov,nv,'variance_percent_vs_Stage7_helper')))
   wr(f'validation/native_{v}_weights.tsv',[dict(feature_id=feat,**{f'Factor{j+1}':float(W[j,k]) for j in range(4)}) for k,feat in enumerate(f)])
   wr(f'validation/native_{v}_data.tsv',[dict(sample_id=s,**{feat:float(Y[n,k]) for k,feat in enumerate(f)}) for n,s in enumerate(samples)])
   for j in range(4):
    ranks=np.argsort(-abs(W[j]),kind='stable');rank={int(k):q+1 for q,k in enumerate(ranks)}
    for k,feat in enumerate(f):weightlong.append(dict(view=v,native_factor=f'Factor{j+1}',feature_id=feat,weight=float(W[j,k]),rank_abs=rank[k]))
  ck(num,'Z_allclose',close(Z0,aligned,'W_Z'),f'max abs delta={np.max(abs(Z0-aligned))}')
  wr('results/W_Z_full_comparison.tsv',wz);wr('results/variance_comparison.tsv',varcmp);wr('results/all_weights.tsv',weightlong);wr('results/top20_weights_per_view_factor.tsv',[r for r in weightlong if r['rank_abs']<=20])
  ck(num,'all_W_Z_entries',all(r['PASS'] for r in wz),'3600 W + 96 Z entries')
  ck(num,'all_variance_entries',all(r['PASS'] for r in varcmp),'native vs Stage7 helper, frozen percent tolerance')
  objective=json.loads((A/'results/final_objective_readonly.json').read_text())['terms']['total'];old_final=float(oldtrace[-1]['ELBO_log_rounded'])
  ck(num,'final_ELBO_vs_Stage7_printed',close(objective,old_final,'ELBO_vs_printed_log'),f'full final new {objective} vs old printed {old_final}; native stats last finite index {expected_indices[-1]}')
 wr('validation/archive_checks.tsv',archive);wr('results/numerical_checks_prephenotype.tsv',num)
 state=dict(archive_native_checks_PASS=all(r['PASS'] for r in archive),numerical_prephenotype_PASS=all(r['PASS'] for r in num),standard_loader_pending=True,phenotype_comparison_pending=True)
 (O/'validation/INDEPENDENT_HDF5_STATUS.json').write_text(json.dumps(state,indent=2))
 assert state['archive_native_checks_PASS'] and state['numerical_prephenotype_PASS'],'Archive/numerical validation failed; stop replacement'
 print(json.dumps(state,indent=2));print('W/Z max delta',max(r['abs_delta'] for r in wz));print('variance max delta',max(r['abs_delta'] for r in varcmp));print('final ELBO',objective)
except BaseException:
 if archive:wr('validation/archive_checks.tsv',archive)
 if num:wr('results/numerical_checks_prephenotype.tsv',num)
 (O/'FAILURE_independent_validation.txt').write_text(traceback.format_exc(),encoding='utf-8');raise
