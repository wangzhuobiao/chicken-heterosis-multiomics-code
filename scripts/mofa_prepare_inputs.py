from pathlib import Path
import csv,json,hashlib,os
import numpy as np,h5py
O=Path(os.environ["A1_OUTPUT"]);R=Path(os.environ["A1_INPUT_ROOT"]);S=Path(os.environ["A1_MOFA_SELECTION"]);P=Path(os.environ["A1_REPAIR"]);I=R/"Article1_multiomics/Article1_multiomics_scripts"
def read(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t' if Path(p).suffix=='.tsv' else ','))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def decode(a):return [x.decode() if isinstance(x,bytes) else str(x) for x in a]
for d in ['inputs','model','results','logs']:(O/d).mkdir(exist_ok=True)
model=I/'output_MOFA2_sensitivity/model/Article1_MOFA2_model.hdf5';sources=[model,P/'matrix/processed_log2_5272x24.csv',P/'matrix/feature_key_mapping_5272.tsv',P/'validation/sample_order_24.tsv',S/'E_MOFA/repaired_top300.tsv']
mapping=read(sources[2]);samples=read(sources[3]);raw=np.loadtxt(sources[1],delimiter=',',skiprows=1);sampleids=[r['sample_id'] for r in samples];ids=[r['legacy_feature_id'] for r in mapping];keys={r['legacy_feature_id']:r['new_feature_key'] for r in mapping}
qc=[];feature_alias=[];stored={}
with h5py.File(model,'r') as h:
 hs=decode(h['samples/group1'][...]);assert len(set(hs))==24 and set(hs)==set(sampleids);idx=[sampleids.index(x) for x in hs];raw=raw[:,idx];views=['Duodenum_RNA','Plasma_metabolome','Cecal_metagenome']
 variances=np.var(raw,axis=1,ddof=1);top=np.flatnonzero(np.isfinite(variances)&(variances>0)&np.isfinite(raw).all(1));top=top[np.argsort(-variances[top],kind='stable')[:300]];topids=[ids[i] for i in top];assert topids==[r['legacy_feature_id'] for r in read(sources[4])]
 for view in views:
  oldfeatures=decode(h['features/'+view][...]);olddata=h['data/'+view+'/group1'][...];assert olddata.shape==(24,300)
  if view=='Plasma_metabolome':
   data=raw[top].T;data=(data-data.mean(0))/data.std(0,ddof=1);features=[keys[x] for x in topids];legacy=topids
   assert set(topids)-set(oldfeatures)=={'M336T626.1'} and set(oldfeatures)-set(topids)=={'M517T467'}
  else:data=olddata.copy();features=oldfeatures;legacy=oldfeatures;assert np.array_equal(data,olddata)
  assert data.shape==(24,300) and np.isfinite(data).all() and len(set(features))==300 and (data.std(0,ddof=1)>0).all()
  assert np.abs(data.mean(0)).max()<1e-6 and np.abs(data.std(0,ddof=1)-1).max()<1e-6
  np.save(O/'inputs'/f'{view}.npy',data);(O/'inputs'/f'{view}_features.json').write_text(json.dumps(features),encoding='utf-8');stored[view]=dict(dtype=str(olddata.dtype),reused_exact=(view!='Plasma_metabolome'))
  qc.append(dict(view=view,N_samples=24,N_features=300,max_abs_mean=float(np.abs(data.mean(0)).max()),max_abs_sd_minus1=float(np.abs(data.std(0,ddof=1)-1).max()),zero_variance=0,missing=0,reused_HDF5_exact=view!='Plasma_metabolome',matrix_SHA256=sha(O/'inputs'/f'{view}.npy')))
  feature_alias.extend(dict(view=view,order_1=j+1,feature_id=f,legacy_feature_id=l) for j,(f,l) in enumerate(zip(features,legacy)))
 names=decode(h['training_opts'].attrs['names']);oldopts=dict(zip(names,h['training_opts'][...].tolist()));assert oldopts['seed']==20260726 and oldopts['maxiter']==1000
 oldmodel={k:decode(np.atleast_1d(v[...])) for k,v in h['model_options'].items()}
 oldZ=h['expectations/Z/group1'][...];np.save(O/'inputs/old_Z.npy',oldZ)
(O/'inputs/samples.json').write_text(json.dumps(hs),encoding='utf-8')
for n,rs in [('input_validation.tsv',qc),('feature_order_alias.tsv',feature_alias)]:
 with (O/'inputs'/n).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
