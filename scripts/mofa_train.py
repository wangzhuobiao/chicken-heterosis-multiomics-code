"""ONE author-approved archive repair run. No retry; original calculation path preserved."""
import sys
import os
import importlib.metadata as metadata
sys.dont_write_bytecode=True
from pathlib import Path
import csv,json,hashlib,time,traceback
import numpy as np,h5py
from mofapy2.run.entry_point import entry_point
O=Path(os.environ["A1_OUTPUT"]); inputs=Path(os.environ["A1_MOFA_INPUTS"]); cfgdir=Path(__file__).resolve().parents[1]/"config"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for requirement in (cfgdir.parent/'requirements-mofa.txt').read_text().splitlines():
 name, version= requirement.split('==')
 assert metadata.version(name)==version, 'Version mismatch: '+name
expected=json.loads((cfgdir/'mofa_input_sha256.json').read_text())
identity=json.loads(Path(os.environ['A1_MOFA_IDENTITY_FILE']).read_text(encoding='utf-8-sig'))
assert sha(inputs/'samples.json')==identity['samples_sha256'], 'Private sample order does not match the approved identity file'
for name, value in expected.items():
 assert sha(inputs/name)==value, 'Frozen training input changed: '+name
for folder in ['model','results','validation']:(O/folder).mkdir(exist_ok=True)
assert not (O/'TRAINING_ATTEMPT.json').exists(),'No automatic second attempt'
outfile=O/'model/Article1_MOFA2_Stage7_1_native.hdf5';assert not outfile.exists()
cfg=json.loads((cfgdir/'training_config.json').read_text());samples=json.loads((inputs/'samples.json').read_text());features=[json.loads((inputs/f'{v}_features.json').read_text()) for v in cfg['views']];data=[[np.load(inputs/f'{v}.npy')] for v in cfg['views']]
try:
 ent=entry_point();ent.set_data_options(**cfg['data']);ent.set_data_matrix(data,likelihoods=cfg['likelihoods'],views_names=cfg['views'],groups_names=['group1'],samples_names=[samples],features_names=features);ent.set_model_options(**cfg['model']);ent.set_train_options(**cfg['training']);assert ent.train_opts['stochastic'] is False
 resolved=json.loads(json.dumps(dict(data=ent.data_opts,model=ent.model_opts,training=ent.train_opts),default=str))
 assert len(samples)==24 and len(set(samples))==24 and resolved['data']['samples_names']==[samples]
 # Redact only the comparison copy. Actual backend data retain the private names and order.
 resolved['data']['samples_names']=[[f'sample_{i:02d}' for i in range(1,25)]]
 assert resolved==json.loads((cfgdir/'resolved_backend_options.json').read_text())
 for i, view in enumerate(cfg['views']):
  np.save(O/'validation'/f'{view}_expected_backend_data.npy',ent.data[i])
 if os.environ.get('A1_MOFA_PREFLIGHT')=='1':
  print('Frozen inputs, dependency versions and resolved options verified. Fit not started.')
  sys.exit(0)
 ent.build()
 # Exclusive create prevents a second invocation, including after failed save.
 with (O/'TRAINING_ATTEMPT.json').open('x') as f:json.dump(dict(attempts=1,seed=20260726,start_time=time.time(),config_SHA256=sha(cfgdir/'training_config.json'),tolerances_SHA256=sha(cfgdir/'tolerances_frozen.json')),f)
 ent.run()
 ent.save(str(outfile),save_data=True,save_parameters=False)
 # No manual HDF5 writing; native read-back check before completion marker.
 with h5py.File(outfile,'r') as h:
  for n in ['data','expectations','features','samples','views','groups','intercepts','model_options','training_opts','training_stats','variance_explained']:assert n in h,n
 # Read-only objective evaluation; NOT inserted into or used to patch native HDF5 stats.
 final_elbo=ent.model.calculateELBO()
 (O/'results/final_objective_readonly.json').write_text(json.dumps(dict(source='read-only calculateELBO after native save; no variational update',terms={k:float(v) for k,v in final_elbo.items()},native_stats_caveat='mofapy2 iterate truncates stats with [:i] when converged; final convergence iteration remains in stdout'),indent=2))
 (O/'TRAINING_COMPLETED.json').write_text(json.dumps(dict(status='NATIVE_SAVE_RETURNED_BASIC_READBACK_PASS',attempts=1,seed=20260726,completed_time=time.time(),model_SHA256=sha(outfile),independent_validation_pending=True),indent=2))
except Exception:
 (O/'FAILURE_training.txt').write_text(traceback.format_exc(),encoding='utf-8');raise
