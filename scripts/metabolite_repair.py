"""Rebuild approved supplier rows; no inference, classifications, scoring or fitting."""
from pathlib import Path
from collections import Counter,defaultdict
import csv,json,hashlib,traceback,sys,re,math,time
import os
import numpy as np
import openpyxl
from openpyxl.utils import get_column_letter

O=Path(os.environ["A1_OUTPUT"]);R=Path(os.environ["A1_INPUT_ROOT"]);PRE=R/'09_output/Article1_R2_Stage5_20261029';M=R/'Article1_core_scripts_20260812/output_plasma_metabolome_article1_core'
TOL=1e-10
sys.stdout.reconfigure(encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t' if Path(p).suffix=='.tsv' else ','))
def write(p,rows,fields=None):
 rows=list(rows);fields=fields or list(rows[0])
 with Path(p).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def num(x):
 if x is None:return float('nan')
 if isinstance(x,bool):raise ValueError('Unexpected boolean abundance')
 if isinstance(x,(int,float)):return float(x)
 s=str(x).replace(',','').replace('，','').strip()
 if s in ['', 'NA','NaN','null','NULL','-','--']:return float('nan')
 return float(s)  # unknown malformed strings stop rather than silently coerce
def is_ms2(x):
 s=str(x).strip().lower();return s=='2' or re.search(r'level\s*2|ms2|二级|2级|^2$',s) is not None
def canonical_level(x):
 if isinstance(x,(int,float)) and float(x).is_integer():return str(int(x))
 return str(x).strip()
def main():
 start=time.perf_counter()
 for n in ['matrix','validation','dependencies','logs']:(O/n).mkdir(exist_ok=True)
 cfgp=M/'qc/QC_analysis_settings_used.csv';cfg={r['parameter']:r['value'] for r in read(cfgp)}
 supplied=Path(os.environ["A1_VENDOR_XLSX"]);matp=M/'tables/Table2_log2_processed_abundance_matrix.csv';annp=M/'tables/Table1_feature_annotation_raw.csv';metap=M/'qc/QC_sample_metadata_used.csv';qcp=M/'qc/QC_feature_detection_missingness.csv'
 conflictpath=PRE/'B_metabolomics/confirmed_input_conflict.json';diagpath=PRE/'B_metabolomics/affected_cell_diagnostic_NOT_ANALYSIS_INPUT.tsv'
 inputs=[supplied,matp,annp,metap,qcp,cfgp,conflictpath,diagpath]
 # The exact archived input hashes are published without private filenames.
 expected=json.loads((Path(__file__).resolve().parents[1]/'config/repair_input_sha256.json').read_text())
 assert [sha(p) for p in inputs]==expected, 'Frozen repair input hash changed'
 write(O/'matrix_input_SHA256.tsv',[dict(path=str(p),bytes=p.stat().st_size,SHA256=sha(p)) for p in inputs])
 required={'ZERO_AS_MISSING':'TRUE','DETECTION_RATE_MIN':'0.5','MISSING_RATE_MAX':'0.5','IMPUTE_METHOD':'half_min','NORMALIZATION_METHOD':'none','LOG_TRANSFORM':'TRUE','LOG_OFFSET':'1','EXCLUDE_SAMPLES':''}
 assert all(cfg[k]==v for k,v in required.items()),'Frozen settings differ from approved repair'
 meta=read(metap);samples=[r['sample_id'] for r in meta];group={r['sample_id']:r['group'] for r in meta}
 assert len(samples)==len(set(samples))==24 and Counter(group.values())==Counter({'YY':8,'YX':8,'XX':8})
 assert all(s.startswith(group[s]) for s in samples)
 old=read(matp);ann=read(annp);qc=read(qcp);conf=json.loads(conflictpath.read_text(encoding='utf-8'));badset=set(conf['affected_feature_IDs'])
 assert len(old)==len(ann)==len(qc)==5272 and len(badset)==26
 oldids=[r['feature_id'] for r in old];assert len(set(oldids))==5272
 assert [r['feature_id'] for r in ann]==oldids
 qmap={r['feature_id']:r for r in qc};assert set(qmap)==set(oldids) and all(r['keep_feature']=='TRUE' for r in qc)
 old_headers=list(old[0]);assert [s for s in old_headers if s in samples]==samples
 wb=openpyxl.load_workbook(supplied,read_only=True,data_only=True);ws=wb.worksheets[0];sheet=ws.title;it=ws.iter_rows(values_only=True);headers=list(next(it));source=[]
 assert len(headers)==len(set(headers)) and set(samples)=={h for h in headers if str(h).startswith(('YY','YX','XX'))}
 for rowno,vals in enumerate(it,2):
  if all(v is None for v in vals):continue
  row=dict(zip(headers,vals));assert row['ID'] not in [None,''],'Unnamed supplier feature row'
  rid=str(row['ID']);mode=str(row['posneg']);assert rid==rid.strip() and mode in ['pos','neg'] and '|' not in rid
  row['_excel_row']=rowno;row['_key']=rid+'|'+mode;source.append(row)
 wb.close();assert len(source)==5272 and len({r['_key'] for r in source})==5272,'Composite key collision or dropped supplier row'
 bykey={r['_key']:r for r in source};byid=defaultdict(list)
 for row in source:byid[str(row['ID'])].append(row)
 keys=[];raw=[];mapping=[];annotations=[]
 for i,(legacy,a) in enumerate(zip(old,ann),1):
  assert legacy['feature_id']==a['feature_id'];key=a['feature_id_original']+'|'+a['posneg']
  assert key in bykey,'Supplier composite key missing';src=bykey[key];assert str(src['ID'])==a['ID']==legacy['ID'] and a['posneg']==legacy['posneg']
  # readxl::read_excel defaults trim_ws=TRUE; retain the verbatim supplier
  # string too, while checking this documented reader-level representation.
  supplier_name='' if src['name'] is None else str(src['name'])
  assert supplier_name.strip()==a['name']==legacy['name'],'Non-whitespace annotation/name differs'
  lev=canonical_level(src['level']);assert lev==canonical_level(num(a['annotation_level']))
  ms=is_ms2(lev);assert ms==(a['is_ms2']=='TRUE')==(legacy['is_ms2']=='TRUE')
  keys.append(key);raw.append([num(src[s]) for s in samples])
  mapping.append(dict(matrix_row_1=i,new_feature_key=key,supplier_original_feature_id=str(src['ID']),posneg=a['posneg'],legacy_feature_id=a['feature_id'],supplier_annotation_name_verbatim=supplier_name,annotation_name=a['name'],name_differs_only_by_reader_trim=supplier_name!=a['name'],annotation_level=lev,is_ms2=ms,supplier_workbook_path=str(supplied),supplier_sheet=sheet,supplier_excel_row=src['_excel_row'],historical_first_match_excel_row=byid[str(src['ID'])][0]['_excel_row'],supplier_ID_occurrences=len(byid[str(src['ID'])]),expected_numeric_repair=a['feature_id'] in badset))
  aa=dict(a);aa['new_feature_key']=key;aa['legacy_feature_id']=aa.pop('feature_id');aa['supplier_excel_row']=src['_excel_row'];aa['supplier_sheet']=sheet;annotations.append(aa)
 assert keys==[r['_key'] for r in source],'Supplier feature row order drifted; do not silently reorder'
 assert len(set(keys))==5272 and set(keys)==set(bykey)
 raw=np.array(raw,dtype=float);oldmat=np.array([[num(r[s]) for s in samples] for r in old]);assert raw.shape==oldmat.shape==(5272,24) and np.isfinite(oldmat).all()
 missing=(raw<=0)|np.isnan(raw);assert not np.isinf(raw).any(),'Unexpected infinite raw abundance'
 detection=np.mean((~np.isnan(raw))&(raw>0),axis=1);missing_rate=missing.mean(axis=1);keep=(detection>=.5)&(missing_rate<=.5)
 assert keep.all(),'Own-row frozen QC would drop a feature; cannot force5272'
 floor=np.array([np.min(v[np.isfinite(v)&(v>0)])/2 if np.any(np.isfinite(v)&(v>0)) else 1 for v in raw]);imputed=np.where(missing,floor[:,None],raw);processed=np.log2(imputed+1)
 assert np.isfinite(processed).all() and processed.shape==(5272,24)
 delta=processed-oldmat;absdelta=np.abs(delta);maxerr=absdelta.max(axis=1);actualbad={oldids[i] for i in range(5272) if maxerr[i]>TOL}
 assert actualbad==badset,'Unexpected change set; stop instead of selecting a repair'
 unaffected=[i for i,g in enumerate(oldids) if g not in badset];affected=[i for i,g in enumerate(oldids) if g in badset]
 assert len(unaffected)==5246 and len(affected)==26 and np.all(absdelta[unaffected]<=TOL)
 assert np.count_nonzero(absdelta>TOL)==624
 msidx=[i for i,r in enumerate(mapping) if r['is_ms2']];assert len(msidx)==383 and len({keys[i] for i in msidx})==383
 assert Counter(mapping[i]['posneg'] for i in msidx)==Counter({'pos':269,'neg':114})
 dmap={(r['feature_id'],r['sample_id']):r for r in read(diagpath)};assert len(dmap)==624
 cellrows=[]
 for i in affected:
  for j,s in enumerate(samples):
   prev=dmap[(oldids[i],s)];assert abs(float(prev['frozen_processed'])-oldmat[i,j])<=TOL and abs(float(prev['supplier_own_row_log2_preprocessing_diagnostic'])-processed[i,j])<=TOL
   cellrows.append(dict(legacy_feature_id=oldids[i],new_feature_key=keys[i],annotation_name=mapping[i]['annotation_name'],posneg=mapping[i]['posneg'],supplier_excel_row=mapping[i]['supplier_excel_row'],supplier_sample_column=headers.index(s)+1,supplier_cell=f'{get_column_letter(headers.index(s)+1)}{mapping[i]["supplier_excel_row"]}',sample_id=s,group=group[s],supplier_raw_value=raw[i,j],old_log2_processed=oldmat[i,j],repaired_log2_processed=processed[i,j],repaired_minus_old=delta[i,j],matches_prior_diagnostic=True))
 rowqc=[]
 for i,m in enumerate(mapping):
  rowqc.append(dict(**m,detection_rate=detection[i],missing_rate=missing_rate[i],keep_feature=bool(keep[i]),imputed_cells=int(missing[i].sum()),half_min_fill_value=floor[i],historical_QC_detection_rate=qmap[oldids[i]]['detection_rate'],historical_QC_missing_rate=qmap[oldids[i]]['missing_rate'],max_abs_processed_difference=float(maxerr[i]),changed_cells_above_tolerance=int(np.count_nonzero(absdelta[i]>TOL)),absolute_tolerance=TOL,relative_tolerance=0,validation='REPAIRED_TO_OWN_SUPPLIER_ROW' if oldids[i] in badset else 'UNCHANGED_WITHIN_TOLERANCE'))
 # All gates passed before publishing any usable matrix.
 np.savetxt(O/'matrix/processed_log2_5272x24.csv',processed,delimiter=',',header=','.join(samples),comments='',fmt='%.17g')
 np.savetxt(O/'matrix/raw_supplier_5272x24.csv',raw,delimiter=',',header=','.join(samples),comments='',fmt='%.17g')
 np.savetxt(O/'matrix/imputed_normalized_none_5272x24.csv',imputed,delimiter=',',header=','.join(samples),comments='',fmt='%.17g')
 # Keyed alternative: one text index + exactly24 numeric sample columns.
 with (O/'matrix/processed_log2_keyed_5272x24.tsv').open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['feature_key']+samples)
  for i,key in enumerate(keys):w.writerow([key]+[format(v,'.17g') for v in processed[i]])
 write(O/'matrix/feature_key_mapping_5272.tsv',mapping);write(O/'matrix/feature_annotations_5272.tsv',annotations)
 write(O/'validation/all5272_row_validation.tsv',rowqc);write(O/'validation/unaffected5246_consistency.tsv',[rowqc[i] for i in unaffected]);write(O/'validation/repaired26_feature_summary.tsv',[rowqc[i] for i in affected]);write(O/'validation/repaired624_values.tsv',cellrows)
 write(O/'validation/MS2_383_feature_validation.tsv',[rowqc[i] for i in msidx])
 write(O/'validation/annotation_reader_whitespace_records.tsv',[r for r in mapping if r['name_differs_only_by_reader_trim']])
 duprows=[r for r in rowqc if r['supplier_ID_occurrences']>1];write(O/'validation/cross_ion_duplicate_IDs_52_rows.tsv',duprows)
 samplemap=[dict(matrix_column_1=j+1,sample_id=s,group=group[s],group_rank_YY_YX_XX=['YY','YX','XX'].index(group[s])+1,supplier_column_1=headers.index(s)+1,supplier_column_letter=get_column_letter(headers.index(s)+1),supplier_header=headers[headers.index(s)],legacy_CSV_column_1=old_headers.index(s)+1,metadata_CSV_line=j+2,supplier_sheet=sheet,verified_order=True) for j,s in enumerate(samples)]
 write(O/'validation/sample_order_24.tsv',samplemap)
 settings=dict(scope='approved input reconstruction and read-only dependency tracing only',key='literal supplier_original_feature_id + "|" + posneg',row_order='frozen matrix / supplier row order, verified identical',sample_order=samples,numeric_columns_only=24,zero_as_missing=True,detection_rate_min=.5,missing_rate_max=.5,impute_method='feature-wise positive minimum /2; fallback1 only if no positive values (none retained in that case)',normalization='none',log_transform='log2(imputed+1)',absolute_tolerance=TOL,relative_tolerance=0,seed='not applicable: deterministic; no RNG calls',input_versions='matrix_input_SHA256.tsv',software={'Python':sys.version,'NumPy':np.__version__,'openpyxl':openpyxl.__version__})
 (O/'reconstruction_settings.json').write_text(json.dumps(settings,ensure_ascii=False,indent=2),encoding='utf-8')
 result=dict(status='MATRIX_VALIDATED_PASS_AUTHOR_REVIEW',shape=[5272,24],feature_keys_unique=5272,raw_rows_consumed=5272,unaffected_rows=5246,unaffected_max_abs_difference=float(maxerr[unaffected].max()),repaired_rows=26,repaired_cells=624,MS2_features=383,MS2_affected=[keys[i] for i in msidx if oldids[i] in badset],all_sample_order_checks=True,all_rows_pass_frozen_QC=True,imputed_cells=int(missing.sum()),changed_QC_rows=sum(abs(float(qmap[g]['detection_rate'])-detection[i])>1e-12 or abs(float(qmap[g]['missing_rate'])-missing_rate[i])>1e-12 for i,g in enumerate(oldids)),no_rows_merged=True,no_statistics_run=True,elapsed_seconds=time.perf_counter()-start)
 (O/'validation/reconstruction_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)),encoding='utf-8');print(json.dumps(result,ensure_ascii=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
if __name__=='__main__':
 try:main()
 except Exception as e:
  (O/'STOPPED_RECONSTRUCTION_ERROR.txt').write_text(traceback.format_exc(),encoding='utf-8')
  (O/'TODO_CONFIRM.md').write_text('STOPPED: '+str(e)+'\nDo not use unfinished matrices; source reconciliation required.\n',encoding='utf-8')
  raise
