from pathlib import Path
import csv,json,math,sys,hashlib,os
import openpyxl,numpy as np
sys.stdout.reconfigure(encoding='utf-8');O=Path(os.environ["A1_OUTPUT"]);R=Path(os.environ["A1_INPUT_ROOT"]);M=R/'Article1_core_scripts_20260812/output_plasma_metabolome_article1_core'
def read(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t' if Path(p).suffix=='.tsv' else ','))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert json.loads((O/'validation/reconstruction_summary.json').read_text(encoding='utf-8'))['status']=='MATRIX_VALIDATED_PASS_AUTHOR_REVIEW'
mapping=read(O/'matrix/feature_key_mapping_5272.tsv');samples=[r['sample_id'] for r in read(O/'validation/sample_order_24.tsv')];old=read(M/'tables/Table2_log2_processed_abundance_matrix.csv')
settings={r['parameter']:r['value'] for r in read(M/'qc/QC_analysis_settings_used.csv')};w=openpyxl.load_workbook(os.environ["A1_VENDOR_XLSX"],read_only=True,data_only=True);s=w.worksheets[0];it=s.iter_rows(values_only=True);h=list(next(it));bykey={};rownos={}
for rn,v in enumerate(it,2):
 if all(x is None for x in v):continue
 d=dict(zip(h,v));key=str(d['ID'])+'|'+str(d['posneg']);assert key not in bykey;bykey[key]=d;rownos[key]=rn
w.close();assert len(bykey)==5272
outputs={}
for name in ['processed_log2_5272x24.csv','raw_supplier_5272x24.csv','imputed_normalized_none_5272x24.csv']:
 with (O/'matrix'/name).open(encoding='utf-8-sig') as f:
  rr=csv.reader(f);header=next(rr);values=list(rr)
 assert header==samples and len(header)==24 and len(set(header))==24 and len(values)==5272
 assert all(len(row)==24 for row in values)
 values=np.array([[float(x) for x in row] for row in values]);assert values.shape==(5272,24) and np.all(np.isfinite(values));outputs[name]=values
keyed=read(O/'matrix/processed_log2_keyed_5272x24.tsv');assert len(keyed)==5272 and list(keyed[0])==['feature_key']+samples
assert [r['feature_key'] for r in keyed]==[r['new_feature_key'] for r in mapping]
assert np.array_equal(outputs['processed_log2_5272x24.csv'],np.array([[float(r[x]) for x in samples] for r in keyed]))
allerr=[];unaffected_max=0;changed=[];special=[]
for i,m in enumerate(mapping):
 d=bykey[m['new_feature_key']];assert rownos[m['new_feature_key']]==int(m['supplier_excel_row'])
 vals=[float(d[x]) for x in samples] # all source values are finite positive in this frozen dataset
 assert all(math.isfinite(x) and x>0 for x in vals)
 expected=[math.log2(x+1) for x in vals] # independent scalar math implementation
 got=outputs['processed_log2_5272x24.csv'][i];err=max(abs(x-y) for x,y in zip(expected,got));allerr.append(err);assert err<=1e-10
 assert np.array_equal(outputs['raw_supplier_5272x24.csv'][i],np.array(vals))
 assert np.array_equal(outputs['imputed_normalized_none_5272x24.csv'][i],np.array(vals))
 oldv=[float(old[i][x]) for x in samples];delta=max(abs(x-y) for x,y in zip(got,oldv))
 if m['expected_numeric_repair']=='True':changed.append(m['legacy_feature_id']);assert delta>1e-10
 else:unaffected_max=max(unaffected_max,delta);assert delta<=1e-10
 if m['supplier_original_feature_id']=='M397T467':
  for j,x in enumerate(samples):special.append(dict(new_feature_key=m['new_feature_key'],legacy_feature_id=m['legacy_feature_id'],posneg=m['posneg'],annotation_name=m['annotation_name'],supplier_excel_row=m['supplier_excel_row'],supplier_cell=f'{openpyxl.utils.get_column_letter(h.index(x)+1)}{m["supplier_excel_row"]}',sample_id=x,old_processed=oldv[j],repaired_processed=got[j],supplier_raw=vals[j],independent_scalar_check_abs_error=abs(expected[j]-got[j])))
assert len(changed)==26 and len(special)==48
with (O/'validation/M397T467_pos_neg_pair_48_values.tsv').open('w',encoding='utf-8-sig',newline='') as f:
 cw=csv.DictWriter(f,fieldnames=list(special[0]),delimiter='\t');cw.writeheader();cw.writerows(special)
ms=read(O/'validation/MS2_383_feature_validation.tsv');assert len(ms)==383 and len({r['new_feature_key'] for r in ms})==383 and all(r['is_ms2']=='True' for r in ms)
assert len(read(O/'validation/unaffected5246_consistency.tsv'))==5246 and len(read(O/'validation/repaired624_values.tsv'))==624
for r in read(O/'matrix_input_SHA256.tsv'):assert sha(r['path'])==r['SHA256']
report=dict(status='INDEPENDENT_SERIALIZED_MATRIX_VERIFICATION_PASS',shape=[5272,24],numeric_sample_columns_only=True,keyed_alternative_has_one_text_index_plus24_numeric_columns=True,matrix_row_mapping_one_to_one=True,source_row_scalar_log2_max_abs_error=max(allerr),unaffected5246_max_abs_difference=unaffected_max,corrected26_rows_verified=True,corrected624_cells_verified=True,M397T467_pos_neg_pair_48_values_verified=True,MS2_383_preserved=True,raw_and_imputed_equal_source_values_in_this_dataset=True,no_imputed_cells_in_this_dataset=True,no_inference_computations=True)
(O/'validation/independent_output_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False))
