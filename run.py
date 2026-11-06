"""Run one analysis task in a fresh output directory using authorized local inputs."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

PACKAGE = Path(__file__).resolve().parent
TASKS = {
    '01-phenotype': {
        'population': ['phenotype_population.R'],
        'abdominal-fat': ['phenotype_abdominal_fat.R'],
        'positive-and-subset': ['phenotype_positive_and_subset.py'],
        'cw-correction': ['phenotype_CW_correction.R', 'phenotype_CW_intervals.py'],
    },
    '02-rna': {
        'deseq-gsea': ['rna_deseq_gsea.R'],
        'f1-baseline': ['rna_f1_baseline.R'],
        'bootstrap100': ['rna_f1_bootstrap100.R'],
        'audit-existing100': ['rna_f1_stability_audit.R'],
    },
    '03-metabolomics': {
        'repair': ['metabolite_repair.py', 'metabolite_verify.py'],
        'statistics': ['metabolite_statistics.R'],
    },
    '04-associations': {
        'main': ['module_and_associations.R'],
        'residual160': ['residual_spearman160.R'],
    },
    '05-mofa-training': {
        'select': ['mofa_select_top300.R'],
        'prepare': ['mofa_prepare_inputs.py'],
        'preflight': ['mofa_train.py'],
        'train-once': ['mofa_train.py'],
    },
    '06-mofa-downstream': {'archive-and-phenotype': ['mofa_native_downstream.py', 'mofa_phenotype32.R']},
}
REQUIRED = {
    'population': ['A1_PHENOTYPE_XLSX'], 'abdominal-fat': ['A1_PHENOTYPE_XLSX'],
    'positive-and-subset': ['A1_INPUT_ROOT','A1_PHENOTYPE_XLSX'],
    'cw-correction': ['A1_INPUT_ROOT','A1_CW_RECORD_FILE'],
    'deseq-gsea': ['A1_RNA_COUNTS','META_FILE','ANNOT_FILE'],
    'f1-baseline': ['RNA_CORE_OUTDIR','META_FILE'],
    'bootstrap100': ['A1_INPUT_ROOT','A1_RNA_COUNTS'],
    'audit-existing100': ['A1_INPUT_ROOT','A1_RNA_COUNTS','A1_F1_RUN'],
    'repair': ['A1_INPUT_ROOT','A1_VENDOR_XLSX'],
    'statistics': ['A1_INPUT_ROOT','A1_REPAIR'],
    'main': ['A1_INPUT_ROOT','A1_REPAIR'],
    'residual160': ['A1_INPUT_ROOT','A1_ASSOCIATIONS'],
    'select': ['A1_INPUT_ROOT','A1_REPAIR'],
    'prepare': ['A1_INPUT_ROOT','A1_REPAIR','A1_MOFA_SELECTION'],
    'train-once': ['A1_MOFA_INPUTS','A1_MOFA_IDENTITY_FILE'], 'preflight': ['A1_MOFA_INPUTS','A1_MOFA_IDENTITY_FILE'],
    'archive-and-phenotype': ['A1_INPUT_ROOT','A1_MOFA_ARCHIVE','A1_MOFA_REFERENCE'],
}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('module', choices=TASKS)
    parser.add_argument('action')
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.action not in TASKS[args.module]:
        parser.error('Actions: '+', '.join(TASKS[args.module]))
    config = json.loads(args.config.read_text(encoding='utf-8-sig'))
    env = os.environ.copy()
    env.update({k: str(v) for k,v in config['inputs'].items() if v})
    optional_paths = {'archive-and-phenotype': ['A1_MOFA_INPUTS']}
    input_keys = REQUIRED[args.action] + [key for key in optional_paths.get(args.action, []) if env.get(key)]
    for key in input_keys:
        if not env.get(key):
            parser.error(f'Missing authorized input configuration: {key}')
        input_path = Path(env[key]).expanduser().resolve()
        if not input_path.exists():
            parser.error(f'Input unavailable: {key} ({input_path})')
        env[key] = str(input_path)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env.update(A1_OUTPUT=str(output), A1_CODE=str(PACKAGE/'scripts'),
               PHENO_OUTDIR=str(output), AF_OUTDIR=str(output),
               RNA_OUTDIR=str(output), RNA_F1_OUTDIR=str(output),
               PYTHONDONTWRITEBYTECODE='1')
    env['A1_MOFA_PREFLIGHT']='1' if args.action=='preflight' else '0'
    if 'A1_RNA_COUNTS' in env:
        env['COUNT_FILE'] = env['A1_RNA_COUNTS']
    if 'A1_PHENOTYPE_XLSX' in env:
        env['PHENO_FILE'] = env['A1_PHENOTYPE_XLSX']
    r_launcher = ('if (.Platform$OS.type=="windows") '
                  'Sys.setlocale("LC_CTYPE","English_United States.utf8"); '
                  'source(Sys.getenv("A1_SCRIPT"), encoding="UTF-8")')
    for name in TASKS[args.module][args.action]:
        script = PACKAGE/'scripts'/name
        env['A1_SCRIPT'] = str(script)
        if script.suffix == '.R':
            cmd = [config.get('rscript','Rscript'),'--vanilla','-e',r_launcher,str(output)]
        else:
            cmd = [config.get('python',sys.executable),'-X','utf8','-B',str(script)]
        log = output/(name+'.log')
        if args.action=='train-once':
            (output/'logs').mkdir(exist_ok=True)
            log = output/'logs/training_once.log'
        print(f'{name} -> {log}', flush=True)
        with log.open('w',encoding='utf-8') as f:
            subprocess.run(cmd,env=env,cwd=output,stdout=f,stderr=subprocess.STDOUT,check=True)

if __name__ == '__main__':
    main()
