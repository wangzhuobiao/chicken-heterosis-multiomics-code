# Article 1 R2 analysis code

This repository candidate contains the analysis code for a three-group chicken multi-omics study. The CW correction branch awaits author acceptance.

## Data access

Sequence accessions and release dates are pending verification. Obtain authorized input matrices and sample metadata from the data custodian. `input_file_contract.tsv` records the input layouts. `metadata/feature_keys_5272.tsv` contains the complete key map for the analytical feature universe. `metadata/plasma_module_membership_1050.tsv` contains all frozen plasma module membership rows.

## Runtime and commands

Use Python 3.12.14 with `requirements.txt` for the repair and phenotype Python steps. Use a separate Python 3.13.9 environment with `requirements-mofa.txt` for MOFA preparation, preflight, training and native HDF5 comparison. `environment.yml` defines the general Python environment. The locally verified R runtime is R 4.6.0. `r_package_versions.tsv` records observed versions and missing packages. R dependencies are installed by the user in advance. Historical annotation databases and pathway membership snapshots are required to claim an exact historical GSEA reconstruction. A fresh download may differ from the historical database.

Copy `config/example.json` to a private `local.json` and enter authorized paths. Configuration paths resolve relative to the command's working directory. Use a fresh output directory per command. The launcher keeps the six module entry points explicit and runs the original scientific calculations in the listed script units. It creates per-task logs. Where the archived analysis reads another frozen stage, supply the relevant directory using the configured input key. For a local rerun use newly verified outputs in place of that directory as shown below.

| Module | Actions and results |
|---|---|
| `01-phenotype` | `population` and `abdominal-fat` reproduce the pre-CW phenotype, direction-aware MPH/BPH, 2000-draw bootstrap and five zero-handling branches. `positive-and-subset` consumes the frozen phenotype and AF tables and reproduces the positive-only mean intervals, pairwise tests and subset estimates. `cw-candidate` requires `--allow-unfrozen-cw` and the authorized private `A1_CW_RECORD_FILE`. It reproduces the pending CW correction, 36-trait BH, AF/CW and CW residual updates, and dependent CI rows. |
| `02-rna` | `deseq-gsea` is the historical source candidate for initial DESeq2 and GSEA. `f1-baseline` uses archived VST and pairwise tables. `bootstrap100` runs the accepted full refit protocol and `audit-existing100` reproduces the stability summaries from stored results. For an independent run provide its new `F1_trial100` directory as `A1_F1_RUN`. |
| `03-metabolomics` | `repair` reconstructs and independently verifies the 5,272 × 24 matrix. Set `A1_REPAIR` to that output directory. `statistics` generates three contrasts, the two BH families and the F1 position classifications. |
| `04-associations` | `main` updates frozen module scores and the 160 linear and 132 cross-omics tests. Set `A1_ASSOCIATIONS` to its output directory. `residual160` then reproduces the independent 160 residual Spearman sensitivity results. Corrected RNA and cecal module scores enter as frozen inputs. |
| `05-mofa-training` | `select` reconstructs the variance-ranked feature list from the old model and repaired matrix. Set `A1_MOFA_SELECTION` to its output directory. `prepare` creates three 24 × 300 input arrays. Set `A1_MOFA_INPUTS` to `prepare/inputs`. `preflight` verifies exact input hashes, versions, backend options and data conversion. `train-once` runs the archived protocol once and saves the native HDF5. Use a Python executable from the MOFA environment in the configuration. |
| `06-mofa-downstream` | `archive-and-phenotype` reads the native archive or regenerated training output through `A1_MOFA_ARCHIVE`, compares to the Stage 7 alignment reference through `A1_MOFA_REFERENCE`, and calculates the 32 phenotype associations. Supply `A1_MOFA_INPUTS` when reading a new training output. |

```sh
python run.py 01-phenotype population --config local.json --output runs/population
python run.py 01-phenotype abdominal-fat --config local.json --output runs/af
python run.py 01-phenotype positive-and-subset --config local.json --output runs/positive
python run.py 02-rna audit-existing100 --config local.json --output runs/f1-audit
python run.py 03-metabolomics repair --config local.json --output runs/repair
# Set A1_REPAIR to runs/repair.
python run.py 03-metabolomics statistics --config local.json --output runs/metabolites
python run.py 04-associations main --config local.json --output runs/associations
# Set A1_ASSOCIATIONS to runs/associations.
python run.py 04-associations residual160 --config local.json --output runs/sensitivity
python run.py 05-mofa-training select --config local.json --output runs/selection
# Set A1_MOFA_SELECTION to runs/selection.
python run.py 05-mofa-training prepare --config local.json --output runs/preparation
# Set A1_MOFA_INPUTS to runs/preparation/inputs and choose the MOFA Python runtime.
python run.py 05-mofa-training preflight --config local.json --output runs/preflight
python run.py 05-mofa-training train-once --config local.json --output runs/training
# Set A1_MOFA_ARCHIVE to runs/training and retain the Stage 7 reference configuration.
python run.py 06-mofa-downstream archive-and-phenotype --config local.json --output runs/mofa-downstream
```

The optional `bootstrap100` action uses the same command pattern and performs repeated model fits. Its input gate requires the approved RNA internal identity record and the same raw counts. The `f1-baseline` and `deseq-gsea` actions require `RNA_CORE_OUTDIR` or `ANNOT_FILE`, respectively, plus `META_FILE`. For the R3 candidate the explicit invocation is `python run.py 01-phenotype cw-candidate --allow-unfrozen-cw --config local.json --output runs/cw-candidate`. The historical intermediate CSVs consumed by later stages must retain their original row order, IDs, fields and test family. This workflow is a replay from the archived inputs. For initial outputs generated from the raw phenotype XLSX, the archive also provides cleaned matrices and the fixed selection key used by the later steps.

## Evidence and release status

See `analysis_code_source_mapping.tsv` for original script hashes, run evidence, linked figures and tables, and per-step replay status. `reproduction_QA.md` records the fresh checks and their limitations. `excluded_scripts_manifest.tsv` describes the bounded source inventory and why excluded code is not in the release. `additional_source_requirements.tsv` lists the minimum separate cecal and figure dependencies. These supplemental source leads have individual evidence statuses. R3 artwork is still a candidate and Figure S8 remains unlocated.

The release awaits CW acceptance, author approval of the license and data release authorization. See `LICENSE` for the review status.

## Private inputs and outputs

- `A1_MOFA_INPUTS/samples.json` is the authorized ordered sample list. The data custodian supplies `A1_MOFA_IDENTITY_FILE` as JSON with the `samples_sha256` key. The public backend configuration uses `sample_01` through `sample_24` as order placeholders. Training still uses the authorized labels and order locally.
- `A1_RNA_LABEL_MAP` is an optional private TSV with the header `sample_id`, `original_sample_id`, `plot_label`. Provide it to reproduce archived display aliases. With the path omitted, display labels come from the local input metadata. The aliases do not determine sample inclusion or statistical group membership.
- `A1_CW_RECORD_FILE` is a private TSV with one data row under `wing_id`, `group`, `expected_cw_g`. The CW task verifies the row and expected value before applying the correction.
- **Run outputs can contain original sample identifiers. Keep output directories, logs, HDF5/NPY files and private input/configuration files out of public commits and deposits.** Archive only the reviewed release files. `.gitignore` includes common local run and configuration paths.
