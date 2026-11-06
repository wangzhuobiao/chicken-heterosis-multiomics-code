# Article 1 analysis scripts

Use Python 3.12 with `requirements.txt`, R 4.6.0 with the packages listed in `r_package_versions.tsv`, and a separate Python 3.13 environment with `requirements-mofa.txt` for MOFA. The RNA enrichment step also uses `clusterProfiler`.

Copy `config/example.json` to `local.json` and set local input paths. Run a task in a fresh output directory.

```sh
python run.py <module> <task> --config local.json --output runs/<task>
```

| Module | Tasks in execution order |
|---|---|
| `01-phenotype` | `population`, `abdominal-fat`, `positive-and-subset`, `cw-correction` |
| `02-rna` | `deseq-gsea`, `f1-baseline`, `bootstrap100`, `audit-existing100` |
| `03-metabolomics` | `repair`, `statistics` |
| `04-associations` | `main`, `residual160` |
| `05-mofa-training` | `select`, `prepare`, `preflight`, `train-once` |
| `06-mofa-downstream` | `archive-and-phenotype` |

Set each downstream input path to the corresponding local result directory. The CW task reads a local TSV with `wing_id`, `group`, and `expected_cw_g` columns. MOFA reads its ordered input files locally and checks the sample list against `samples_sha256` in the local integrity JSON.

**Keep input data, local configurations, logs, and run outputs out of public uploads. Run outputs may contain sample identifiers.**
