# epitopes-world-pred

CAMAP-based MHC-I peptide presentation predictions for viral proteomes, behind
[epitopes.world](https://www.epitopes.world).

The repository holds two independent pipelines. Both score peptides with CAMAP and join
NetMHCpan binding predictions, but they answer different questions:

| Folder | Question | Output |
|---|---|---|
| `vaccine/` | Which peptides of the SARS-CoV reference proteomes are likely presented? | vaccine candidates |
| `gisaid/` | How do those predictions vary across GISAID SARS-CoV-2 strains? | the portal's strain tables |

`vaccine/` is the analysis code for the preprint *Codon arrangement modulates MHC-I peptides
presentation: implications for a SARS-CoV-2 peptide-based vaccine* (Daouda, Dumont-Lagace,
Feghaly, Villani; bioRxiv 2021.02.04.429819).

## Layout

```
data/         inputs (not tracked): GISAID download, genome records, peptide and allele lists
validation/   known epitopes: SARS-CoV-1, HIV-1, allele target list (tracked)
output/       everything generated (not tracked)
gisaid/       strain pipeline: step1-step10, job1-job7
vaccine/      reference-proteome pipeline: tools.py, model_tools.py, bootstrap.py
```

`output/` subfolders: `sequences/` and `scored/` (per-protein contexts), `pickles/`,
`tables/`, `netmhc/gisaid/` and `netmhc/sarscov/` (per-allele NetMHCpan predictions, one folder
per pipeline), `predictions/`, `figures/`, `pytorch-sgd/` (the CAMAP model).

## Running

**All scripts expect the repository root as the working directory**, including the job
scripts, which submit with `--chdir $PWD`. Submit them as `gisaid/job1.parse.sh`, not from
inside `gisaid/`.

## Dependencies

`camaptools` and `camap` (the CAMAP_Tools and CAMAP repositories), `pyGeno`, Biopython,
pandas, numpy, `mpi4py`, `parse`, plus NetMHCpan 4.0a and MUSCLE as external tools.
