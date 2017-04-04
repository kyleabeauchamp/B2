# B2: better benchmarking

1.  Install requirements (not 100% working yet, possible issues with packmol gfort and openmm etc.)

```bash
conda env create -f environment.yml
```

2.  Activate environment and export environment variables:

```bash
source activate py35
export PYTHONPATH=$PYTHONPATH:./code/
export PATH=$PATH:./code/
```


3.  Run pipeline via snakemake:

```bash
# Inspect commands to be run via dry-run:
snakemake -np

# Run the pipeline (local mode):
snakemake -p
```

Cluster mode is also available via, e.g., `snakemake --cluster="qsub ......"`
