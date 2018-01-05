# SolutionFFBench

SolutionFFBench is a streamlined version of https://github.com/choderalab/LiquidBenchmark
that has been updated to support modern miniconda packaging, python 3.6, Snakemake,
and the SMIRNOFF forcefield.  


## Instructions

0.  Install miniconda if you haven't already:

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

1.  Install requirements, including python 3.6.  Note that on some systems, this
may not 100% automated due to linker issues with packmol, gfort, and openmm.

```bash
conda env create -f environment.yml
source activate benchmark
pip install -i https://pypi.anaconda.org/OpenEye/simple OpenEye-toolkits
# Depending on your system, you may need to do some other manual installations at this point.
# In particular, I ended up having to manually delete openmm and manually build :(
```

2.  Activate environment and export environment variables:

```bash
source activate benchmark
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

Cluster based parallelization is also available via, e.g., `snakemake --cluster="qsub ......"`
The details of the pipeline are the rules in `Snakefile`.
