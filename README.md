# dalitz-studies

Little framework designed to check for differences between $B$ and $D$ mesons
due to the phase space size

# Install instructions

Create a conda environment, activate it and install snakemake

```bash
conda env create --file=environment.yml
conda activate dalitzstudies
python -m pip install snakemake
```

Install qft++ inside `$CONDA_PREFIX`
```bash
bash install_qft.sh
```

# Usage instructions

Here are some examples for generating, fitting and plotting signal and normalisation toy events

## Compilation

The first step is to compile all c++ files. You can do it with SCons 
```bash
scons -j$(python -c "import os; print(os.cpu_count())")
```

TODO: add CMake 

## Generation

The first step is to obtain a normalisation set for computing the normalisation integral during the fit. This integral is achieved by MC (the only practical solution for 4-body problems).
You can modify the number of events within the main `Snakefile`

```bash
snakemake tuples/normalisation/Bs/normalisation_set.root -j7
```

(scale the number of jobs to your available resources)

Then generate the signal set
```bash
snakemake tuples/signal/Bs/signal_set.root -j7
```

(Optional) check the number of events
```bash
python -c "import uproot; print('Sig:', uproot.open('tuples/signal/Bs/signal_set.root')['fitTree'].num_entries)"
python -c "import uproot; print('Norm:', uproot.open('tuples/normalisation/Bs/normalisation_set.root')['fitTree'].num_entries)"
```

## Fit
TODO



