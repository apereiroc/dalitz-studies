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

The first step is to obtain a normalisation set to compute the normalisation integral during the fit. This integral is achieved by MC (the only practical solution for 4-body problems).
You can modify the number of events within the main `Snakefile`

```bash
snakemake tuples/normalisation/Bs/normalisation_set.root -j7
```

```bash
snakemake tuples/normalisation/Du/normalisation_set.root -j7
```

(scale the number of jobs to your available resources)

Then generate the signal set
```bash
snakemake tuples/signal/Bs/signal_set.root -j7
```

```bash
snakemake tuples/signal/Du/signal_set.root -j7
```

(Optional) check the number of events
```bash
python -c "import uproot; print('Sig:  ', uproot.open('tuples/signal/Bs/signal_set.root')['fitTree'].num_entries)"
python -c "import uproot; print('Norm: ', uproot.open('tuples/normalisation/Bs/normalisation_set.root')['fitTree'].num_entries)"
```

```bash
python -c "import uproot; print('Sig:  ', uproot.open('tuples/signal/Du/signal_set.root')['fitTree'].num_entries)"
python -c "import uproot; print('Norm: ', uproot.open('tuples/normalisation/Du/normalisation_set.root')['fitTree'].num_entries)"
```

## Fit

Let's fit our model to the signal data generated in the previous step, using the normalisation events to normalise the PDF.
The -2logL function is computed in parallel, with the number of jobs read by the Snakemake workflow

```bash
snakemake output/Bs/pars.json -j7
```

```bash
snakemake output/Du/pars.json -j7
```

## Plot 

Let's now plot the fit results 

```bash
snakemake output/Bs/plot.pdf -j7
```

```bash
snakemake output/Du/plot.pdf -j7
```

## Scan

Wanna run a 1D likelihood scan for some specific parameter? You could run

```bash
snakemake output/{meson}/scan/plot-{param}.pdf -j7
```

Example:

```bash
snakemake output/Du/scan/plot-arg_VV_D.pdf -j7
open output/Du/scan/plot-arg_VV_D.pdf
```

