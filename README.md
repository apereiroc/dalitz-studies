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

Install qft++
```bash
bash install_qft.sh
```

You probably have to create a symlink
```bash
ln -s $CONDA_PREFIX/external/qft/lib/libqft++* /usr/local/lib
```
# Usage instructions

Here are some examples for generating, fitting and plotting signal and normalisation toy events

## Generation

The first step is to obtain a normalisation set for computing the normalisation integral during the fit. This integral is achieved by MC (the only practical solution for 4-body problems).
You can modify the number of events within the main `Snakefile`

```bash
scons btokpikpi_CovGen_norm

snakemake tuples/normalisation/Bs/normalisation_set.root -j7
```

(scale the number of jobs to your available resources)

Then generate the signal set
```bash
scons btokpikpi_CovGen_sig get_hmax_Bs -j2

snakemake tuples/signal/Bs/signal_set.root -j7
```








