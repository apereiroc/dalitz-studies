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
