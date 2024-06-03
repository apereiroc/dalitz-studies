# dalitz-studies

Little framework designed to check for differences in dynamics between $B$ and $D$ mesons
due to the phase space size, as well as to run some studies of interest at physics level.

# Install instructions

Create a `conda` environment, activate it and install `Snakemake`

```bash
conda env create --file=environment.yml
conda activate dalitzstudies
python -m pip install snakemake plottergeist
```

Install qft++ inside `$CONDA_PREFIX`
```bash
bash install_qft.sh
```

# Usage instructions

Here are some examples for generating, fitting and plotting signal and normalisation toy events

## Compilation

The first step is to compile all `C++` files. 

You can do it with `SCons` 
```bash
scons -j$(python -c "import os; print(os.cpu_count())")
```

Or with `CMake`
```bash
cmake -B build -S . && cmake --build build -- -j$(python -c "import os; print(os.cpu_count())")
```


## Generation

### Normalisation events
The first step is to obtain a normalisation set to compute the normalisation integral during the fit, required by the maximum likelihood principle

$$\int  d\mathbf{x} \\  \mathcal{P}(\mathbf{x}) = 1$$

This integral is achieved by Monte Carlo, being the only practical solution for 4-body problems

$$\int  d\mathbf{x} \\  \mathcal{P}(\mathbf{x}) \approx \sum_{i=1}^{N} \mathcal{P} (\mathbf{x}_i)$$

You can modify the number of events in the main `Snakefile`. Also you could need to scale the number of jobs to your available resources.

#### - $B^0_s$ normalisation events
```bash
snakemake tuples/normalisation/Bs/normalisation_set.root -j7
```

#### - $D^0$ normalisation events
```bash
snakemake tuples/normalisation/Du/normalisation_set.root -j7
```


### Signal events
Then generate the Monte Carlo signal set using the accept-reject method.

#### - $B^0_s$ signal events
```bash
snakemake tuples/signal/Bs/signal_set.root -j7
```

#### - $D^0$ signal events
```bash
snakemake tuples/signal/Du/signal_set.root -j7
```

(Optional) check the number of events

#### - For $B^0_s$
```bash
python -c "import uproot; print('Sig:  ', uproot.open('tuples/signal/Bs/signal_set.root')['fitTree'].num_entries)"
python -c "import uproot; print('Norm: ', uproot.open('tuples/normalisation/Bs/normalisation_set.root')['fitTree'].num_entries)"
```

#### - For $D^0$
```bash
python -c "import uproot; print('Sig:  ', uproot.open('tuples/signal/Du/signal_set.root')['fitTree'].num_entries)"
python -c "import uproot; print('Norm: ', uproot.open('tuples/normalisation/Du/normalisation_set.root')['fitTree'].num_entries)"
```

## Fit

Let's fit our model to the signal data generated in the previous step, using the normalisation events to normalise the PDF. The likelihood function,

$$-2\log\mathcal{L} = -2 \\ \sum_{i=1}^{\rm Data} \log \mathcal{P}(\mathbf{x}_i)$$

as well as the normalisation integral

$$\int  d\mathbf{x} \\  \mathcal{P}(\mathbf{x})$$

are computed in parallel, with the number of jobs read by the Snakemake workflow

#### - $B^0_s$ fit
```bash
snakemake output/Bs/pars.json -j7
```

#### - $D^0$ fit
```bash
snakemake output/Du/pars.json -j7
```

## Plot 

Let's now plot the fit results 

#### - $B^0_s$ plot
```bash
snakemake output/Bs/plot.pdf -j7
```

#### - $D^0$ plot
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

