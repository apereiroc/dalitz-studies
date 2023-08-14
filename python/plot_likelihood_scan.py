import uproot
import argparse
import json
import numpy as np
import plottergeist

plottergeist.fig_creator.FIGSIZE = None

p = argparse.ArgumentParser()
p.add_argument('--input-file')
p.add_argument('--param')
p.add_argument('--input-pars')
p.add_argument('--output-fig')
args = vars(p.parse_args())


with open(args['input_pars'], 'r') as file:
    pars = json.load(file)


f = uproot.open(args["input_file"])
plot_tree = f["plot"]
data = plot_tree.arrays(library="pd").to_numpy()

x = data[:,0]
y = data[:,1]

y -= np.min(y)


fig, ax = plottergeist.make_plot(ndim=1, pull=False)

ax.plot(x, y, linewidth=2.0)
ax.set_xlabel(pars[args['param']]['latex'])
ax.set_ylabel(r'$-2\log \mathcal{L}$')

ax.set_title('Likelihood scan with {} fits'.format(len(x)))

# Calculate minimum x,y-values
y_min = np.min(y)
x_min = x[(y == y_min)][0]
x_true = pars[args['param']]['val']
y_true = np.interp(x_true, x, y)


ax.axvline(x=x_min, color='#ff2d55', linestyle='--', label='Min x={:.3f} -2logL={:.3f}'.format(x_min, y_min))
ax.axvline(x=x_true, color='#ff9500', linestyle='--', label='Gen x={:.3f} -2logL={:.3f}'.format(x_true, y_true))


ax.legend(loc='best')

fig.savefig(args['output_fig'])
