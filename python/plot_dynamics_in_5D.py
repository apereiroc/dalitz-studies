import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.interpolate as interpol
import matplotlib as mpl
import argparse
import uproot
import plottergeist as pg

p = argparse.ArgumentParser(description="pipas")
p.add_argument("--input-file", required=True, type=str)
p.add_argument("--output-fig", required=True, type=str)
p.add_argument("--bins", required=True, type=int)
args = vars(p.parse_args())

plt.rcParams.update({"text.latex.preamble": ""})

f = uproot.open(args["input_file"])
vars_tree = f["variables"]
vars = vars_tree.arrays(library="pd").to_numpy()

VV_tree = f["VV"]
df_VV = VV_tree.arrays(library="pd")

fig, ax, axleg = pg.make_plot(ndim=5, pull=False, figsize=(10, 6))

units = [
    "",
    "",
    r"rad",
    r"GeV/c$^2$",
    r"GeV/c$^2$",
]

xlabels = [
    r"$\cos \theta_1$",
    r"$\cos \theta_2$",
    r"$\varphi$",
    r"m($K^+ \pi^-$)",
    r"m($K^- \pi^+$)",
]

bins = args["bins"]

density = True

for i in range(5):
    var = vars[:, i]

    hdata = pg.make_hist(var, bins=bins)

    ax[i].hist(
        var,
        weights=df_VV["spin_S"],
        align="mid",
        linewidth=1.5,
        histtype="step",
        bins=bins,
        density=density,
        label="Spin L=0",
    )

    ax[i].hist(
        var,
        weights=df_VV["spin_P"],
        align="mid",
        linewidth=1.5,
        histtype="step",
        bins=bins,
        density=density,
        label="Spin L=1",
    )

    ax[i].hist(
        var,
        weights=df_VV["spin_D"],
        align="mid",
        linewidth=1.5,
        histtype="step",
        bins=bins,
        density=density,
        label="Spin L=2",
    )

    ax[i].hist(
        var,
        weights=df_VV["bf_S"],
        align="mid",
        linewidth=1.5,
        histtype="step",
        bins=bins,
        density=density,
        label="Barrier factor L=0",
    )

    ax[i].hist(
        var,
        weights=df_VV["bf_P"],
        align="mid",
        linewidth=1.5,
        histtype="step",
        bins=bins,
        density=density,
        label="Barrier factor L=1",
    )

    ax[i].hist(
        var,
        weights=df_VV["bf_D"],
        align="mid",
        linewidth=1.5,
        histtype="step",
        bins=bins,
        density=density,
        label="Barrier factor L=2",
    )

    ax[i].set_ylim(ymin=0)
    ax[i].set_ylabel(
        "Candidates / {:.2f} {}".format(hdata.bins[1] - hdata.bins[0], units[i]),
    )

    u = ""
    if units[i] != "":
        u = "[{}]".format(units[i])

    ax[i].set_xlabel(r"{} {}".format(xlabels[i], u), fontsize=12)

ymax_cos = max(ax[0].get_ylim()[1], ax[1].get_ylim()[1])
ymax_mass = max(ax[3].get_ylim()[1], ax[4].get_ylim()[1])

ax[0].set_ylim(ymax=ymax_cos)
ax[1].set_ylim(ymax=ymax_cos)
ax[3].set_ylim(ymax=ymax_mass)
ax[4].set_ylim(ymax=ymax_mass)

handles, labels = ax[0].get_legend_handles_labels()

axleg.axis("off")
axleg.legend(handles, labels, loc="center")

# fig.tight_layout()

fig.savefig(args["output_fig"])
