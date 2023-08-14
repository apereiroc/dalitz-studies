import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.interpolate as interpol
import matplotlib as mpl
import argparse
import uproot
import plottergeist

p = argparse.ArgumentParser(description="pipas")
p.add_argument('--input-file', required=True, type=str)
p.add_argument('--output-fig', required=True, type=str)
p.add_argument('--bins', required=True, type=int)
args = vars(p.parse_args())



f = uproot.open(args["input_file"])
data_tree = f["data"]
data = data_tree.arrays(library="pd").to_numpy()

fit_tree = f["fit"]
fit = fit_tree.arrays(library="pd").to_numpy()


fig, ax, ax6, axpull = plottergeist.make_plot(ndim=5, pull=True)

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


for i in range(5):

  hdata, hfit, hpull = plottergeist.compare_hist(data=data[:,i], fit=fit[:,i],
                                    density=False,
                                    data_weights=data[:,-1], fit_weights=fit[:,-1],
                                    bins=bins)


  # hdata = make_hist(data[:,i],bins=bins)

  # f = interpol.UnivariateSpline(x=hfit.bins, y=hfit.counts, k=3, s=300)
  # x = np.linspace(hfit.bins[0], hfit.bins[-1], 1000)


  ax[i].errorbar(hdata.bins, hdata.counts, xerr=hdata.xerr, yerr=hdata.yerr, fmt=".", color="black", label="Data")
  # ax[i].plot(x, f(x), color="dodgerblue", label="Model")
  # ax[i].hist(fit[:,i], weights=fit[:,-1]*len(data)*(hdata.bins[1]-hdata.bins[0])/hfit.norm, align="mid", histtype="step", bins=hdata.edges, color="dodgerblue", label="Fit")
  ax[i].hist(fit[:,i], weights=fit[:,-1]*hdata.norm/hfit.norm, align="mid", linewidth=1.5, histtype="step", bins=hdata.edges, label="Fit")
  # ax[i].hist(fit[:,i], weights=fit[:,-2]*hdata.norm/hfit.norm, align="mid", linewidth=0.5, histtype="step", bins=hdata.edges, color="#bdfe86", label="SS")
  # ax[i].hist(fit[:,i], weights=fit[:,-3]*hdata.norm/hfit.norm/np.trapz(fit[:,-3]), align="mid", linewidth=0.5, histtype="step", bins=hdata.edges, color="#5ac8fa", label="VS minus")
  # ax[i].hist(fit[:,i], weights=fit[:,-4]*hdata.norm/hfit.norm/np.trapz(fit[:,-4]), align="mid", linewidth=0.5, histtype="step", bins=hdata.edges, color="#ff8388", label="VS plus")
  # ax[i].hist(fit[:,i], weights=fit[:,-5]*hdata.norm/hfit.norm, align="mid", linewidth=0.5, histtype="step", bins=hdata.edges, color="#ff2d55", label="VV D")
  # ax[i].hist(fit[:,i], weights=fit[:,-6]*hdata.norm/hfit.norm, align="mid", linewidth=0.5, histtype="step", bins=hdata.edges, color="#626cf0", label="VV P")
  # ax[i].hist(fit[:,i], weights=fit[:,-7]*hdata.norm/hfit.norm, align="mid", linewidth=0.5, histtype="step", bins=hdata.edges, color="#3bef56", label="VV S")
  # ax[i].hist(fit[:,i], weights=fit[:,-4]*hdata.norm/hfit.norm, align="mid", linewidth=0.5, histtype="step", bins=hdata.edges, color="#626cf0", label="Vector [S,P,D]")
  # ax[i].hist(fit[:,i], weights=fit[:,-3]*hdata.norm/hfit.norm, align="mid", linewidth=0.5, histtype="step", bins=hdata.edges, color="#3bef56", label="Scalar [VS+, VS-, SS]")
  # ax[i].hist(fit[:,i], weights=fit[:,-2]*hdata.norm/hfit.norm, align="mid", linewidth=0.5, histtype="step", bins=hdata.edges, color="#ff2d55", label="Interference")

  ax[i].set_ylim(ymin=0)
  ax[i].set_ylabel("Candidates / {:.2f} {}".format(hdata.bins[1] - hdata.bins[0], units[i]), fontsize=14)

  axpull[i].fill_between(hdata.bins, hpull, 0, facecolor="dodgerblue", alpha=1.0)
  u = ""
  if units[i] !="":
    u = "[{}]".format(units[i])

  axpull[i].set_xlabel(r"{} {}".format(xlabels[i], u), fontsize=12)

  xmin, xmax = axpull[i].get_xlim()

  axpull[i].hlines(y=3.0, xmin=xmin, xmax=xmax, linewidth=1.0, color="lightsteelblue", linestyle="dotted")
  axpull[i].hlines(y=-3.0, xmin=xmin, xmax=xmax, linewidth=1.0, color="lightsteelblue",linestyle="dotted")
  axpull[i].hlines(y=5.0, xmin=xmin, xmax=xmax, linewidth=1.0, color="lightsteelblue", linestyle="dashed")
  axpull[i].hlines(y=-5.0, xmin=xmin, xmax=xmax, linewidth=1.0, color="lightsteelblue", linestyle="dashed")

  # ax[i].set_yscale('log')

  # axpull[i].fill_between([xmin, xmax], 0, 3, facecolor="lightsteelblue", alpha=0.3)
  # axpull[i].fill_between([xmin, xmax], 3, 5, facecolor="lightsteelblue", alpha=0.1)
  # axpull[i].fill_between([xmin, xmax], 0, -3, facecolor="lightsteelblue", alpha=0.3)
  # axpull[i].fill_between([xmin, xmax], -5, -3, facecolor="lightsteelblue", alpha=0.1)

ymax_cos = max(ax[0].get_ylim()[1], ax[1].get_ylim()[1])
ymax_mass = max(ax[3].get_ylim()[1], ax[4].get_ylim()[1])

ax[0].set_ylim(ymax=ymax_cos)
ax[1].set_ylim(ymax=ymax_cos)
ax[3].set_ylim(ymax=ymax_mass)
ax[4].set_ylim(ymax=ymax_mass)


handles, labels = ax[0].get_legend_handles_labels()

ax6.axis("off")
ax6.legend(handles, labels, loc='center')

# fig.tight_layout()

fig.savefig(args["output_fig"])
