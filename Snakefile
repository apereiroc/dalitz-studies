import os

"""
  Define all paths from current directory
"""
code_path = os.getcwd() + '/'
tuples_path = "tuples/"

phsp_mc_Bs_path = tuples_path + "normalisation/Bs/"
phsp_mc_Du_path = tuples_path + "normalisation/Du/"
data_Bs_path = tuples_path + "signal/Bs/"
data_Du_path = tuples_path + "signal/Du/"


def create_paths(path):
    return f"[[ -d {path} ]] || mkdir -p {path}"


os.system(create_paths(tuples_path))
os.system(create_paths(phsp_mc_Bs_path))
os.system(create_paths(phsp_mc_Du_path))
os.system(create_paths(data_Bs_path))
os.system(create_paths(data_Du_path))

"""
  Simulation config
"""
nsig_Bs = 3_000  # number of signal events for Bs0
nsig_Du = 3_000  # number of signal events for D0
nint_Bs = 1_000_000  # number of integration events for Bs0
nint_Du = 1_000_000  # number of integration events for D0

"""
  Include the snakefiles
"""
include: "snakefiles/functions.snake"
include: "snakefiles/integration_events.snake"
include: "snakefiles/generate_signal.snake"
include: "snakefiles/fit_signal.snake"
include: "snakefiles/plot_signal.snake"
include: "snakefiles/scan_signal.snake"
