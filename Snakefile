import os

"""
  Define all paths from current directory
"""
code_path = ""
tuples_path = os.path.join(code_path, "tuples/")
bin_path = os.path.join(code_path, "bin/")
python_path = os.path.join(code_path, "python/")

norm_path = os.path.join(tuples_path, "normalisation/")
norm_Bs_path = os.path.join(norm_path, "Bs/")
norm_Du_path = os.path.join(norm_path, "Du/")

data_path = os.path.join(tuples_path, "signal/")
data_Bs_path = os.path.join(tuples_path, "signal/Bs/")
data_Du_path = os.path.join(tuples_path, "signal/Du/")


def create_paths(path):
    return f"[[ -d {path} ]] || mkdir -p {path}"


os.system(create_paths(tuples_path))
os.system(create_paths(norm_Bs_path))
os.system(create_paths(norm_Du_path))
os.system(create_paths(data_Bs_path))
os.system(create_paths(data_Du_path))

include: "snakefiles/functions.snake"

"""
  Simulation config
"""
nsig_Bs = 3_000  # number of signal events for Bs0
nsig_Du = 3_000  # number of signal events for D0
nint_Bs = 1_000_000  # number of integration events for Bs0
nint_Du = 1_000_000  # number of integration events for D0

available_cpu = workflow.cores

events_per_job_for_norm_Bs = get_events_per_job(nint_Bs, available_cpu)
events_per_job_for_norm_Du = get_events_per_job(nint_Du, available_cpu)

norm_events_per_job = {
    'Bs': events_per_job_for_norm_Bs, 
    'Du': events_per_job_for_norm_Du, 
}

events_per_job_for_sig_Bs = get_events_per_job(nsig_Bs, available_cpu)
events_per_job_for_sig_Du = get_events_per_job(nsig_Du, available_cpu)

sig_events_per_job = {
    'Bs': events_per_job_for_sig_Bs, 
    'Du': events_per_job_for_sig_Du, 
}

"""
  Wildcard constraints
"""
wildcard_constraints:
    meson="(Bs|Du)"


"""
  Include the snakefiles
"""
include: "snakefiles/integration_events.snake"
include: "snakefiles/generate_signal.snake"
include: "snakefiles/fit_signal.snake"
include: "snakefiles/plot_signal.snake"
include: "snakefiles/scan_signal.snake"
