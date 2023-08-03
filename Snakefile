import yaml
import os

# def get_user():
#     current_dir = os.getcwd()
#     if "asier.pereiro" in current_dir or "apereiro" in current_dir:
#         return "asier"
#     else:
#         raise NotImplementedError("User not defined")
#
# user = get_user()

user = "asier"

with open("config.yaml", "r") as file:
  config = yaml.load(file, Loader=yaml.FullLoader)[user]

"""
  Define all paths from code and tuples root paths
"""

if config['code_path'][-1] != '/':
    config['code_path'] += '/'

if config['tuples_path'][-1] != '/':
    config['tuples_path'] += '/'

fitter_path = config["code_path"]
tuples_path = config["tuples_path"]
phsp_mc_Bs_path = config["tuples_path"] + "normalisation/Bs/"
phsp_mc_Du_path = config["tuples_path"] + "normalisation/Du/"
data_Bs_path = config["tuples_path"] + "signal/Bs/"
data_Du_path = config["tuples_path"] + "signal/Du/"

def create_paths(path):
    return f"[[ -d {path} ]] || mkdir -p {path}"

os.system(create_paths(tuples_path))
os.system(create_paths(phsp_mc_Bs_path))
os.system(create_paths(phsp_mc_Du_path))
os.system(create_paths(data_Bs_path))
os.system(create_paths(data_Du_path))

"""
  Include the snakefiles
"""
include: "snakefiles/integration_events.snake"
