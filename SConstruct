import os
import subprocess
from sys import platform

env = Environment(ENV=os.environ)
env.Tool("compilation_db")
env.CompilationDatabase()

# Set paths to qft, eigen
EIGEN_PATH = os.path.join(os.environ["CONDA_PREFIX"], "include/eigen3")
QFT_PATH = os.path.join(os.environ["CONDA_PREFIX"], "external/qft")
QFT_INC = QFT_PATH + "/include"
QFT_LIB = QFT_PATH + "/lib"
BOOST_INC = os.path.join(os.environ["CONDA_PREFIX"], "include/boost")
BOOST_LIB = os.path.join(os.environ["CONDA_PREFIX"], "lib")


if platform == "darwin":
    env.Append(
        CXXFLAGS=[
            "-O3",
            "-Wall",
            "-Wextra",
            "-pedantic",
            "-std=c++20",
            "-Xpreprocessor",
            "-fopenmp",
        ]
    )
    env.Append(LIBS=["Minuit2", "MathMore", "qft++", "omp", "boost_program_options"])
    env.Append(LIBPATH=[QFT_LIB, BOOST_LIB])
else:
    env.Append(
        CXXFLAGS=[
            "-O3",
            "-Wall",
            "-Wextra",
            "-pedantic",
            "-std=c++20",
            "-fopenmp",
        ]
    )
    env.Append(LIBS=["Minuit2", "MathMore", "qft++", "boost_program_options"])

    env.Append(LIBPATH=[QFT_LIB, BOOST_LIB])
    env.Append(LINKFLAGS=["-fopenmp"])

# Set environment
env.AppendUnique(CCFLAGS=("-isystem", EIGEN_PATH))
env.AppendUnique(CCFLAGS=("-isystem", QFT_INC))
env.AppendUnique(CCFLAGS=("-isystem", BOOST_INC))
env.Append(CPPPATH=["#", "#framework"])
env.ParseConfig("root-config --libs --cflags --ldflags")

# Export environment
Export("env")

# Run SConscript files
env.SConscript("src/SConscript")
env.SConscript("test/SConscript")
