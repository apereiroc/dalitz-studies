CURRENT_PATH=`pwd`

mkdir -p $CONDA_PREFIX/external &&
git clone https://github.com/jdalseno/qft.git $CONDA_PREFIX/external/qft &&
cd $CONDA_PREFIX/external/qft  &&
scons -j7 &&
ln -sf $CONDA_PREFIX/external/qft/lib/libqft++.* $CONDA_PREFIX/lib &&

cd $CURRENT_PATH
