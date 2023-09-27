Import('*')

import glob
import os

for i in glob.glob('*.cpp'):
    bin_name, bin_ext = os.path.splitext(i)
    env.Program(bin_name, [i])
