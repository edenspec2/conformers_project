#!/gpfs0/gaus/users/itamarwa/Python-3.7.10/python

import os
from enum import Enum
import numpy as np
import pandas as pd
import multiprocessing as mp
import zipfile
path=r'C:\Users\edens\Documents\GitHub\conformers_project\python_code\crest'
os.chdir(path)

molecules_dir=[molecule_dir for molecule_dir in os.listdir() if os.path.isdir(molecule_dir)] 

for molecule_dir in molecules_dir:
    molecule_name=molecule_dir
    path=os.path.abs(molecule_dir)
    os.mkdir(path)
    os.chdir(path)
    xyz_filename=molecule_name+'.xyz'
    os.system('cp '+Constant.HOME_DIR.value+'/fixed_locations.inp '+Constant.HOME_DIR.value+'/'+molecule_name+'/fixed_locations.inp')
    os.system('cp '+Constant.HOME_DIR.value+'/'+smile_filename+' '+Constant.HOME_DIR.value+'/'+molecule_name+'/'+smile_filename)
    os.system(CrestConstants.OBABEL.value+smile_filename+CrestConstants.OBABEL_XYZ_SETTINGS_1.value+xyz_filename+CrestConstants.OBABEL_XYZ_SETTINGS_2.value)
    os.system('echo >> '+xyz_filename)
    os.system('echo >> '+xyz_filename)
    os.system('echo "\$write" >> '+xyz_filename)
    os.system('echo "   output file=properties.out" >> '+xyz_filename)
    os.chdir(Constant.HOME_DIR.value)
    os.system(CrestConstants.REMOVEL_COMMAND.value+smile_filename)
    xyz_filenames.append(xyz_filename)

