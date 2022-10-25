# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 14:45:24 2022

@author: edens
"""
from openbabel import pybel
from openbabel import openbabel as ob
import help_functions
import os
import numpy as np
import pandas as pd
from biopandas.mol2 import PandasMol2
from modified_project import *

class Constant(Enum):
    # HOME_DIR=r'/gpfs0/gaus/users/edenspec/crest_runs'
    HOME_DIR=r'/gpfs0/gaus/users/edenspec/molecules'
    
if __name__ == '__main__':
    
    xyz_files=[filename for filename in os.listdir() if 'xyz' in filename]
    os.chdir(Constant.HOME_DIR.value)
    for xyz_file in xyz_files:
        # molecule_name=xyz_file.split('.')[0]
        # path=os.path.join(Constant.HOME_DIR.value+'/', molecule_name)
        # os.mkdir(path)
        # os.chdir(path)
        xyz_filename=xyz_file
        os.system('obabel {} -O {} --conformer --nconf 30--writeconformers'.format(*xyz_filename,('obconformers_'+xyz_filename)))
        # os.system('cp '+Constant.HOME_DIR.value+'/fixed_locations.inp '+Constant.HOME_DIR.value+'/'+molecule_name+'/fixed_locations.inp')
        # os.system('cp '+Constant.HOME_DIR.value+'/'+smile_filename+' '+Constant.HOME_DIR.value+'/'+molecule_name+'/'+smile_filename)
        # os.system(CrestConstants.OBABEL.value+smile_filename+CrestConstants.OBABEL_XYZ_SETTINGS_1.value+xyz_filename+CrestConstants.OBABEL_XYZ_SETTINGS_2.value)
        # os.system('echo >> '+xyz_filename)
        # os.system('echo >> '+xyz_filename)
        # os.system('echo "\$write" >> '+xyz_filename)
        # os.system('echo "   output file=properties.out" >> '+xyz_filename)
        # os.chdir(Constant.HOME_DIR.value)
        # os.system(CrestConstants.REMOVEL_COMMAND.value+smile_filename)
        # xyz_filenames.append(xyz_filename)