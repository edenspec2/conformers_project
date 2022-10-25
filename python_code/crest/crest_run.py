#!/gpfs0/gaus/users/itamarwa/Python-3.7.10/python

import os
from enum import Enum
import numpy as np
import pandas as pd
import multiprocessing as mp
import zipfile
path=r'C:\Users\edens\Documents\GitHub\conformers_project\python_code\crest'
os.chdir(path)
from XTB.xtb_output_converter import *
from CREST.crest_operation_constants import *
from CREST.crest_output_converter import *

class Constant(Enum):
    # HOME_DIR=r'/gpfs0/gaus/users/edenspec/crest_runs'
    HOME_DIR=r'C:\Users\edens\Documents\GitHub\conformers_project\python_code\crest'

def process_xyz_file(xyz_filename):
    molecule_name=xyz_filename.split('.')[0]
    path=os.path.abspath(molecule_name)
    os.chdir(path)
    crest_handler=CrestHandler(molecule_name, xyz_filename)
    crest_handler.init_crest_calculation()
    full_zip_file_name, csv_file_name=crest_handler.run_xtb_calculations()
    os.chdir(Constant.HOME_DIR.value)
    return path, full_zip_file_name, csv_file_name

os.chdir(Constant.HOME_DIR.value)

csv_file_name='smiles_to_run.csv'
smiles_df=pd.read_csv(csv_file_name)

smiles_array=smiles_df.smiles.to_numpy()
names_array=smiles_df.name.to_numpy()

for smile_index, smile in enumerate(smiles_array):
    smile_filename=names_array[smile_index]+FileExtensions.SMI.value
    try:
        with open(smile_filename,'a') as outfile:
            outfile.write(smile)
    except:
        pass   

smi_filenames=[filename for filename in os.listdir() if filename.endswith(FileExtensions.SMI.value)]
xyz_filenames=[]

for smile_filename in smi_filenames:
    molecule_name=smile_filename.split('.')[0]
    path=os.path.join(Constant.HOME_DIR.value+'/', molecule_name)
    os.mkdir(path)
    os.chdir(path)
    xyz_filename=molecule_name+FileExtensions.XYZ.value
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

pool=mp.Pool(mp.cpu_count())
results=pool.map(process_xyz_file, xyz_filenames)
pool.close() 
seprated_data=list(map(list, zip(*results)))
path_list=seprated_data[0]
zip_file_name_list=seprated_data[1]
csv_file_name_list=seprated_data[2]

##    for xyz_filename in xyz_filenames:
##        os.system(CrestConstants.REMOVEL_COMMAND.value+xyz_filename) 

with zipfile.ZipFile('all_files_zip.zip', 'w', zipfile.ZIP_DEFLATED) as batch_zip:
    for path_index, path in enumerate(path_list):
        os.chdir(path)
        batch_zip.write(zip_file_name_list[path_index])
        batch_zip.write(csv_file_name_list[path_index])
        os.chdir(Constant.HOME_DIR.value)
        os.rmdir(path)
