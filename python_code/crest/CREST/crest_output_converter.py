import sys
path_to_add=r'C:\Users\edens\Documents\GitHub'
sys.path.insert(0, path_to_add)

import os
import numpy as np
import pandas as pd

import zipfile

from XTB.xtb_operation_constants import *
from XTB.xtb_output_converter import *
from CREST.crest_operation_constants import *

##def get_file_lines(filename, encoding=None):
##    with open(filename, 'r', encoding=encoding) as f:
##        lines=f.readlines()
##    return lines

def get_xyz_from_crest_zipfile(zip_filename):
    with zipfile.ZipFile(zip_filename, 'r', zipfile.ZIP_DEFLATED) as conformer_zip:
        small_zipfile_path=conformer_zip.extract('0_'+zip_filename) # this is only the lowest energy
        with zipfile.ZipFile(small_zipfile_path, 'r', zipfile.ZIP_DEFLATED) as xtb_zip:
            conformer_xyz_filepath=xtb_zip.extract(XTBOutputFilenames.OPT_XYZ.value)
        new_filename=conformer_xyz_filepath.split('.')[0]+zip_filename.split('_')[0]+'.xyz'
        os.rename(conformer_xyz_filepath, new_filename)
    os.remove(small_zipfile_path)
    return new_filename

class CrestHandler():
    def __init__(self, csd_id, xyz_filename):
        self.csd_id=csd_id
        self.xyz_filename=xyz_filename

    def init_crest_calculation(self):
        run_crest_calculation(self.xyz_filename)
        self.file_lines=get_file_lines(CrestConstants.CREST_OUTPUT_NAME.value)
        self.number_of_conformers=len(self.file_lines)//self.number_of_single_molecule_lines

    def run_xtb_calculations_from_crest(self):
        dict_list=[]
        zip_file_namelist=[]
        for conformer_index in range(self.number_of_conformers):
            single_xyz_filename=str(conformer_index)+CrestConstants.SINGLE_CONFORMER_XYZ.value
            extract_single_xyz_file_from_ensamble(ensamble_xyz_file=CrestConstants.CREST_OUTPUT_NAME.value,
                                                  file_index=conformer_index,
                                                  output_filename=single_xyz_filename,
                                                  )
            try:
                run_xtb_calculation(single_xyz_filename, detailed_input='xtb_di.inp')
                xtb_handler=XTBHandler(XTBOutputFilenames.PROPERTIES.value, encode_data=True, data_container=XTBDataContainer(), csd_id=self.csd_id)
                xtb_handler.encode_data()
                xtb_handler.encode_molecule_volume(single_xyz_filename)
                dict_list.append(xtb_handler.__dict__)
                zip_filename=str(conformer_index)+'_'+self.csd_id+'_crest.zip'
                save_to_new_zipfile(zip_filename, filenames_to_zip=[single_xyz_filename]+[member.value for name, member in XTBOutputFilenames.__members__.items()])
                zip_file_namelist.append(zip_file_name)
            except:
                pass
        results_df=pd.DataFrame(dict_list)
        csv_file_name=self.csd_id+CrestConstants.SINGLE_CONFORMER_XYZ.value
        results_df.to_csv(csv_file_name)
        full_zip_file_name=self.csd_id+'_crest.zip'
        with zipfile.ZipFile(full_zip_file_name, 'w', zipfile.ZIP_DEFLATED) as full_zip_file:
            for filename in zip_file_namelist:
                full_zip_file.write(filename)
                os.system(LinuxCommands.REMOVEL_COMMAND.value+filename)
            for name, member in CrestOutputFilenames.__members__.items():
                full_zip_file.write(member.value)
        return full_zip_file_name, csv_file_name


if __name__=='__main__':
    print('fill in data!')
    
