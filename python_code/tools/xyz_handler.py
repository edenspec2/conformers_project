import sys
path_to_add=r'C:\Users\\itaro\OneDrive\Documents\GitHub'
sys.path.insert(0, path_to_add)

from Crystal_structure.tools.general_constants import *
from Crystal_structure.tools.file_handlers import *
from Crystal_structure.tools.rmsd_wrappers import *
from Crystal_structure.tools.tab_data import *

import numpy as np
import os

class XYZFilesCompare():
    def __init__(self):
        pass
    
    def rmsd_align_xyz(self, **kwargs):
        rotated_coordinates_1, *_=rmsd_align_xyz(**kwargs)
        return rotated_coordinates_1

    def rmsd_get_RMSD_score(self, **kwargs): #save_new_coords=False, symbols_1=None, export_filename=None, save_new_coords_plot=False
        return rmsd_get_RMSD_score(**kwargs)

    def save_tab_data(self, xyz_filename, attr_name_prefix):
        molecule_tab, data_dict, plot_filename=get_tab_data(xyz_filename)
        setattr(self, attr_name_prefix+'_tab', molecule_tab)
        setattr(self, attr_name_prefix+'_data_dict', data_dict)
        setattr(self, attr_name_prefix+'_plot_name', plot_filename)        

    def set_first_molecule(self, xyz_filename, get_tab_data=False, center_molecule=True):
        symbols, coordinates=load_single_xyz_file(xyz_filename)
        if center_molecule:
            coordinates=rmsd_center_coords(coordinates)
            save_single_xyz_file(symbols, coordinates, output_filename='centered_'+xyz_filename)
            self.first_molecule_df=get_xyz_df_from_file('centered_'+xyz_filename)
        else:
            self.first_molecule_df=get_xyz_df_from_file(xyz_filename)
        self.first_molecule_symbols=symbols
        self.first_molecule_coordinates=coordinates
        self.first_molecule_bond_pairs=get_bond_pairs_from_xyz_df(self.first_molecule_df)
        if get_tab_data:
            self.save_tab_data(xyz_filename, 'first_molecule')

    def align_according_to_first_molecule(self, xyz_df):
        alignment_mapping=get_alignment_mapping_from_xyz_df(xyz_df, self.first_molecule_bond_pairs)
        aligned_xyz_df=align_xyz_df_mapping(xyz_df, alignment_mapping)
        aligned_symbols=aligned_xyz_df['element'].to_numpy()
        aligned_coordinates=aligned_xyz_df[['x', 'y', 'z']].to_numpy()
        return aligned_symbols, aligned_coordinates

    def get_rmsd_results(self, coordinates_1, coordinates_2):
        rmsd_score=self.rmsd_get_RMSD_score(coordinates_1=coordinates_1, coordinates_2=coordinates_2)
        return rmsd_score

    def process_rmsd_results(self, rmsd_values):
        data=RMSDData()
        for rmsd_index, rmsd_value in enumerate(rmsd_values):
            setattr(data, 'rmsd_'+str(rmsd_index), rmsd_value)
        setattr(data, 'min_rmsd', round(min(rmsd_values), RMSDConstants.ROUND_DIGITS.value))
        setattr(data, 'max_rmsd', round(max(rmsd_values), RMSDConstants.ROUND_DIGITS.value))
        setattr(data, 'avg_rmsd', round(np.mean(rmsd_values), RMSDConstants.ROUND_DIGITS.value))
        return data        
        
    def compare_molecules(self, xyz_filenames, get_tab_data=False):
        rmsd_values=[]
        for xyz_filename in xyz_filenames:
            xyz_df=get_xyz_df_from_file(xyz_filename)
            aligned_symbols, aligned_coordinates=self.align_according_to_first_molecule(xyz_df)
            adjusted_coordinates=self.rmsd_align_xyz(coordinates_1=aligned_coordinates, coordinates_2=self.first_molecule_coordinates)
            rmsd_score=self.get_rmsd_results(adjusted_coordinates, self.first_molecule_coordinates)
            rmsd_values.append(rmsd_score)
            if get_tab_data:
                aligned_xyz_filename='aligned_'+xyz_filename
                save_single_xyz_file(aligned_symbols, aligned_coordinates, output_filename=aligned_xyz_filename)
                self.save_tab_data(aligned_xyz_filename, str(len(rmsd_values))+'_rmsd_molecule')
                os.remove(aligned_xyz_filename)
        if len(rmsd_values)>0:           
            return self.process_rmsd_results(rmsd_values)
        return False

    def collect_molecules_attr(self, attr_name):
        attr_list=list(self.__dict__.keys())
        found_attr_names=[possible_attr for possible_attr in attr_list if attr_name in possible_attr]
        return found_attr_names

    def collect_molecules_data_dict_values(self, data_dict_key_name):
        found_data_dict_values=[]
        data_dict_attr_names=self.collect_molecules_attr('data_dict')
        for data_dict_attr in data_dict_attr_names:
            data_dict=getattr(self, data_dict_attr)
            if data_dict_key_name in data_dict.keys():
                found_data_dict_values.append(data_dict.get(data_dict_key_name))
        return found_data_dict_values

if __name__=='__main__':
    home_path=r'C:\Users\\itaro\OneDrive\Documents\GitHub'
##        os.chdir(home_path)
    sys.path.insert(1, path_to_add)
    xyz_filename_1='0_1106253.xyz'
    xyz_filename_2='1_1106253.xyz'
    xyz_filename_3='2_1106253.xyz'
    xyz_file_comparer=XYZFilesCompare()
    xyz_file_comparer.set_first_molecule(xyz_filename_1, get_tab_data=False, center_molecule=True)
    rmsd_results=xyz_file_comparer.compare_molecules([xyz_filename_2, xyz_filename_3], get_tab_data=False)
    print(rmsd_results.__dict__)
