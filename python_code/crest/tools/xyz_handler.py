import sys
path_to_add=r'C:\Users\\itaro\OneDrive\Documents\GitHub'
sys.path.insert(0, path_to_add)

from Crystal_structure.tools.general_constants import *
from Crystal_structure.tools.file_handlers import *
from Crystal_structure.tools.rmsd_wrappers import *
from Crystal_structure.tools.tab_data import *

import numpy as np
import os

def get_tab_data(xyz_filename):
    molecule_tab=TabDataAnalyzer(set_tab_parser(), xyz_filename)
    data_dict, plot_filename=molecule_tab.get_tab_data()
    return molecule_tab, data_dict, plot_filename

class XYZFilesCompare():
    def __init__(self):
        pass
    
    def transform_coords_for_RMSD_calc(self, coordinates_1, coordinates_2):
        coordinates_1, coordinates_2=rmsd_translate_molecules(coordinates_1, coordinates_2)
        coordinates_1, coordinates_2=rmsd_rotate_molecules(coordinates_1, coordinates_2)
        return coordinates_1

    def calculate_molecule_RMSD(self, coordinates_1, coordinates_2): #save_new_coords=False, symbols_1=None, export_filename=None, save_new_coords_plot=False
        RMSD=rmsd_get_RMSD_score(coordinates_1,coordinates_2)
        return round(RMSD, RMSDConstants.ROUND_DIGITS.value)

    def save_tab_data(self, xyz_filename, attr_name_prefix):
        molecule_tab, data_dict, plot_filename=get_tab_data(xyz_filename)
        setattr(self, attr_name_prefix+'_tab', molecule_tab)
        setattr(self, attr_name_prefix+'_data_dict', data_dict)
        setattr(self, attr_name_prefix+'_plot_name', plot_filename)        

    def set_first_molecule(self, xyz_filename, get_plot=False, center_molecule=True):
        symbols, coordinates=load_single_xyz_file(xyz_filename)
        self.first_molecule_symbols=symbols
        if center_molecule:
            coordinates=rmsd_center_coords(coordinates)
            save_single_xyz_file(symbols, coordinates, output_filename='centered_'+xyz_filename)
        self.first_molecule_coordinates=coordinates
        if get_plot:
            self.save_tab_data(xyz_filename, 'first_molecule')

    def align_according_to_first_molecule(self, symbols, coordinates):
        try:
            aligned_symbols, aligned_coordinates=rmsd_renumbering_coordinates(symbols, self.first_molecule_symbols, coordinates, self.first_molecule_coordinates)
        except:
            aligned_symbols, aligned_coordinates=symbols, coordinates
        return aligned_symbols, aligned_coordinates

    def get_rmsd_results(self, coordinates_1, coordinates_2):
        try:
            rmsd_score=self.calculate_molecule_RMSD(coordinates_1, coordinates_2)
            return rmsd_score
        except:
            return None

    def process_rmsd_results(self, rmsd_values):
        data=RMSDData()
        for rmsd_index, rmsd_value in enumerate(rmsd_values):
            setattr(data, 'rmsd_'+str(rmsd_index), rmsd_value)
        setattr(data, 'min_rmsd', round(min(rmsd_values), RMSDConstants.ROUND_DIGITS.value))
        setattr(data, 'max_rmsd', round(max(rmsd_values), RMSDConstants.ROUND_DIGITS.value))
        setattr(data, 'avg_rmsd', round(np.mean(rmsd_values), RMSDConstants.ROUND_DIGITS.value))
        return data        
        
    def compare_molecules(self, xyz_filenames, get_plot=False):
        rmsd_values=[]
        for xyz_filename in xyz_filenames:
            compare_symbols, compare_coordinates=load_single_xyz_file(xyz_filename)
            aligned_symbols, aligned_coordinates=self.align_according_to_first_molecule(compare_symbols, compare_coordinates)         
            adjusted_coordinates=self.transform_coords_for_RMSD_calc(aligned_coordinates, self.first_molecule_coordinates)
            adjusted_xyz_filename='adjusted_'+xyz_filename
            save_single_xyz_file(aligned_symbols, adjusted_coordinates, output_filename=adjusted_xyz_filename)
            possible_score=self.get_rmsd_results(adjusted_coordinates, self.first_molecule_coordinates)
            if possible_score:
                rmsd_values.append(possible_score)
                if get_plot:
                    self.save_tab_data(adjusted_xyz_filename, str(len(rmsd_values))+'_rmsd_molecule')
        if len(rmsd_values)>0:
            return self.process_rmsd_results(rmsd_values)
        else:
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
    print('All Good!')
