import sys
path_to_add=r'C:\Users\\itaro\OneDrive\Documents\GitHub'
sys.path.insert(0, path_to_add)

import numpy as np
import os

from XTB.xtb_operation_constants import *

class XTBDataContainer():
    pass
             
class XTBHandler():
    def __init__(self, xtb_filename, encode_data=False, data_container=None, csd_id=None):
        self.xtb_filename=xtb_filename
        if encode_data:
            self.data_container=data_container
            if csd_id:
                setattr(self.data_container, 'CSD_ID', csd_id)

    def get_file_striped_lines(self, filename, encoding=None):
        f=open(filename, 'r', encoding=encoding)
        lines=f.readlines()
        f.close()
        strip_lines=[line.strip().rstrip('\n') for line in lines]
        return strip_lines

    def extract_numeric_values(self, list_var):
        return list(map(float, filter(is_number, list_var)))

    def extract_negative_values(self, list_var):
        num_list=self.extract_numeric_values(list_var)
        return list(filter(lambda x: x<0, num_list))

    def save_features(self, features_dict):
        for feature_name, feature_value in features_dict.items():
            setattr(self.data_container, feature_name, feature_value)

    def save_energy_levels_features(self, file_lines):
        for index, line in enumerate(file_lines):
            if line.endswith('(HOMO)'):
                energies=[self.extract_negative_values(file_lines[index-2].split())[-1],
                          self.extract_negative_values(file_lines[index-1].split())[-1],
                          self.extract_negative_values(file_lines[index].split())[-1],
                          self.extract_negative_values(file_lines[index+1].split())[-1],
                          self.extract_numeric_values(file_lines[index+2].split())[-1],
                          self.extract_numeric_values(file_lines[index+3].split())[-1],
                          ]
                self.save_features(dict(zip(XTBFeaturesNames.ENERGY_LEVELS,value, energies)))

    def save_dispersion_polarizability_features(self, file_lines):
        start_index=[index for index, line in enumerate(file_lines) if line.startswith('Mol. C6AA')][0]
        values=[self.extract_numeric_values(file_lines[start_index].split())[-1],
                self.extract_numeric_values(file_lines[start_index+1].split())[-1],
                self.extract_numeric_values(file_lines[start_index+2].split())[-1],
                ]
        self.save_features(dict(zip(XTBFeaturesNames.DISPERSION_POLARIZABILITY.value, values)))

    def save_dipole_features(self, file_lines):
        start_index=[index for index, line in enumerate(file_lines) if line.startswith('molecular dipole:')][0]
        dipole_q_line=self.extract_numeric_values(file_lines[start_index+2].split())
##        dipole_q_line.append(0)
        dipole_q_line.append(((dipole_q_line[0]**2)+(dipole_q_line[1]**2)+(dipole_q_line[2]**2)**0.5))
        self.save_features(dict(zip(XTBFeaturesNames.DIPOLE_Q.value, dipole_q_line)))
        dipole_line=self.extract_numeric_values(file_lines[start_index+3].split())
        self.save_features(dict(zip(XTBFeaturesNames.DIPOLE.value, dipole_line)))

    def save_quadrapole_features(self, file_lines):
        start_index=[index for index, line in enumerate(file_lines) if line.startswith('molecular quadrupole')][0]
        quad_q_line=self.extract_numeric_values(file_lines[start_index+2].split())
        self.save_features(dict(zip(XTBFeaturesNames.QUAD_Q.value, quad_q_line)))
        quad_q_dip_line=self.extract_numeric_values(file_lines[start_index+3].split())
        self.save_features(dict(zip(XTBFeaturesNames.QUAD_Q_DIP.value, quad_q_dip_line)))
        quad_full_line=self.extract_numeric_values(file_lines[start_index+4].split())
        self.save_features(dict(zip(XTBFeaturesNames.QUAD_FULL.value, quad_full_line)))

    def save_interia_features(self, file_lines):
        start_index=[index for index, line in enumerate(file_lines) if line.startswith('moments of inertia')][0]
        inertia_line=self.extract_numeric_values(file_lines[start_index].split())
        self.save_features(dict(zip(XTBFeaturesNames.INERTIA_MOMENTS.value, inertia_line)))

    def save_single_bond_details(self, bond_line):
        split_result=bond_line.split()
        base_feature_name=split_result[1]+'-'+split_result[3]+' '
        final_features_names=[base_feature_name+sufix for sufix in XTBFeaturesNames.BONDS_STATS.value]
        self.save_features(dict(zip(final_features_names, split_result[5:])))

    def save_bonds_stats_features(self, file_lines):
        for index, line in enumerate(file_lines):
            split_result=line.split()
            if len(split_result)>=4:
                if split_result[3].startswith('av.'):
                    start_index=index
        number_of_bonds=int(file_lines[start_index-2].split()[1])
        for bond_number in range(1, number_of_bonds+1):
            self.save_single_bond_details(file_lines[start_index+bond_number])

    def save_thermo_features(self, file_lines):
        start_index=[index for index, line in enumerate(file_lines) if line.startswith('298.15  VIB')][0]
        for type_index, line_type in enumerate(XTBFeaturesNames.THERMO_TYPES.value):
            final_features_names=[line_type+'_'+sufix for sufix in XTBFeaturesNames.THERMO_NAMES.value]
            if line_type=='TOT':
                self.save_features(dict(zip(final_features_names, self.extract_numeric_values(file_lines[start_index+type_index].split()))))
            elif line_type=='VIB':
                self.save_features(dict(zip(final_features_names, self.extract_numeric_values(file_lines[start_index+type_index].split())[2:])))
            else:
                self.save_features(dict(zip(final_features_names, self.extract_numeric_values(file_lines[start_index+type_index].split())[1:])))

    def save_zero_point_energy_features(self, file_lines):
        start_index=[index for index, line in enumerate(file_lines) if line.startswith(':: zero point energy')][0]
        zero_point_line=self.extract_numeric_values(file_lines[start_index].split())
        self.save_features(dict(zip(XTBFeaturesNames.ZERO_ENERGY.value, zero_point_line)))
    
    def encode_data(self):
        file_lines=get_file_striped_lines(self.xtb_filename, encoding='utf-8')
        self.save_energy_levels_features(file_lines)
        self.save_dispersion_polarizability_features(file_lines)
        self.save_dipole_features(file_lines)
        self.save_quadrapole_features(file_lines)
        self.save_interia_features(file_lines)
        self.save_bonds_stats_features(file_lines)
        self.save_thermo_features(file_lines)
        self.save_zero_point_energy_features(file_lines)

    def encode_molecule_volume(self, xyz_filename):
        file_lines=get_file_striped_lines(xyz_filename, encoding='utf-8')
        xyz_table=[self.extract_numeric_values(line.split()) for line in file_lines[2:]]
        xyz_array=np.array(xyz_table, dtype=float)
        molecule_volume=(xyz_array.max(axis=0)-xyz_array.min(axis=0)).prod()
        self.save_features({'molecule_volume': molecule_volume})
        
if __name__=='__main__':
    my_data_container=XTBDataContainer()
    xtb_filename=XTBOutputFilenames.PROPERTIES.value
    xtb_handler=XTBHandler(xtb_filename, encode_data=True, data_container=my_data_container)
    xtb_handler.encode_data()
    xtb_handler.encode_molecule_volume(XTBOutputFilenames.OPT_XYZ.value)
    print(xtb_handler.data_container.__dict__)
