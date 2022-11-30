import sys
path_to_add=r'C:\Users\\itaro\OneDrive\Documents\GitHub'
sys.path.insert(0, path_to_add)

from Crystal_structure.tools.general_constants import *
from Crystal_structure.tools.file_handlers import *
from Crystal_structure.tools.rmsd_wrappers import *
from Crystal_structure.tools.tab_data import *

import numpy as np
import os
import re

ESCALTION_INSTRUCTIONS=[['+x', '-x', '+y', '-y', '+z', '-z'],
                        ['+x', '-x', '+y', '-y', '+z', '-z', '-x+y', '-y-x', '-z+x', '-z-x', '-z+y', '-z-y'],
                        ['+x', '-x', '+y', '-y', '+z', '-z', '+x+y', '+x-y','-x+y', '-x-y','+x+z', '+x-z', '-x+z', '-x-z', '+z+y', '-z+y', '+z-y', '-z-y'],
                        ['+x', '-x', '+y', '-y', '+z', '-z', '+x+y', '+x-y','-x+y', '-x-y','+x+z', '+x-z', '-x+z', '-x-z', '+z+y', '-z+y', '+z-y', '-z-y',
                         '+x+y+z', '+x+y-z', '+x-y+z', '+x-y-z', '-x+y+z', '-x+y-z', '-x-y+z', '-x-y-z'],
                        ]

def calc_dist(location_1, location_2): #Utilize?
    distance=sum((np.array(location_1)-np.array(location_2))**2)
    return distance

def get_number_of_molecules_from_cif_labels(atom_labels): # ['C_1', 'C_2', ..]
    possible_number=atom_labels[-1].split('_')[-1]
    if possible_number.isdigit():
        return int(possible_number)
    else:
        return 1

def get_molecule_indices(xyz_df, radii_ext=0.08):
        bonded_atoms_dist_df=connect_the_dots(xyz_df, radii_ext)
        bond_pairs=get_bond_pairs(bonded_atoms_dist_df)
        atom_indices=get_molecule_connections(bond_pairs)
        return atom_indices

def check_for_csd_id_from_cif(cif_filename):
    csd_index=int(cif_filename.split('.')[0])
    from ccdc.search import TextNumericSearch
    searcher=TextNumericSearch()
    searcher.add_ccdc_number(csd_index)
    hits=searcher.search()
    for hit in hits:
        if csd_index==hit.entry.ccdc_number:
            return hit.identifier
    return False   

def diffpy_load_crystal_from_cif(cif_filename):
    from diffpy.structure import loadStructure
    diffpy_structure=loadStructure(cif_filename)
    return diffpy_structure

def diffpy_save_single_xyz_file(diffpy_structure, output_filename):
    diffpy_structure.write(output_filename, format='xyz')

def diffpy_structure_to_xyz_df(diffpy_structure):
    xyz_df=pd.concat([pd.DataFrame(diffpy_structure.element), pd.DataFrame(diffpy_structure.xyz_cartn)], axis=1)
    xyz_df.columns=['element', 'x', 'y', 'z']
    return xyz_df

def diffpy_mutate_single_location(frac_coords, atom_index, mutate_option='+x'): #Change to do x-y-z corralted translations
    orignal_value=frac_coords[atom_index]
    if mutate_option=='+x': # [x, y, z] -> [x+1, y, z]
        new_value=(np.array(orignal_value)+np.array([1.0, 0.0, 0.0])).tolist()
    if mutate_option=='-x': # [x, y, z] -> [x-1, y, z]
        new_value=(np.array(orignal_value)+np.array([-1.0, 0.0, 0.0])).tolist()
    if mutate_option=='+y': # [x, y, z] -> [x, y+1, z]
        new_value=(np.array(orignal_value)+np.array([0.0, 1.0, 0.0])).tolist()
    if mutate_option=='-y': # [x, y, z] -> [x, y-1, z]
        new_value=(np.array(orignal_value)+np.array([0.0, -1.0, 0.0])).tolist()
    if mutate_option=='+z': # [x, y, z] -> [x, y, z+1]
        new_value=(np.array(orignal_value)+np.array([0.0, 0.0, 1.0])).tolist()
    if mutate_option=='-z': # [x, y, z] -> [x, y, z-1]
        new_value=(np.array(orignal_value)+np.array([0.0, 0.0, -1.0])).tolist()
    frac_coords[atom_index]=new_value
    return frac_coords

def diffpy_mutate_all_locations(diffpy_structure, mutate_option='+x'):
    new_diffpy_structure=diffpy_structure.copy()
    frac_coords=new_diffpy_structure.xyz
    for atom_index in range(len(frac_coords)):
        frac_coords=diffpy_mutate_single_location(frac_coords, atom_index, mutate_option)
    new_diffpy_structure.xyz=frac_coords
    return new_diffpy_structure

def diffpy_append_two_structures(diffpy_structure1, diffpy_structure2):
    new_diffpy_structure=diffpy_structure1.copy()
    for atom_index in range(len(diffpy_structure2.xyz)):
        new_diffpy_structure.addNewAtom(atype=diffpy_structure2.element[atom_index],
                                        xyz=diffpy_structure2.xyz[atom_index],
                                        )
    return new_diffpy_structure

def convert_diffpy_structure_to_supercell(diffpy_structure,
                                          mutate_instructions=['+x', '-x', '+y', '-y', '+z', '-z', '+x+y', '+x+z']): 
    new_diffpy_structure=diffpy_structure.copy()
    for mutate_option in mutate_instructions:
        if len(mutate_option)==2:
            mutated_unit_cell=diffpy_mutate_all_locations(diffpy_structure, mutate_option)
        if len(mutate_option)==4:
            mutated_unit_cell=diffpy_mutate_all_locations(diffpy_structure, mutate_option[:2])
            mutated_unit_cell=diffpy_mutate_all_locations(mutated_unit_cell, mutate_option[2:])
        if len(mutate_option)==6:
            mutated_unit_cell=diffpy_mutate_all_locations(diffpy_structure, mutate_option[:2])
            mutated_unit_cell=diffpy_mutate_all_locations(mutated_unit_cell, mutate_option[2:4])
            mutated_unit_cell=diffpy_mutate_all_locations(mutated_unit_cell, mutate_option[4:])
        new_diffpy_structure=diffpy_append_two_structures(new_diffpy_structure, mutated_unit_cell)
    return new_diffpy_structure

class CIFToXYZ():
    def __init__(self):
        pass

    def get_number_of_molecules(self, **kwargs):
        return get_number_of_molecules_from_cif_labels(**kwargs)
##        possible_number_of_molecules=atom_labels[-1].split('_')[-1]
##        if possible_number_of_molecules.isdigit():
##            return int(possible_number_of_molecules)
##        else:
##            return 1

    def record_diffpy_structure_properties(self, diffpy_structure):
        self.number_of_molecules=self.get_number_of_molecules(atom_labels=diffpy_structure.label)
        self.original_size=len(diffpy_structure.label)
        self.num_atoms_in_molecule=self.original_size//self.number_of_molecules 
    
    def get_molecule_indices(self, **kwargs):
        return get_molecule_indices(**kwargs)

    def save_full_xyz(self, diffpy_structure, xyz_filename):
        diffpy_structure.write(xyz_filename, format='xyz')

    def save_individual_xyz(self, diffpy_structure, individual_filename_suffix='_mod.xyz', center_individual_xyz=True, radii_ext=0.08):
        from operator import itemgetter
        self.individual_xyz_filenames=[]
        xyz_df=diffpy_structure_to_xyz_df(diffpy_structure)
        molecule_indices_list=self.get_molecule_indices(xyz_df=xyz_df, radii_ext=radii_ext)
        self.molecule_number=0
        for single_molecule_indices in molecule_indices_list:
            if len(single_molecule_indices)==self.num_atoms_in_molecule:
                coordinates=np.array(itemgetter(*single_molecule_indices)(diffpy_structure.xyz_cartn))
                symbols=np.array(itemgetter(*single_molecule_indices)(diffpy_structure.element))
                if center_individual_xyz:
                    coordinates=rmsd_center_coords(coordinates)
                individual_xyz_filename=str(self.molecule_number)+'_'+individual_filename_suffix
                save_single_xyz_file(symbols, coordinates, output_filename=individual_xyz_filename)
                self.individual_xyz_filenames.append(individual_xyz_filename)
                self.molecule_number+=1
  
    def load_crystal_from_cif(self, cif_filename, mutate_instructions=['+x', '-x', '+y', '-y', '+z', '-z'],
                              record_full_xyz=False, record_individual_xyz=False, center_individual_xyz=True, remove_cif=False):
        from diffpy.structure import loadStructure
        diffpy_structure=loadStructure(cif_filename)
        self.record_diffpy_structure_properties(diffpy_structure)
##        if super_cell:
        new_diffpy_structure=convert_diffpy_structure_to_supercell(diffpy_structure, mutate_instructions)
        clean_name=cif_filename.split('.')[0]
        if record_full_xyz:
            diffpy_save_single_xyz_file(new_diffpy_structure, output_filename=clean_name+'_full'+FileExtensions.XYZ.value)
        if record_individual_xyz:
            self.save_individual_xyz(new_diffpy_structure, individual_filename_suffix=clean_name+FileExtensions.XYZ.value, center_individual_xyz=center_individual_xyz)
            if self.molecule_number==0:
                for new_mutate_instructions in ESCALTION_INSTRUCTIONS:
                    new_diffpy_structure=convert_diffpy_structure_to_supercell(diffpy_structure, new_mutate_instructions)
                    self.save_individual_xyz(new_diffpy_structure, individual_filename_suffix=clean_name+FileExtensions.XYZ.value, center_individual_xyz=center_individual_xyz)
            else:
                print('saved '+str(self.molecule_number)+' individual xyz files')
        if remove_cif:
            os.remove(cif_filename)
            
    def load_crystal_from_csd_id(self, csd_id, mutate_instructions=['+x', '-x', '+y', '-y', '+z', '-z'],
                                 record_full_xyz=False, record_individual_xyz=False, center_individual_xyz=True, remove_cif=False):
        from ccdc.io import CrystalReader, CrystalWriter
        cif_filename=csd_id+FileExtensions.CIF.value
        crystal_reader=CrystalReader('CSD')
        crystal=crystal_reader.crystal(csd_id)
        with CrystalWriter(cif_filename) as crystal_writer: # append=True
            crystal_writer.write(crystal)
        self.load_crystal_from_cif(cif_filename, mutate_instructions=mutate_instructions,
                                   record_full_xyz=record_full_xyz, record_individual_xyz=record_individual_xyz,
                                   center_individual_xyz=center_individual_xyz, remove_cif=remove_cif)

if __name__=='__main__':
    cif_filename='1106253.cif'
    cif_to_xyz=CIFToXYZ()
    cif_to_xyz.load_crystal_from_cif(cif_filename,
##                                     super_cell=True,
                                     mutate_instructions=['+x', '-x', '+y', '-y', '+z', '-z', '+x+y', '+x+z'],
                                     record_full_xyz=False,
                                     record_individual_xyz=True,
                                     center_individual_xyz=True,
                                     remove_cif=False)
