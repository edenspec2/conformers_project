from openbabel import openbabel as ob
from openbabel import pybel

import pandas as pd
from biopandas.mol2 import PandasMol2

# import psi4
from modified_project_conformers import *

os.chdir(r'C:\Users\edens\Documents\GitHub\conformers_project\python_code')

import tab_data

"""
pybel ob diffrences

## pbmol -- pybel molecule // obmol--openbabel molecule
##  pbmol = pybel.Molecule(mol) # to go from OBMol to pybel
##  obmol = ob.OBMol(pbmol.OBMol)  # to go from pybel to open babel  

molecule.obmol.GetFormula()
Out[6]: 'C132H112N12P4S4'

molecule.pbmol.formula
Out[7]: 'C132H112N12P4S4'

  """


# path_to_add=r'C:\Users\edens\Documents\GitHub\conformers_project\python_code\molecules'
# os.chdir(path_to_add)




def split_molecule_file(file_name):  ##for crest output
    file_type = file_name.split('.')[-1]
    molecules = [pbmol.write(file_type) for pbmol in pybel.readfile(file_type, file_name)]
    return molecules


def calc_pbmol_energy(pbmol):
    obmol = ob.OBMol(pbmol.OBMol)
    pff = ob.OBForceField_FindType("mmff94")
    # pff.SetLogLevel(ob.OBFF_LOGLVL_HIGH)
    pff.SetLogToStdErr()
    pff.Setup(obmol)
    # pff.SetLogLevel(ob.OBFF_LOGLVL_NONE)
    return pff.Energy()

def get_pbmol_energy_list(pbmol_list):
    energy_list=pd.DataFrame([calc_pbmol_energy(pbmol) for pbmol in pbmol_list], columns=['energy'])
    return energy_list


def calc_pbmol_charge(pbmol):
    obmol = ob.OBMol(pbmol.OBMol)
    ob_charge_model = ob.OBChargeModel.FindType("eem2015bn")
    ob_charge_model.ComputeCharges(obmol)
    return ob_charge_model.GetPartialCharges()
def get_pbmol_charge_list(pbmol_list:list)->list:
    charge_list=[calc_pbmol_charge(pbmol) for pbmol in pbmol_list]
    return charge_list



#TODO: change so it wont need to write a temp file.
def pbmol_to_mol2_df(pbmol):
    pbmol.write('mol2', 'temp.mol2')
    pmol = PandasMol2()
    pmol.read_mol2(os.path.abspath('temp.mol2'))
    os.remove(os.path.abspath('temp.mol2'))
    return pmol.df.reset_index(drop=True)

def pmbol_list_to_mol2_df(pbmol_list):
    mol2_list=[pbmol_to_mol2_df(pbmol) for pbmol in pbmol_list]
    return mol2_list

def get_pbmol_list_from_string(molecule_string_list, molecule_format='xyz'):
    pbmol_list=[pybel.readstring(molecule_format, molecule_string) for molecule_string in molecule_string_list]
    return pbmol_list
def get_atom_type_list(mol2_df_list):
    atom_type_list=[mol2['atom_type'] for mol2 in mol2_df_list]
    return atom_type_list

def get_bonds_df(pbmol):
    data_string = pbmol.write('mol2')
    bonds = data_string.split('BOND')[1]
    bonds_array = np.array(bonds.split()).reshape(-1, 4)
    bonds_df = pd.DataFrame(bonds_array[:, 1:3])
    return bonds_df.astype(int)
def get_bonds_list(pbmol_list):
    bonds_list=[get_bonds_df(pbmol) for pbmol in pbmol_list]
    return bonds_list

def get_coordinates_array_list(mol2_df_list):
    coordinates_array_list=[np.array(coordinates[['x', 'y', 'z']].astype(float)) for coordinates in mol2_df_list]
    return coordinates_array_list


def array_list_to_dict(coordinates_array_list):
    dicts = {}
    for idx, array in enumerate(coordinates_array_list):
        dicts['conformer_num_{}'.format(idx)] = {'conformer_coordinates': array}
    return dicts

def coordinates_dfs_for_analyze(mol2_df_list):
    coordinates_df_list = [coordinates[['atom_type', 'x', 'y', 'z']] for coordinates in mol2_df_list]
    for df in coordinates_df_list:
        df['atom_type'].replace(GeneralConstants.REGULAR_BOND_TYPE.value, inplace=True)
        df.rename(columns={'atom_type': 'atom'}, inplace=True)
        df.reset_index(drop=True)
    return coordinates_df_list

def get_dipole_list(coordinates_array_list, charges_list):
    dipole_list=pd.concat([calc_dipole_charges(coordinates, charges) for coordinates, charges in
               zip(coordinates_array_list, charges_list)]).reset_index(drop=True)
    return dipole_list

#TODO: fix the problam with direction so its automatic
def get_sterimol_list(coordinates_df_list, bonds_list, atom_type_list, direction=None):
    # if direction==None:
    #     direction=list(bonds_list[0].iloc[0].values)[::-1]
    sterimol_list=pd.concat(
        [get_sterimol(coordinates_df, bonds_df, atom_type, direction) for coordinates_df, bonds_df, atom_type in
         zip(coordinates_df_list, bonds_list, atom_type_list)]).reset_index(drop=True)
    return sterimol_list




def calculate_energy(xyz_coordinates):
    # Create a Psi4 molecule object from the coordinates
    molecule = psi4.geometry(xyz_coordinates)

    # Set up the calculation
    psi4.set_options({'basis': 'sto-3g', 'scf_type': 'pk'})

    # Run the energy calculation
    energy = psi4.energy('scf')

    return energy

def filter_unique(mols, crit=0.3):
    """
    Remove structures that are very similar.
        Remove unconverged structures.

    Arguments:
        mols: pybel mol objects
    Returns:
        unique_mols: unique molecules
    """

    # Remove similar structures
    unique_mols = []
    aligner = pybel.ob.OBAlign()
    for mol_i in mols:
        aligner.SetRefMol(mol_i.OBMol)
        unique = True
        for mol_j in unique_mols:
            aligner.SetTargetMol(mol_j.OBMol)
            aligner.Align()
            rmsd = aligner.GetRMSD()
            if rmsd < crit:
                unique = False
                break
        if unique:
            unique_mols.append(mol_i)
    return unique_mols
def get_rmsd_df(pbmol_list, ref_idx):
    rmsd_df = pd.DataFrame()
    aligner = pybel.ob.OBAlign()
    aligner.SetRefMol(pbmol_list[ref_idx].OBMol)
    rmsd = []
    for pbmol_i in pbmol_list:
        aligner.SetTargetMol(pbmol_i.OBMol)
        aligner.Align()
        rmsd.append(aligner.GetRMSD())
    rmsd_df['conformer_{}'.format(ref_idx)] = rmsd
    return rmsd_df


# def rmsd_between_pbmol_lists(pbmol_list_1,pbmol_list_2):
#     rmsd_df = pd.DataFrame()
#     aligner = pybel.ob.OBAlign()
#     if len(pbmol_list_1)<len(pbmol_list_2):
#     pbmol_list_2=
#         for index,pbmol_1 in enumerate(pbmol_list_1):
#             aligner.SetRefMol(pbmol_1.OBMol)
#             rmsd=[]
#             for pbmol_2 in  pbmol_list_2 :
#                 aligner.SetTargetMol(pbmol_2.OBMol)
#                 aligner.Align()
#                 rmsd.append(aligner.GetRMSD())
#             rmsd_df['conformer_{}'.format(index)]=rmsd
#     return rmsd_df
def calc_dipole_charges(coordinates_array, charges_array, sub_atoms=None):  ##added option for subunits
    """
    a function that recives coordinates and npa charges, transform the coordinates
    by the new base atoms and calculates the dipole in each axis

    Parameters
    ---------
    coordinates_array: np.array
        contains x y z atom coordinates

    charges_array: np.array
        array of npa charges
    optional-sub_atoms:list
        calculate npa charges from a set of sub_atoms instead of all atoms.
    Returns:
    -------
    dip_df=calc_npa_charges(coordinates_array,charges,base_atoms_indexes,sub_atoms)
    Output:
    dip_df : pd.DataFrame
        output:            dip_x     dip_y     dip_z     total
                       0  0.097437 -0.611775  0.559625  0.834831
    """
    if sub_atoms == None:
        dip_xyz = np.vstack([(row[0] * row[1]) for row in list(zip(coordinates_array, charges_array))])
    else:
        dip_xyz = np.vstack([(row[0] * row[1]) for row in
                             np.array(list(zip(coordinates_array, charges_array)), dtype=object)[sub_atoms, :]])
    dip_vector = np.sum(dip_xyz, axis=0)
    array_dipole = np.hstack([dip_vector, np.linalg.norm(dip_vector)])

    dip_df = pd.DataFrame(array_dipole, index=['dip_x', 'dip_y', 'dip_z', 'total_dipole']).T
    return dip_df['total_dipole']


def get_sterimol(coordinates_df, bonds_df, atom_type, base_atoms, radii='bondi'):
    coordinates_array = np.array(coordinates_df[['x', 'y', 'z']].astype(float))
    bonds_direction = direction_atoms_for_sterimol(bonds_df, base_atoms)
    transformed_coordinates = calc_coordinates_transformation(coordinates_array, base_atoms)
    connected_from_direction = get_molecule_connections(bonds_df, bonds_direction[0], bonds_direction[1])
    ## creating df of atoms and bonds
    bonded_atoms_df = get_specific_bonded_atoms_df(bonds_df, connected_from_direction, coordinates_df)
    ### filtering nof bonds and H bonds *****deleting H only from the left column
    edited_coordinates_df = filter_atoms_for_sterimol(bonded_atoms_df, coordinates_df)
    ##adding colums
    extended_df = get_extended_df_for_sterimol(edited_coordinates_df, atom_type, radii)
    ###calculations
    try:
        b1s, b1s_loc = get_b1s_list(extended_df)
    except ValueError:
        return pd.DataFrame([float('NaN'), float('NaN'), float('NaN')], index=['B1', 'B5', 'L']).T
    B1 = min(b1s[b1s >= 0])
    # loc_B1=max(b1s_loc[np.where(b1s[b1s>=0]==min(b1s[b1s>=0]))])
    B5 = max(extended_df['B5'].values)
    L = max(extended_df['L'].values + 0.4)
    # try:
    #     loc_B5=min(extended_df['y'].iloc[np.where(extended_df['B5'].values==B5)[0][0]])
    # except TypeError:
    #     loc_B5=min(extended_df['y'].iloc[np.where(extended_df['B5'].values==B5)])
    df = pd.DataFrame([B1, B5, L], index=['B1', 'B5', 'L'])
    return df.T
class Conformers():

    def __init__(self, conformers_filename, molecule_format='xyz'):
        # self.path=os.path.abspath(conformers_filename.split('.')[0])
        self.conformers_filename = conformers_filename
        # os.chdir(self.path)
        self.conformers_string_list = split_molecule_file(conformers_filename)
        self.molecule_format = molecule_format
        self.pbmol_list =get_pbmol_list_from_string(self.conformers_string_list, self.molecule_format)
        self.mol2_df_list = pmbol_list_to_mol2_df(self.pbmol_list)
        self.charge_list = get_pbmol_charge_list(self.pbmol_list)
        self.atom_type_list = get_atom_type_list(self.mol2_df_list)
        self.bonds_list = get_bonds_list(self.pbmol_list)
        self.coordinates_array_list = get_coordinates_array_list(self.mol2_df_list)
        self.coordinates_df_list = coordinates_dfs_for_analyze(self.mol2_df_list)
        self.conformers_dict = array_list_to_dict(self.coordinates_array_list)
        self.dipole_list = get_dipole_list(self.coordinates_array_list, self.charge_list)
        self.sterimol_list = get_sterimol_list(self.coordinates_df_list, self.bonds_list, self.atom_type_list, [2,1])
        self.energy_list = get_pbmol_energy_list(self.pbmol_list)
        self.conformers_summary_df =self.get_conformers_summary_df()

    def write_all_conformers(self):
        for index, conformer in enumerate(self.conformers_list):
            with open('{}_{}.{}'.format(self.conformers_filename.split('.')[0], index,
                                        self.conformers_filename.split('.')[1]), 'w') as xyz_file:
                xyz_file.write(conformer)

    def get_conformers_summary_df(self):
        summary_df=pd.concat(
            [pd.DataFrame(self.conformers_string_list, columns=['conformer']), self.energy_list, self.dipole_list,
             self.sterimol_list], axis=1)  # ,self.conformers_energy_list
        return summary_df


    def compare_specific_conformers(self):
        ref_conformer = self.coordinates_df_list[0]

        max_dip_index = \
        np.where(self.extended_conformers['total_dipole'] == self.extended_conformers['total_dipole'].max())[0][0]

        max_dip_conformer = self.coordinates_df_list[max_dip_index]
        # print(max_dip_conformer)
        min_dip_index = \
        np.where(self.extended_conformers['total_dipole'] == self.extended_conformers['total_dipole'].min())[0][0]
        min_dip_conformer = self.coordinates_df_list[min_dip_index]
        # print(min_dip_conformer)
        min_B1_index = np.where(self.extended_conformers['B1'] == self.extended_conformers['B1'].min())[0][0]
        # print(min_B1_index)
        min_B1_conformer = self.coordinates_df_list[min_B1_index]
        try:
            overlay_analyzer = tab_data.OverlayAnalyzer(parser=tab_data.set_overlay_parser(),
                                                        xyz_filenames=['parameters', 'comparison'],
                                                        xyz_dfs=[ref_conformer, max_dip_conformer, min_dip_conformer,
                                                                 min_B1_conformer], fit_mode='all', )
        except IndexError:
            return print('one or more conformers are the same')



    def get_parameter_index(self, parameter, index_num=None, ascending=True):

        if index_num == None:
            ordered_df = self.conformers_summary_df.sort_values(parameter, axis=0, ascending=ascending)
        else:
            ordered_df = self.conformers_summary_df.sort_values(parameter, axis=0, ascending=ascending).head(index_num)
        return ordered_df.index.values.tolist()


class ConformersFieldAnalayzer():

    def __init__(self, conformers_object):
        self.conformers_class = conformers_object
        self.num_conformers = len(self.conformers_class.conformers_string_list)
        self.rmsd_energy_dict = self.get_rmsd_dict('energy')

    def get_rmsd_dict(self, parameter='energy', ref_idx=0):
        order_mask = self.conformers_class.get_parameter_index(parameter)
        ordered_pbmol_list = np.array(self.conformers_class.pbmol_list)[order_mask]
        rmsd_df = get_rmsd_df(ordered_pbmol_list, ref_idx=0)
        rmsd_mean = rmsd_df.mean()
        rmsd_varience = rmsd_df.var(ddof=1)
        rsmd_dict = {'rmsd_df': rmsd_df, 'rmsd_mean': rmsd_mean, 'rmsd_varience': rmsd_varience}
        return rsmd_dict


# class ConformersCompare():

#     def __init__(self, conformers_field_1, conformers_field_2):

#         self.rmsd_diff_df=rmsd_between_pbmol_lists(conformers_field_1.conformers_class.pbmol_list, conformers_field_2.conformers_class.pbmol_list)

if __name__ == '__main__':
    rdkit_conformers = Conformers('rdkitconformers_gauss_compare.sdf', 'sdf')
    ob_conformers = Conformers('obconformers_gauss_compare.xyz', 'xyz')
    # crest_conformers=Conformers('crest_conformers.xyz')
    conformers_field_1 = ConformersFieldAnalayzer(rdkit_conformers)
    # conformers_field_2 = ConformersFieldAnalayzer(ob_conformers)

    # compare=ConformersCompare(conformers_field_1, conformers_field_2)

#     data=tab_data.TabDataAnalyzer(parser=tab_data.set_tab_parser(),origin_df=crest_conformers.coordinates_df_list[1],xyz_filename='conformer',get_plot=True)
#     # overlay_analyzer=tab_data.OverlayAnalyzer(parser=tab_data.set_overlay_parser(),xyz_filenames=['0','1'] ,xyz_dfs=[crest_conformers.coordinates_df_list[0], crest_conformers.coordinates_df_list[1],crest_conformers.coordinates_df_list[2]], fit_mode='all',)
# ##                                         atoms=['30', '20', '16', '0', '24', '8', '5', '33'],
# 'all'
