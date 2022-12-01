from openbabel import pybel
from openbabel import openbabel as ob
import help_functions
import os
from enum import Enum
import numpy as np
import pandas as pd
from biopandas.mol2 import PandasMol2
from modified_project import *
import tab_data

"""
pybel ob diffrences

## pbmol -- pybel molecule // obmol--openbabel molecule
##pybelmol = pybel.Molecule(mol) # to go from OBMol to pybel
##obmol=mol.OBMol # to go from pybel to open babel  

molecule.obmol.GetFormula()
Out[6]: 'C132H112N12P4S4'

molecule.pbmol.formula
Out[7]: 'C132H112N12P4S4'

  """
# path_to_add=r'C:\Users\edens\Documents\GitHub\conformers_project\python_code\molecules'
# os.chdir(path_to_add)

 
class GeneralConstants(Enum):
    """
    Holds constants for calculations and conversions
    1. covalent radii from Alvarez (2008) DOI: 10.1039/b801115j
    2. atomic numbers
    2. atomic weights
    """
    COVALENT_RADII= {
            'H': 0.31, 'He': 0.28, 'Li': 1.28,
            'Be': 0.96, 'B': 0.84, 'C': 0.76, 
            'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
            'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
            'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
            'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
            'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
            'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
            'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
            'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
            'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
            'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
            'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
            'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
            'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
            'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
            'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
            'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
            'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
            'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
            'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
            'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
            'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
            'Am': 1.80, 'Cm': 1.69
    }
    
    CPK_RADII={
        'H': 1.10, 'C': 1.70, 'F': 1.47,
        'S': 1.80, 'B': 1.92, 'I': 1.98, 
        'N': 1.55, 'O': 1.52, 'Co': 2.00, 
        'Br': 1.83, 'Si': 2.10,'Ni': 2.00,
        'P': 1.80, 'Cl': 1.75, 
    }
    
    BONDI_RADII={
        'C':1.50,   'H':1.00,   'S.O':1.70,  'Si':2.10,
        'C2':1.60,  'N':1.50,   'S1':1.00,   'Co':2.00,
        'C3':1.60,  'C66':1.70, 'F':1.35,    'Ni':2.00,
        'C4':1.50,  'N4':1.45,  'Cl':1.80,
        'C5/N5':1.70, 'O':1.35, 'S4':1.40,
        'C6/N6':1.70, 'O2':1.35, 'Br':1.95,
        'C7':1.70,    'P':1.40,  'I':2.15,
        'C8':1.50,    'S':1.70,  'B':1.92,
    
    }
    
    BOND_TYPE={
        
        'O.2':'O2', 'N.2':'C6/N6','S.3':'S4',
        'O.3':'O', 'N.1':'N', 'S.O2':'S1',
        'O.co2':'O', 'N.3':'C6/N6','P.3':'P',
        'C.1':'C3', 'N.ar':'C6/N6',
        'C.2':'C2', 'N.am':'C6/N6',
        "C.cat":'C3', 'N.pl3':'C6/N6',
        'C.3':'C', 'N.4':'N',
        'C.ar':'C6/N6', 'S.2':'S',
        }
    REGULAR_BOND_TYPE={
        
        'O.2':'O', 'N.2':'N','S.3':'S',
        'O.3':'O', 'N.1':'N', 'S.O2':'S',
        'O.co2':'O', 'N.3':'N','P.3':'P',
        'C.1':'C', 'N.ar':'N',
        'C.2':'C', 'N.am':'N',
        "C.cat":'C', 'N.pl3':'N',
        'C.3':'C', 'N.4':'N',
        'C.ar':'C', 'S.2':'S',
        }
    
    ATOMIC_NUMBERS ={
    '1':'H', '5':'B', '6':'C', '7':'N', '8':'O', '9':'F', '14':'Si',
             '15':'P', '16':'S', '17':'Cl', '35':'Br', '53':'I', '27':'Co', '28':'Ni'}
        

    ATOMIC_WEIGHTS = {
            'H' : 1.008,'He' : 4.003, 'Li' : 6.941, 'Be' : 9.012,
            'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,
            'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,
            'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,
            'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,
            'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,
            'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,
            'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,
            'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,
            'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,
            'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,
            'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,
            'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,
            'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,
            'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,
            'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'Gd' : 157.25,
            'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,
            'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,
            'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,
            'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,
            'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,
            'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,
            'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,
            'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247
    }

 


def pbmol_from_coordinates(file_name,output_type='mol2'):
    file_type=file_name.split('.')[-1]
    pbmol=pybel.readfile(file_type, file_name)
    return pbmol

def split_molecule_file(file_name): ##for crest output
    file_type=file_name.split('.')[-1]
    molecules=[pbmol.write(file_type) for pbmol in pybel.readfile(file_type,file_name)]
    return molecules

def smile_to_coordinates(smile):
    pbmol=pybel.readstring('smi', smile)
    pbmol.make3D(forcefield="gaff",steps=100)
    pbmol.localopt()    
    return pbmol.write('pdb')

def smiles_to_coordinates(smiles):
    return [smile_to_coordinates(smile) for smile in smiles]


def freeze_atoms_for_confab(obmol,atoms_to_freeze):
    constraints = ob.OBFFConstraints()
    for atom in ob.OBMolAtomIter(obmol):
        atom_id = atom.GetIndex() 
        if atom_id in atoms_to_freeze:
            constraints.AddAtomConstraint(atom_id)
    return constraints

def confab_search(obmol,set_constraints=False,atoms_to_freeze=None,output_format='xyz'):
    pff = ob.OBForceField_FindType( "mmff94" )
    if set_constraints :
        constraints=freeze_atoms_for_confab(obmol,atoms_to_freeze)
        pff.SetConstraints(constraints)
    pff.DiverseConfGen(0.5, 1000, 50.0, True) #allow change
    pff.Setup(obmol)
    pff.GetConformers(obmol)
    obconversion = ob.OBConversion()
    obconversion.SetOutFormat(output_format)
    output_strings = []
    for conf_num in range(obmol.NumConformers()):
        obmol.SetConformer(conf_num)
        output_strings.append(obconversion.WriteString(obmol))
    return output_strings

# def calc_energy(obmol):
#     pff = ob.OBForceField_FindType( "mmff94" )
#     # pff.SetLogLevel(ob.OBFF_LOGLVL_HIGH)
#     pff.SetLogToStdErr()
#     pff.Setup(obmol)
#     # pff.SetLogLevel(ob.OBFF_LOGLVL_NONE)
#     return pff.Energy()

def calc_dipole_charges(coordinates_array,charges_array,sub_atoms=None):##added option for subunits
    """
    a function that recives coordinates and npa charges, transform the coordinates
    by the new base atoms and calculates the dipole in each axis
    
    Parameters
    ---------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    charges_array: np.array
        array of npa charges
    base_atoms_indexes:list
        3/4 atom indexes for coordinates transformation
        
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
    if sub_atoms==None:
        dip_xyz=np.vstack([(row[0]*row[1]) for row in list(zip(coordinates_array,charges_array))])
    else:
        dip_xyz=np.vstack([(row[0]*row[1]) for row in np.array(list(zip(coordinates_array,charges_array)),dtype=object)[sub_atoms,:]])
    dip_vector=np.sum(dip_xyz,axis=0)
    array_dipole=np.hstack([dip_vector,np.linalg.norm(dip_vector)])
    
    dip_df=pd.DataFrame(array_dipole,index=['dip_x','dip_y','dip_z','total_dipole']).T
    return dip_df['total_dipole']

               
def get_obmol_charge(obmol):
    ob_charge_model = ob.OBChargeModel.FindType("eem2015bn")
    ob_charge_model.ComputeCharges(obmol)
    return ob_charge_model.GetPartialCharges()





def filter_unique(mols, crit=0.3):
    """Remove structures that are very similar.
        Remove unconverged structures.
    
    Arguments:
        mols: pybel mol objects
    Returns:
        unique_mols: unique molecules
    """

    #Remove similar structures
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
 



def pbmol_to_mol2_df(pbmol):
    pbmol.write('mol2','temp.mol2')
    pmol = PandasMol2()
    pmol.read_mol2(os.path.abspath('temp.mol2'))
    os.remove(os.path.abspath('temp.mol2'))
    return pmol.df.reset_index(drop=True)


def get_pbmol_bonds_df(pbmol) :
    data_string=pbmol.write('mol2')
    bonds=data_string.split('BOND')[1]
    bonds_array=np.array(bonds.split()).reshape(-1,4)
    bonds_df=pd.DataFrame(bonds_array[:,1:3])
    return bonds_df.astype(int)


def get_rmsd_df(pbmol_list):
    rmsd_df = pd.DataFrame()
    aligner = pybel.ob.OBAlign()
    for index,pbmol_i in enumerate(pbmol_list):
        aligner.SetRefMol(pbmol_i.OBMol)
        rmsd=[]
        for pbmol_j in pbmol_list:
            aligner.SetTargetMol(pbmol_j.OBMol)
            aligner.Align()
            rmsd.append(aligner.GetRMSD())
        rmsd_df['conformer_{}'.format(index)]=rmsd
    return rmsd_df

def RMSD_between_pbmol_lists(pbmol_list_1,pbmol_list_2):
    rmsd_df = pd.DataFrame()
    aligner = pybel.ob.OBAlign()
    for index,pbmol_1 in enumerate(pbmol_list_1):
        aligner.SetRefMol(pbmol_1.OBMol)  
        rmsd=[]
        for pbmol_2 in  pbmol_list_2 :
            aligner.SetTargetMol(pbmol_2.OBMol)
            aligner.Align()
            rmsd.append(aligner.GetRMSD())
        rmsd_df['conformer_{}'.format(index)]=rmsd
    return rmsd_df
            
def atomType(pbmol, atomIdx):
    """get the atomic type given an atom index, both in pybel mol object"""
    atom_num = pbmol.OBMol.GetAtom(atomIdx).GetAtomicNum()
    return atom_num

def check_bonds(pbmol):
    """Check if a pybel mol object has valid bond lengths"""
    # Initialize UFF bond radii (Rappe et al. JACS 1992)
    # Units of angstroms 
    # These radii neglect the bond-order and electronegativity corrections in the original paper. Where several values exist for the same atom, the largest was used. 
    Radii = {1:0.354, 
             5:0.838, 6:0.757, 7:0.700,  8:0.658,  9:0.668,
             14:1.117, 15:1.117, 16:1.064, 17:1.044,
             32: 1.197, 33:1.211, 34:1.190, 35:1.192,
             51:1.407, 52:1.386,  53:1.382}

    for bond in ob.OBMolBondIter(pbmol.OBMol):
        length = bond.GetLength()
        begin = atomType(pbmol, bond.GetBeginAtomIdx())
        end = atomType(pbmol, bond.GetEndAtomIdx())
        reference_length = (Radii[begin] + Radii[end]) * 1.25
        if length > reference_length:
            return False
    return True

    

    
def get_sterimol(coordinates_df,bonds_df,atype,base_atoms,radii='bondi'):
    coordinates_array=np.array(coordinates_df[['x','y','z']].astype(float))
    bonds_direction=direction_atoms_for_sterimol(bonds_df,base_atoms)
    transformed_coordinates=calc_coordinates_transformation(coordinates_array,base_atoms)
    connected_from_direction=get_molecule_connections(bonds_df,bonds_direction[0],bonds_direction[1])
    ## creating df of atoms and bonds
    bonded_atoms_df=get_specific_bonded_atoms_df(bonds_df,connected_from_direction,coordinates_df)     
      ### filtering nof bonds and H bonds *****deleting H only from the left column       
    edited_coordinates_df=filter_atoms_for_sterimol(bonded_atoms_df,coordinates_df)
    ##adding colums
    extended_df=get_extended_df_for_sterimol(edited_coordinates_df,atype,radii)
    ###calculations
    try:
        b1s,b1s_loc=get_b1s_list(extended_df)
    except ValueError:
        return pd.DataFrame([float('NaN'),float('NaN'),float('NaN')],index=['B1','B5','L']).T
    B1=min(b1s[b1s>=0])
    # loc_B1=max(b1s_loc[np.where(b1s[b1s>=0]==min(b1s[b1s>=0]))])
    B5=max(extended_df['Bs'].values)
    L=max(extended_df['L'].values+0.4)
    # try:
    #     loc_B5=min(extended_df['y'].iloc[np.where(extended_df['Bs'].values==B5)[0][0]])
    # except TypeError:
    #     loc_B5=min(extended_df['y'].iloc[np.where(extended_df['Bs'].values==B5)])
    df=pd.DataFrame([B1,B5,L],index=['B1','B5','L']) 
    return df.T
  
# """

# """
# from pyscf import gto, scf
# mol = gto.M(atom="ethene.xyz")

# class ConformerCompare():
    
#     def __init__(self,compare_list):
#         self.compare_list=compare_list
#         self.ob_conformers,self.rdkit_conformers,self.crest_conformers=(*compare_list)
    
#     def sort_conformers_by_parameter(self,parameter,ascending=True):
#         sorted_compare_list=[]
#         for conformers in compare_list:
#             sortedconformers.extended_conformers.sort_values(by=parameter)
            
#     def compare_conformers_by_parameter(self,parameter,ascending=True):
#         rmsd_list=[]
       
#             for conformers_j in compare_list:
        
class Conformers():
    
    def __init__(self,conformers_filename,molecule_format='xyz'):
        # self.path=os.path.abspath(conformers_filename.split('.')[0])
        self.conformers_filename=conformers_filename
        # os.chdir(self.path)
        self.conformers_list=split_molecule_file(conformers_filename)
        self.molecule_format=molecule_format
        self.pbmol_list=[pybel.readstring(molecule_format, molecule_string) for molecule_string in self.conformers_list]
        self.obmol_list=[pbmol.OBMol for pbmol in self.pbmol_list]
        self.mol2_df_list=[pbmol_to_mol2_df(pbmol) for pbmol in self.pbmol_list]
        self.charges_list=[get_obmol_charge(obmol) for obmol in self.obmol_list]
        self.atype_list=[mol2['atom_type'] for mol2 in self.mol2_df_list]
        self.conformers_bonds_df_list=[get_pbmol_bonds_df(pbmol) for pbmol in self.pbmol_list]
        self.coordinates_array_list=[np.array(coordinates[['x','y','z']].astype(float)) for coordinates in self.mol2_df_list]
        self.coordinates_df_list=self.coordinates_dfs_for_analyze()
        self.conformers_dict=self.array_list_to_dict()
        self.dipole_list=pd.concat([calc_dipole_charges(coordinates,charges) for coordinates,charges in zip(self.coordinates_array_list,self.charges_list)]).reset_index(drop=True)
        self.sterimol_list=pd.concat([get_sterimol(coordinates_df,bonds_df,atype,[2,1]) for coordinates_df,bonds_df,atype in zip(self.coordinates_df_list, self.conformers_bonds_df_list,self.atype_list)]).reset_index(drop=True)
       # self.conformers_energy_list=pd.DataFrame([calc_energy(obmol) for obmol in self.obmol_list],columns=['energy'])
        self.extended_conformers=pd.concat([pd.DataFrame(self.conformers_list,columns=['conformer']),self.dipole_list,self.sterimol_list],axis=1) #,self.conformers_energy_list
        
    def write_all_conformers(self):
        for index,conformer in enumerate(self.conformers_list):
            with open('{}_{}.{}'.format(self.conformers_filename.split('.')[0],index,self.conformers_filename.split('.')[1]), 'w') as xyz_file:
                xyz_file.write(conformer)
                
    def coordinates_dfs_for_analyze(self)    :
        coordinates_df_list=[coordinates[['atom_type','x','y','z']] for coordinates in self.mol2_df_list]
        for df in coordinates_df_list:
            df['atom_type'].replace(GeneralConstants.REGULAR_BOND_TYPE.value, inplace=True)
            df.rename(columns={'atom_type': 'element'},inplace=True)
            df.reset_index(drop=True)
        return coordinates_df_list
           
    def compare_specific_conformers(self):
        ref_conformer=self.coordinates_df_list[0]
        
        max_dip_index=np.where(self.extended_conformers['total_dipole']==self.extended_conformers['total_dipole'].max())[0][0]
      
        max_dip_conformer=self.coordinates_df_list[max_dip_index]
        # print(max_dip_conformer)
        min_dip_index=np.where(self.extended_conformers['total_dipole']==self.extended_conformers['total_dipole'].min())[0][0]
        min_dip_conformer=self.coordinates_df_list[min_dip_index]
        # print(min_dip_conformer)
        min_B1_index=np.where(self.extended_conformers['B1']==self.extended_conformers['B1'].min())[0][0]
        # print(min_B1_index)
        min_B1_conformer=self.coordinates_df_list[min_B1_index]
        try:
            overlay_analyzer=tab_data.OverlayAnalyzer(parser=tab_data.set_overlay_parser(),xyz_filenames=['parameters','comparison'] ,xyz_dfs=[ref_conformer,max_dip_conformer, min_dip_conformer,min_B1_conformer], fit_mode='all',)
        except IndexError:
            return print('one or more conformers are the same')
        
    def array_list_to_dict(self):
        dicts={}
        
        for idx,array in enumerate(self.coordinates_array_list):
            dicts['conformer_num_{}'.format(idx)]={'conformer_coordinates':array}
            
        return dicts

        
        os.chdir('../')
     
if __name__ == '__main__':
            
    # help_functions.create_molecule_directories()
    # help_functions.move_files_directory('.xyz')
    
    # os.chdir(r'C:\Users\edens\Documents\GitHub\conformers_project\python_code\molecules\crest_conformers')
    # # balloon_conformers_list=split_molecule_file('balloon_coor_conformers.xyz')
    # # balloon_conformers=Conformers(balloon_conformers_list)
    # rdkit_conformers=Conformers('rdkitconformers_origin_molecule.sdf','sdf')
    # # ob_conformers=Conformers('obconformers_origin_molecule.sdf','sdf')
    crest_conformers=Conformers('crest_conformers.xyz')
    # conformers_list=[rdkit_conformers, ob_conformers, crest_conformers]
#     data=tab_data.TabDataAnalyzer(parser=tab_data.set_tab_parser(),origin_df=crest_conformers.coordinates_df_list[1],xyz_filename='conformer',get_plot=True)
#     # overlay_analyzer=tab_data.OverlayAnalyzer(parser=tab_data.set_overlay_parser(),xyz_filenames=['0','1'] ,xyz_dfs=[crest_conformers.coordinates_df_list[0], crest_conformers.coordinates_df_list[1],crest_conformers.coordinates_df_list[2]], fit_mode='all',)
# ##                                         atoms=['30', '20', '16', '0', '24', '8', '5', '33'],
                                       #'all'
