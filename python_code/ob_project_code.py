from openbabel import pybel
from openbabel import openbabel as ob
import help_functions
import os
import numpy as np
import pandas as pd
from biopandas.mol2 import PandasMol2
from modified_project import *
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
path_to_add=r'C:\Users\edens\Documents\GitHub\conformers_project\python_code\molecules'
os.chdir(path_to_add)

  
def obmol_to_3d(obmol,output_name):
    pbmol = pybel.Molecule(obmol)
    pbmol.addh()
    pbmol.make3D(forcefield="gaff",steps=100)
    pbmol.localopt()
    pbmol.write('pdb',help_functions.change_filetype(output_name,'pdb'),overwrite=True)
    return 

def obmol_from_coordinates(file_name,output_type='mol2'):
    obConversion = ob.OBConversion()
    file_type=file_name.split('.')[-1]
    obConversion.SetInAndOutFormats(file_type, output_type)
    obmol = ob.OBMol()
    obConversion.ReadFile(obmol, file_name)
    return obmol

def create_3d_from_smiles(smiles_code,out_put_name): #'O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N'
    obmol = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInFormat('smi')
    conv.ReadString(obmol, smiles_code)
    return obmol_to_3d(obmol,out_put_name)

def create_3d_from_file(file_name):
    file_type=file_name.split('.')[-1]
    obmol = obmol_from_coordinates(file_type,file_name)
    return obmol_to_3d(obmol,help_functions.change_filetype(file_name,'pdb'))

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

def calc_energy(obmol):
    pff = ob.OBForceField_FindType( "mmff94" )
    pff.SetLogLevel(ob.OBFF_LOGLVL_HIGH)
    pff.SetLogToStdErr()
    pff.Setup(obmol)
    pff.SetLogLevel(ob.OBFF_LOGLVL_NONE)
    return pff.Energy()

def calc_dipole_charges(coordinates_array,charges_array,base_atoms_indexes,sub_atoms=None):##added option for subunits
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
    transformed_coordinates_array=calc_coordinates_transformation(coordinates_array,base_atoms_indexes)
    if sub_atoms==None:
        dip_xyz=np.vstack([(row[0]*row[1]) for row in list(zip(transformed_coordinates_array,charges_array))])
    else:
        dip_xyz=np.vstack([(row[0]*row[1]) for row in np.array(list(zip(transformed_coordinates_array,charges_array)),dtype=object)[sub_atoms,:]])
    dip_vector=np.sum(dip_xyz,axis=0)
    array_dipole=np.hstack([dip_vector,np.linalg.norm(dip_vector)])
    dip_df=pd.DataFrame(array_dipole,index=['dip_x','dip_y','dip_z','total']).T
    return dip_df   

               
def get_obmol_charge(obmol):
    ob_charge_model = ob.OBChargeModel.FindType("eem2015bn")
    ob_charge_model.ComputeCharges(obmol)
    return ob_charge_model.GetPartialCharges()

# def get_molecule_dipole(mol):
#     ob_charge_model = ob.OBChargeModel.FindType("eem2015bn")
#     return ob_charge_model.GetDipoleMoment(mol) ##cant open vector3 object

# def get_molecule_angle_data(mol):
#     ob_angle_data=ob.OBAngleData()
#     return ob_angle_data.GetData(mol)


  
def get_molecule_string(obmol,output_type='mol2'):
    obconversion = ob.OBConversion()
    obconversion.SetOutFormat(output_type)
    return obconversion.WriteString(obmol)

def split_molecule_file(file_name): ##for crest output
    file_type=file_name.split('.')[-1]
    molecules=[pbmol.write(file_type) for pbmol in pybel.readfile(file_type,file_name)]
    return molecules

        #mol.write("sdf", "%s.sdf" % mol.title) to create file
        
def molecules_file_to_3d(file_name): ##from a file containing many molecules
    file_type=file_name.split('.')[-1]
    for pbmol in pybel.readstring(file_type,file_name):
        obmol=pbmol.OBMol
        obmol_to_3d(obmol,"%s.pdb" % pbmol.title)

def obmol_to_mol2_df(pbmol):
    pbmol.write('mol2','temp.mol2')
    pmol = PandasMol2()
    pmol.read_mol2(os.path.abspath('temp.mol2'))
    os.remove(os.path.abspath('temp.mol2'))
    return pmol.df


def get_pbmol_bonds_df(pbmol) :
    data_string=pbmol.write('mol2')
    bonds=data_string.split('BOND')[1]
    bonds_array=np.array(bonds.split()).reshape(-1,4)
    bonds_df=pd.DataFrame(bonds_array[:,1:3])
    return bonds_df.astype(int)


def mol_string_list_to_3d(mol_format,mol_string_list):
    """
    Parameters
    ----------
    mol_format : str
        type of molecule input -mol2 xyz smi ...
    mol_string_list : list 
        a list of molecule strings as recived from confab_search.

    Returns
    -------
    create 3d file stractures

    """
    obconversion = ob.OBConversion()
    obconversion.SetInFormat(mol_format)
    counter=0
    for mol in mol_string_list:
        counter+=1
        obmol = ob.OBMol()
        obconversion.ReadString(obmol, mol)
        obmol_to_3d(obmol,"%s.pdb" % counter )

def RMSD_between_mol2_dfs(mol2_df_1,mol2_df_2):
    r_heavy = PandasMol2.rmsd(mol2_df_1, mol2_df_2)
    r_all  = PandasMol2.rmsd(mol2_df_1, mol2_df_2, heavy_only=False)
    print('Heavy-atom RMSD: {.4f} Angstrom' .format(r_heavy))
    print('All-atom RMSD: {.4f} Angstrom'.format( r_all))

def RMSD_between_mol2_lists(mol2_list_1,mol2_list_2):
    for mol1 in mol2_list_1:
        for mol2 in mol2_list_2:
            RMSD_between_mol2_dfs(mol1,mol2)

def RMSD_between_obmol(obmol_1,obmol_2):
    obalign=ob.OBAlign(obmol_1,obmol_2)
    obalign.Align()
    return obalign.GetRMSD()


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

def sort_output(pbmols, output_name):
    """Sort molecules based on energies
    
    Arguments:
        pbmol: a list of pybel mol objects
        output_name: a SDF file to store the molecules
    """
    l = []
    for pbmol in pbmols:
        name = pbmol.data['ID'].split('_')[0].strip()
        e = pbmol.data['E_tot']
        e_rel = pbmol.data['E_relative']

        l.append((name, e_rel, pbmol))

    l = sorted(l)

    f = pybel.Outputfile('sdf', output_name)
    for n_er_m in l:
        name, e_relative, mol = n_er_m
        f.write(mol)
    f.close()
    
# def unique_conformers(files, crit=0.5): #need editing
#     """Removing conformers whose RMSD is within crit
    
#     Arguments:
#         files: sdf files
#     """
#     unique_files = []
#     duplicate_files = []
#     aligner = pybel.ob.OBAlign()
#     for file in files:
#         mol = next(pybel.readfile("sdf", file))
#         aligner.SetRefMol(mol.OBMol)
#         unique = True
#         for f in unique_files:
#             mol_j = next(pybel.readfile("sdf", f))
#             aligner.SetTargetMol(mol_j.OBMol)
#             aligner.Align()
#             rmsd = aligner.GetRMSD()
#             if rmsd < crit:
#                 unique = False
#                 break
#         if unique:
#             unique_files.append(file)
#         else:
#             duplicate_files.append(file)
#     c = len(unique_files) + len(duplicate_files)
#     assert(c == len(files))
#     for file in duplicate_files:
#         os.remove(file)
    
def get_sterimol(coordinates_array,bonds_df,atype,base_atoms,radii='bondi'):
 
    bonds_direction=direction_atoms_for_sterimol(bonds_df,base_atoms)
    transformed_coordinates=calc_coordinates_transformation(coordinates_array,base_atoms)
    connected_from_direction=get_molecule_connections(coordinates_array,bonds_direction[0],bonds_direction[1])
    ## creating df of atoms and bonds
    bonded_atoms_df=get_specific_bonded_atoms_df(bonds_df,connected_from_direction,coordinates_array)     
      ### filtering nof bonds and H bonds *****deleting H only from the left column       
    edited_coordinates_df=filter_atoms_for_sterimol(bonded_atoms_df,coordinates_df)
    ##adding colums
    extended_df=get_extended_df_for_sterimol(edited_coordinates_df,atype,radii)
    ###calculations
    b1s,b1s_loc=get_b1s_list(extended_df)
    B1=min(b1s[b1s>=0])
    loc_B1=max(b1s_loc[np.where(b1s[b1s>=0]==min(b1s[b1s>=0]))])
    B5=max(extended_df['Bs'].values)
    L=max(extended_df['L'].values+0.4)
    loc_B5=min(extended_df['y'].iloc[np.where(extended_df['Bs'].values==B5)[0][0]])
    df=pd.DataFrame([B1,B5,L,loc_B5,loc_B1],index=[['B1','B5','L','loc_B5','loc_B1']]) 
    
class Conformers():
    
    def __init__(self,conformer_list,molecule_format='xyz'):
        
        self.comformers_list=conformer_list
        self.molecule_format=molecule_format
        self.pbmol_list=[pybel.readstring(molecule_format, molecule_string) for molecule_string in conformer_list]
        self.obmol_list=[pbmol.OBMol for pbmol in self.pbmol_list]
        self.mol2_df_list=[obmol_to_mol2_df(pbmol) for pbmol in self.pbmol_list]
        self.charges_list=[get_obmol_charge(obmol) for obmol in self.obmol_list]
        self.atype_list=[mol2['atom_type'] for mol2 in self.mol2_df_list]
        # mol_string_list_to_3d(self.molecule_format,self.comformer_list)
        self.conformers_bonds_df_list=[get_pbmol_bonds_df(pbmol) for pbmol in self.pbmol_list]
        self.coordinates_array_list=[np.array(coordinates[['x','y','z']].astype(float)) for coordinates in self.mol2_df_list]
        self.dipole_list=[calc_dipole_charges(coordinates,charges,[1,2,3]) for coordinates,charges in zip(self.coordinates_array_list,self.charges_list)]
        # self.sterimol_list=[get_sterimol(coordinates_array,bonds_df,atype,[2,1]) for coordinates_array,bonds_df,atype in zip(self.coordinates_array_list, self.conformers_bonds_df_list,self.atype_list)]
        self.conformers_energy_list=[calc_energy(obmol) for obmol in self.obmol_list]
        
class Molecule():
    """
    """
    def __init__(self,molecule_dir_name): 
        """
        

        Parameters:
        ----------
        molecule_dir_name : str
            the name of the molecule directory.

        Description:
        -------
  

        """
        self.molecule_path=os.path.abspath(molecule_dir_name)
        os.chdir(self.molecule_path)
        self.origin_file_name=(os.listdir(self.molecule_path)[0])
        self.obmol=obmol_from_coordinates(self.origin_file_name)
        self.pbmol=pybel.Molecule(self.obmol)
        self.mol2_df=obmol_to_mol2_df(self.pbmol)
        self.charges=get_obmol_charge(self.obmol)
        
        os.chdir('../')
     
if __name__ == '__main__':        
    # help_functions.create_molecule_directories()
    # help_functions.move_files_directory('.xyz')
    
    os.chdir(r'C:\Users\edens\Documents\GitHub\conformers_project\python_code\molecules')
    # # molecule_1=Molecule('ABAKAR')
    # molecule_2=Molecule('gauss_compare') ## only 2 conformers
    # # mol2_df=get_molecule_string(molecule_2.obmol)
    molecule_3=Molecule('origin_molecule')
    x=confab_search(molecule_3.obmol)
    #creates conformers and turns each one to a 3d model
    # rdkit_conformers_list=split_molecule_file('rdkitconformers_origin_molecule.sdf')[0:5]
    # rdkit_conformers=Conformers(rdkit_conformers_list,'sdf')
    # mol_string_list_to_3d('sdf',rdkit_conformers_list)
  
    ob_conformers=confab_search(molecule_3.obmol)
    ob_conformers=Conformers(x[0:10])
    
    
    # balloon_conformers_list=split_molecule_file('balloon_coor_conformers.xyz')
    # balloon_conformers=Conformers(balloon_conformers_list)
    # crest_conformers_list=split_molecule_file('crest_conformers.xyz')[0:5]
    # crest_conformers=Conformers(crest_conformers_list)
    
    # RMSD_between_mol2_lists(crest_conformers.mol2_df_list,ob_conformers.mol2_df_list)
    
    # mol = create_3d_from_smiles(('O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N'),'3d.pdb')
    # mol_1=create_3d_from_xyz('conf_pyt.xyz','1.pdb')
    # mol_2=create_3d_from_xyz('conf_bash.xyz','2.pdb')

