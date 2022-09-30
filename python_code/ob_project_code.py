from openbabel import pybel
import openbabel as ob
import help_functions
import os
import numpy as np
import pandas as pd
from biopandas.mol2 import PandasMol2
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

  
def molecule_to_3d(obmol,output_name):
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
    return molecule_to_3d(obmol,out_put_name)

def create_3d_from_file(file_name):
    file_type=file_name.split('.')[-1]
    obmol = obmol_from_coordinates(file_type,file_name)
    return molecule_to_3d(obmol,help_functions.change_filetype(file_name,'pdb'))

def confab_search(obmol,output_format='xyz'):
    pff = ob.OBForceField_FindType( "mmff94" )
    pff.DiverseConfGen(0.5, 10000, 50.0, True) #allow change
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
    return pff.Energy()
                   
def get_molecule_charge(obmol):
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
        molecule_to_3d(obmol,"%s.pdb" % pbmol.title)

def obmol_to_mol2_df(pbmol):
    pbmol.write('mol2','temp.mol2')
    pmol = PandasMol2()
    pmol.read_mol2(os.path.abspath('temp.mol2'))
    os.remove(os.path.abspath('temp.mol2'))
    return pmol.df
    
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
        molecule_to_3d(obmol,"%s.pdb" % counter )

def RMSD_between_mol2_dfs(mol2_df_1,mol2_df_2):
    r_heavy = PandasMol2.rmsd(mol2_df_1, mol2_df_2)
    r_all  = PandasMol2.rmsd(mol2_df_1, mol2_df_2, heavy_only=False)
    print('Heavy-atom RMSD: %.4f Angstrom' % r_heavy)
    print('All-atom RMSD: %.4f Angstrom' % r_all)

def RMSD_between_mol2_lists(mol2_list_1,mol2_list_2):
    for mol1 in mol2_list_1:
        for mol2 in mol2_list_2:
            RMSD_between_mol2_dfs(mol1,mol2)
    
    
class Conformers():
    
    def __init__(self,conformer_list,molecule_format='xyz'):
        
        self.comformer_list=conformer_list
        self.molecule_format=molecule_format
        self.pbmol_list=[pybel.readstring(molecule_format, molecule_string) for molecule_string in conformer_list]
        self.mol2_df_list=[obmol_to_mol2_df(pbmol) for pbmol in self.pbmol_list]
        # mol_string_list_to_3d(self.molecule_format,self.comformer_list)
        # self.conformer_energy_list=[calc_energy(obmol) for obmol in self.pbmol_list]
        
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
        self.charges=get_molecule_charge(self.obmol)
        
        os.chdir('../')
     
        
# help_functions.create_molecule_directories()
# help_functions.move_files_directory('.xyz')


molecule_1=Molecule('ABAKAR')
molecule_2=Molecule('gauss_compare') ## only 2 conformers
molecule_3=Molecule('origin_molecule')




#creates conformers and turns each one to a 3d model
x=confab_search(molecule_3.obmol)
# mol_string_list_to_3d('xyz',x)

# ob_conformers=Conformers(x)
# crest_conformers_list=split_molecule_file('crest_conformers.xyz')[0:5]
# crest_conformers=Conformers(crest_conformers_list)

# RMSD_between_mol2_lists(crest_conformers.mol2_df_list,ob_conformers.mol2_df_list)

# mol = create_3d_from_smiles(('O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N'),'3d.pdb')
# mol_1=create_3d_from_xyz('conf_pyt.xyz','1.pdb')
# mol_2=create_3d_from_xyz('conf_bash.xyz','2.pdb')

