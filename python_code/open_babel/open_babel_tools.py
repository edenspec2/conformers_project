from openbabel import pybel
import openbabel as ob
import os
path_to_add=r'C:\Users\edens\Documents\GitHub'
os.chdir(path_to_add)

def read_molecules(filepath, single = False):
    in_format = filepath.strip().split( '.' )[-1]
    obconversion = ob.OBConversion()
    obconversion.SetInFormat( in_format )
    obmol = ob.OBMol()

    molecules = []

    notatend = obconversion.ReadFile( obmol, filepath )
    while notatend:
        molecules.append( obmol )
        obmol = ob.OBMol()
        notatend = obconversion.Read( obmol )

    if single:
        assert( len(molecules) == 1 )
        return molecules[0]
    else:
        return molecules
def molecule_to_3d(mol,out_put_name):
    pbmol = pybel.Molecule(mol)
    pbmol.addh()
    pbmol.make3D(forcefield="gaff",steps=100)
    pbmol.localopt()
    pbmol.write('pdb',out_put_name,overwrite=True)
    return pbmol.OBMol

def molecule_from_coordinates(file_type,file_name):
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats(file_type, "mol2")
    mol = ob.OBMol()
    obConversion.ReadFile(mol, file_name)
    return mol

def create_3d_from_smiles(smiles_code,out_put_name):
    mol = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInFormat('smi')
    conv.ReadString(mol, smiles_code)
    return molecule_to_3d(mol,out_put_name)

def create_3d_from_xyz(xyz_file_name,out_put_name):
    mol = molecule_from_coordinates('xyz',xyz_file_name)
    return molecule_to_3d(mol,out_put_name)

def confab_search(mol,output_format):
    pff = ob.OBForceField_FindType( "mmff94" )
    pff.DiverseConfGen(0.5, 100, 50.0, True) #allow change
    pff.GetConformers(mol)
    obconversion = ob.OBConversion()
    obconversion.SetOutFormat(output_format)
    output_strings = []
    for conf_num in range(mol.NumConformers()):
        mol.SetConformer(conf_num)
        output_strings.append(obconversion.WriteString(mol))
    return output_strings
                
def get_molecule_charge(mol):
    ob_charge_model = ob.OBChargeModel.FindType("eem2015bn")
    ob_charge_model.ComputeCharges(mol)
    return ob_charge_model.GetPartialCharges()

def get_molecule_dipole(mol):
    ob_charge_model = ob.OBChargeModel.FindType("eem2015bn")
    return ob_charge_model.GetDipoleMoment(mol) ##cant open vector3 object

class Molecules():
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
        os.chdir(self.molecules_path)
        self.molecule=ob.OBMol()
        self.file_type=os.listdir(self.molecule_path)[0]
       
        self.list_of_molecules=[file_name for file_name in os.listdir(self.molecule_path)]
        self.origin_coordinates_df=self.get_unedited_df('coor')
        os.chdir('../')
     
    def get_file_name(self,file_type):
        return [file_name for file_name in self.list_of_files if file_type in file_name][0]



# mol = create_3d_from_smiles(('O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N'),'3d.pdb')
# mol_1=create_3d_from_xyz('conf_pyt.xyz','1.pdb')
# mol_2=create_3d_from_xyz('conf_bash.xyz','2.pdb')

##writes xyz file with conformers
mol=molecule_from_coordinates('xyz','coor.xyz')
data=ob.OBGenericData()
print(data.GetValue())
conformers=confab_search(mol,'xyz')
charges = get_molecule_charge(mol)
dipole=get_molecule_dipole(mol)

# obConversion = ob.OBConversion()
# obConversion.SetInAndOutFormats("smi", "mdl")

# mol = ob.OBMol()
# obConversion.ReadString(mol, "C1=CC=CS1")

# print ('Should print 5 (atoms)')
# print (mol.NumAtoms())

# mol.AddHydrogens()
# print ('Should print 9 (atoms) after adding hydrogens')
# print (mol.NumAtoms())

# outMDL = obConversion.WriteString(mol)

# ff = ob.OBForceField.FindForceField("mmff94")
# ff.Setup(mol)
