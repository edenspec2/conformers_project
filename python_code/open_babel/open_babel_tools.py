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
def pbmol_to_3d(pbmol,out_put_name):
    pbmol.write('pdb',out_put_name,overwrite=True)
    return 

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
    return pbmol_to_3d(mol,out_put_name)

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

def smile_to_coordinates(smile):
    pbmol = pybel.readstring('smi', smile)
    pbmol.make3D(forcefield="gaff", steps=100)
    pbmol.localopt()
    return pbmol.write('pdb')


def smiles_to_coordinates(smiles):
    return [smile_to_coordinates(smile) for smile in smiles]