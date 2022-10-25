# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 15:19:49 2022

@author: edens
"""
import numpy as np
from openbabel import pybel
from openbabel import openbabel as ob
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd 
import help_functions
from rdkit.Chem import PandasTools
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D
import rdkit_utils
path=r'C:\Users\edens\Documents\GitHub\conformers_project\python_code\test_run'
os.chdir(path)

"""
general functions
"""
def split_molecule_file(file_name): ##for crest output
    file_type=file_name.split('.')[-1]
    molecules=[pbmol.write(file_type) for pbmol in pybel.readfile(file_type,file_name)]
    return molecules

def xyz_to_sdf(file_name,output_type='sdf'):
    obConversion = ob.OBConversion()
    file_type=file_name.split('.')[-1]
    obConversion.SetInAndOutFormats(file_type, output_type)
    obmol = ob.OBMol()
    obConversion.ReadFile(obmol, file_name)
    pbmol=pybel.Molecule(obmol)
    pbmol.write(output_type,help_functions.change_filetype(file_name,'sdf'))
    return 

# def mol_to_coordinates_array(mol):   
#   editable_mol = Chem.RWMol(mol)
#   conformer = Chem.Conformer(l
#   for i, s in enumerate(symbols):
#       atom = editable_mol.AddAtom(Chem.Atom(_atomic_number[s]))
#       atom_position = Point3D(geometry[i][0], geometry[i][1], geometry[i][2])
#       conformer.SetAtomPosition(atom, atom_position)
"""
openbabel conformers search -
"""
def obmol_from_coordinates(file_name,output_type='sdf'):
    obConversion = ob.OBConversion()
    file_type=file_name.split('.')[-1]
    obConversion.SetInAndOutFormats(file_type, output_type)
    obmol = ob.OBMol()
    obConversion.ReadFile(obmol, file_name)
    return obmol

def freeze_atoms_for_confab(obmol,atoms_to_freeze):
    constraints = ob.OBFFConstraints()
    for atom in ob.OBMolAtomIter(obmol):
        atom_id = atom.GetIndex() 
        if atom_id in atoms_to_freeze:
            constraints.AddAtomConstraint(atom_id)
    return constraints

def confab_search(obmol,file_name,set_constraints=False,atoms_to_freeze=None,output_format='sdf'):
    
    pff = ob.OBForceField_FindType( "mmff94" )
    if set_constraints :
        constraints=freeze_atoms_for_confab(obmol,atoms_to_freeze)
        pff.SetConstraints(constraints)
    pff.DiverseConfGen(0.5, 1000, 50.0, True) #allow change
    pff.SetLogLevel(ob.OBFF_LOGLVL_HIGH)
    pff.Setup(obmol)
    pff.GetConformers(obmol)
    obconversion = ob.OBConversion()
    obconversion.SetOutFormat(output_format)
    output_strings = []
    for conf_num in range(obmol.NumConformers()):
        obmol.SetConformer(conf_num)
        output_strings.append(obconversion.WriteString(obmol))
        
    # with open(('obconformers_'+help_functions.change_filetype(file_name,'sdf')),'a') as f:
    #     for conf_num in range(obmol.NumConformers()):
    #         obmol.SetConformer(conf_num)
    #         # output_strings.append(obconversion.WriteString(obmol))
    #         f.write(obconversion.WriteString(obmol))
    return output_strings


"""
rdkit conformer search-
"""
def rdkit_conformer_write(mol, info): ###not in use at the moment
    writer = Chem.SDWriter('conformers.sdf')
    for cid,energy in zip(range(mol.GetNumConformers()),info):
        mol.SetProp('ID',f'energy_{energy}')
        writer.write(mol, confId=cid)
    writer.close()

def rdkit_conformers_search(file_name, forcefield='UFF',add_ref=True):
    mol =Chem.MolFromMolFile(file_name,removeHs=False)
    refmol = Chem.Mol(mol)
    # Number of conformers to be generated
    num_of_conformer=100
    max_iter=500
    # Default values for min energy conformer
    min_energy=10000
    min_energy_index=0
    cids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_of_conformer,params=rdDistGeom.ETKDGv2())
    ids = list(cids)
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    # print(mol.GetConformer(-1).GetPositions())
    
    if forcefield=='UFF':
        results= AllChem.UFFOptimizeMoleculeConfs(mol,maxIters=max_iter)
    else:
        results = AllChem.MMFFOptimizeMoleculeConfs(mol,maxIters=max_iter,mmffVariant='MMFF94s')
        
       
    w = Chem.SDWriter(('rdkitconformers_'+help_functions.change_filetype(file_name,'sdf')))
    if add_ref:
        mp_origin = AllChem.MMFFGetMoleculeProperties(refmol, mmffVariant='MMFF94s')
        ff = AllChem.MMFFGetMoleculeForceField(refmol,mp_origin)
        e = ff.CalcEnergy()
        refmol.SetProp('CID', '-1')
        refmol.SetProp('Energy', str(e))
        w.write(refmol)
    res = []
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
        e = ff.CalcEnergy()
        res.append((cid, e))
    sorted_res = sorted(res, key=lambda x:x[1])
    # print(sorted_res)
    rdMolAlign.AlignMolConformers(mol)
    for cid, e in sorted_res:
        mol.SetProp('CID', str(cid))
        mol.SetProp('Energy', str(e))
        w.write(mol, confId=cid)
    w.close()
    
    for i in range(mol.GetNumConformers()):
        refmol.AddConformer(mol.GetConformer(i))
    
    # rdkit_conformer_write(mol,results)
    
    return mol,refmol

def create_conformers_from_dir(path):
    os.chdir(path)
    xyz_files=help_functions.get_file_name_list('xyz')
    # help_functions.create_molecule_directories()
    try:
        [xyz_to_sdf(file) for file in xyz_files]
    except OSError:
        pass
    sdf_files=help_functions.get_file_name_list('sdf')
    
    mol_num=0
    df=pd.DataFrame()
    for xyz_file,sdf_file in zip(xyz_files,sdf_files):
        obmol=obmol_from_coordinates(xyz_file)
        # x.append(confab_search(obmol,xyz_file))
        mol_num+=1
        rdmol=rdkit_conformers_search(sdf_file)
        df['rdkitconformer_{}'.format(mol_num)]=split_molecule_file(('rdkitconformers_'+sdf_file))
        print(sdf_file)
        # help_functions.move_files_directory('conformers')
    
    # help_functions.delete_type_files('sdf')
    
    return df
    
if __name__ == '__main__':        
    
    x=create_conformers_from_dir(path)
 
    # obmol=obmol_from_coordinates('origin_molecule.xyz')
    # rdmol=rdkit_conformers_search('origin_molecule.sdf')[0]
    # rdmol_2=rdkit_conformers_search('origin_molecule.sdf')[1]
    # conformer=rdkit_utils.ConformerGenerator()
    # energy_list=conformer.get_conformer_energies(rdmol_2)
    
    # # confab_search(obmol,'origin_molecule.xyz')
    # rdmol_2.GetConformer(0).GetPositions()
    # obconformer_list=split_molecule_file('obconformers_origin_molecule.sdf')
    # rdkitconformer_list=split_molecule_file('rdkitconformers_origin_molecule.sdf')
    # x=Chem.SupplierFromFilename('obconformers_origin_molecule.sdf')
    
  
    # molecules=Chem.SupplierFromFilename('obconformers_origin_molecule.sdf')
    # df = PandasTools.LoadSDF('obconformers_origin_molecule.sdf')
