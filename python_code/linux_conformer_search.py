# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 14:45:24 2022

@author: edens
"""
from rdkit.Chem import PandasTools
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel
from openbabel import openbabel as ob
import os

"""
rdkit conformer search-
"""
def xyz_to_sdf(file_name,output_type='sdf'):
    obConversion = ob.OBConversion()
    file_type=file_name.split('.')[-1]
    obConversion.SetInAndOutFormats(file_type, output_type)
    obmol = ob.OBMol()
    obConversion.ReadFile(obmol, file_name)
    pbmol=pybel.Molecule(obmol)
    pbmol.write(output_type,(file_name.split('.')[0]+'.sdf'))
    return 

def rdkit_conformer_write(mol, info): ###not in use at the moment
    writer = Chem.SDWriter('conformers.sdf')
    for cid,energy in zip(range(mol.GetNumConformers()),info):
        mol.SetProp('ID',f'energy_{energy}')
        writer.write(mol, confId=cid)
    writer.close()

def rdkit_conformers_search(file_name, forcefield='UFF',add_ref=True):
    xyz_to_sdf(file_name)
    mol =Chem.MolFromMolFile((file_name.split('.')[0]+'.sdf'),removeHs=False)
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
        
       
    w = Chem.SDWriter(('rdkitconformers_'+(file_name.split('.')[0]+'.sdf')))
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
    os.remove(os.path.abspath((file_name.split('.')[0]+'.sdf')))
    return 
    
if __name__ == '__main__':
    xyz_files=[filename for filename in os.listdir() if 'xyz' in filename]
    os.chdir(r'/gpfs0/gaus/users/edenspec/crest_runs/test')
    for xyz_file in xyz_files:
        molecule_name=xyz_file.split('.')[0]
        
        
        xyz_filename=xyz_file
        os.system('/gpfs0/gaus/users/edenspec/working_obabel/build/bin/obabel {} -O {} --conformer --nconf 30--writeconformers'.format(*xyz_file,('obconformers_'+xyz_file)))
        os.system('/gpfs0/gaus/projects/crest --input {} -gfn2 -g h2o -T 4 --xnam /gpfs0/gaus/projects/xtb-6.4.1/bin/xtb'.format(xyz_file))
        rdkit_conformers_search(xyz_filename)
        # os.mkdir(path)
        # os.chdir(path)
        # os.system('cp '+Constant.HOME_DIR.value+'/fixed_locations.inp '+Constant.HOME_DIR.value+'/'+molecule_name+'/fixed_locations.inp')
        # os.system('cp '+Constant.HOME_DIR.value+'/'+smile_filename+' '+Constant.HOME_DIR.value+'/'+molecule_name+'/'+smile_filename)
        # os.system(CrestConstants.OBABEL.value+smile_filename+CrestConstants.OBABEL_XYZ_SETTINGS_1.value+xyz_filename+CrestConstants.OBABEL_XYZ_SETTINGS_2.value)
        # os.system('echo >> '+xyz_filename)
        # os.system('echo >> '+xyz_filename)
        # os.system('echo "\$write" >> '+xyz_filename)
        # os.system('echo "   output file=properties.out" >> '+xyz_filename)
        # os.chdir(Constant.HOME_DIR.value)
        # os.system(CrestConstants.REMOVEL_COMMAND.value+smile_filename)
        # xyz_filenames.append(xyz_filename)