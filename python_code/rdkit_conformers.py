"""
Created on Thu Feb  6 16:14:17 2020

@author: Murat Cihan Sorkun

A basic search for the minimum energy conformer
Note: It is more efficient to align similar conformers first, but this example does not include it
"""
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd 

path=r'C:\Users\edens\Documents\GitHub\conformers_project\python_code\test_run'
os.chdir(path)
# Creating mol object from SMILES and add H's
mol = Chem.MolFromSmiles("O=C(NC1CCCCC1)NS(=O)(c1ccc(CCNC(c2cnc(cn2)C)=O)cc1)=O")
mol_h_UFF = Chem.AddHs(mol)
mol_h_MMFF = Chem.AddHs(mol)


def ConfToMol(mol, conf_id):
    conf = mol.GetConformer(conf_id)
    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(Chem.Conformer(conf))
    return new_mol

def file_to_mol(file_name):
    mol=Chem.MolFromMolFile(file_name,removeHs=False) ##sdf file
    return Chem.AddHs(mol)



def GetConformerRMSdf(mol, atomIds=None, prealigned=False): ##edited to return a df with all possible rms for each conformer
    """ Returns the RMS matrix of the conformers of a molecule.
    As a side-effect, the conformers will be aligned to the first
    conformer (i.e. the reference) and will left in the aligned state.

    Arguments:
      - mol:     the molecule
      - atomIds: (optional) list of atom ids to use a points for
                 alingment - defaults to all atoms
      - prealigned: (optional) by default the conformers are assumed
                    be unaligned and will therefore be aligned to the
                    first conformer

    Note that the returned RMS matrix is symmetrical, i.e. it is the
    lower half of the matrix, e.g. for 5 conformers::

      rmsmatrix = [ a,
                    b, c,
                    d, e, f,
                    g, h, i, j]

    where a is the RMS between conformers 0 and 1, b is the RMS between
    conformers 0 and 2, etc.
    This way it can be directly used as distance matrix in e.g. Butina
    clustering.

    """
    # if necessary, align the conformers
    # Note: the reference conformer is always the first one
    rmsvals = []
    df = pd.DataFrame([])
    confIds = [conf.GetId() for conf in mol.GetConformers()]
    if not prealigned:
        if atomIds:
           AllChem.AlignMolConformers(mol, atomIds=atomIds, RMSlist=rmsvals)
        else:
            AllChem.AlignMolConformers(mol, RMSlist=rmsvals)
    else:  # already prealigned
        for i in range(1, len(confIds)):
            rmsvals.append(
              AllChem.GetConformerRMS(mol, confIds[0], confIds[i], atomIds=atomIds, prealigned=prealigned))
    
    df['base_conformer']=pd.DataFrame(rmsvals)
    for i in range(1, len(confIds)):
        cmat=[]
        for j in range(0, len(confIds)):
            if i!=j:
                cmat.append(AllChem.GetConformerRMS(mol, confIds[i],
                            confIds[j], atomIds=atomIds, prealigned=True))
        df['conformer_{}'.format(i)]=pd.DataFrame(cmat)
    return df
if __name__ == '__main__':
    # conformers=rdkit_conformers_search(mol)
    
    # volume=[AllChem.ComputeMolVolume(conformers,confId=cid) for cid in range(conformers.GetNumConformers())]
    
    mol2=Chem.MolFromMolFile('origin_molecule.sdf',removeHs=False) 
    # conformers=rdkit_conformers_search(mol2)
    # rms=AllChem.GetConformerRMSdf(conformers)
    
    
    
    
    # # Write minimum energy conformers into a SDF file
    # w = Chem.SDWriter('minimum-energy-conformer-UFF.sdf')
    # w.write(Chem.Mol(mol_h_UFF,False,min_energy_index_UFF))
    # w.flush()  
    # w.close()
    
    
    # # Write minimum energy conformer into a SDF file
    # w = Chem.SDWriter('minimum-energy-conformer-MMFF.sdf')
    # w.write(Chem.Mol(mol_h_MMFF,False,min_energy_index_MMFF))
    # w.write()
    # w.flush()  
    # w.close()
