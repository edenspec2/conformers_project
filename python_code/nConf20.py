# Calculate nConf20 - from https://pubs.acs.org/doi/10.1021/acs.jcim.6b00565
# Jerome G. P. Wicker and Richard I. Cooper  J. Chem. Inf. Model. 2016 56(12) pp. 2347-2352
# reformatted for Python3 by Geoffrey R. Hutchison https://hutchison.chem.pitt.edu 

from collections import OrderedDict
import fileinput
import sys
print(sys.executable)
print(sys.path)
import numpy as np
import openbabel as ob



from rdkit import Chem
from rdkit.Chem import AllChem

def create_3d_from_smiles(smiles_code): #'O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N'
    obmol = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInFormat('smi')
    conv.ReadString(obmol, smiles_code)
    return obmol

def GenerateConformers(mol, nConfs):
    # Add H atoms to skeleton
    molecule = Chem.AddHs(mol)
    conformerIntegers = []
    # Embed and optimise the conformers with MMFF94
    conformers = AllChem.EmbedMultipleConfs(molecule, nConfs, pruneRmsThresh=0.5)
    optimized_and_energies = AllChem.MMFFOptimizeMoleculeConfs(molecule, maxIters=600, nonBondedThresh=100.0)

    EnergyDictionaryWithIDAsKey = {}
    FinalConformersToUse = {}

    # Only keep the conformers which were successfully fully optimised
    for conformer in conformers:
        optimised , energy = optimized_and_energies[conformer]
        if optimised == 0:
            EnergyDictionaryWithIDAsKey[conformer] = energy
            conformerIntegers.append(conformer)

    # Keep the lowest energy conformer
    lowestenergy = min( EnergyDictionaryWithIDAsKey.values() )

    for k, v in EnergyDictionaryWithIDAsKey.items():
        if v == lowestenergy:
            lowestEnergyConformerID = k

    FinalConformersToUse [ lowestEnergyConformerID ] = lowestenergy

    # Remove H atoms to speed up substructure matching
    molecule = AllChem.RemoveHs(molecule)
    # Find all substructure matches of the molecule with itself , to account for symmetry
    matches = molecule.GetSubstructMatches(molecule, uniquify=False)
    maps = [ list (enumerate(match)) for match in matches]

    # Loop over conformers other than the lowest energy one
    for conformerID in EnergyDictionaryWithIDAsKey.keys() :
        okayToAdd = True

        # Loop over reference conformers already added to list
        for finalconformerID in FinalConformersToUse.keys() :

            # Calculate the best RMS of this conformer with the reference conformer in the list
            RMS = AllChem.GetBestRMS(molecule, molecule , finalconformerID , conformerID , maps)
            # Do not add if a match is found with a reference conformer
            # (i.e., heavy-atom RMS within 1.0 Å)
            if RMS< 1.0:
                okayToAdd = False
                break

        # Add the conformer if the RMS is greater than 1.0 for every reference conformer
        if okayToAdd :
            FinalConformersToUse[conformerID] = EnergyDictionaryWithIDAsKey [conformerID]

    # Sort the conformers by energy
    sortedDictionary = OrderedDict (sorted( FinalConformersToUse.items(), key=lambda t: t[1]))

    energies = [val for val in sortedDictionary.values()]
    return energies

def CalcNConf20(energylist, upper=20): 
    # Count up the number of conformers under 20 kcal/mol
    descriptor = 0
    # Subtract the lowest energy found in the ordered list
    relativeEnergies = np.array(energylist) - energylist[0]

    # Only look at the energies of conformers other than the global minimum
    for energy in relativeEnergies [1:]:
        # Optimized lower and upper energy limits for conformer energy
        if 0 <= energy < upper:
            descriptor += 1
    return descriptor

if __name__ == "__main__":
    # Assume we're reading a bunch of SMILES from the stdin
    smiles = ('O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N')
    molecule = Chem.MolFromSmiles(smiles)
    energyList = GenerateConformers(molecule, 50)
    print(CalcNConf20(energyList))
    
    
    # for line in fileinput.input():
    #     smiles = line.split()[0]
    #     molecule = Chem.MolFromSmiles(smiles)
    #     # maximum of 50 conformers was the limit in the paper
    #     energyList = GenerateConformers(molecule, 50)
    #     print(CalcNConf20(energyList))