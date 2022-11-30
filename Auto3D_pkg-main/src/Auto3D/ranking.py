#!/usr/bin/env python
'''
Finding 3D structures that satisfy the input requirement.
'''
import os
import logging
import pandas as pd
from openbabel import pybel
from .utils import guess_file_type, filter_unique
from .utils import hartree2ev, ev2kcalpermol

# eV2kcalmol = 23.06055  #1 eV = 23.06055 kcal/mol
# hartree2eV = 27.211385
class ranking(object):
    '''
    Finding 3D structures that satisfy the user-defined requirements.

    Arguments:
        input_path: An SDF file that contains all isomers.
        energies: A txt file that contains the IDs and energies.
        out_path: An SDF file that stores the optimical structures.
        k: Outputs the top-k structures for each SMILES.
        window: Outputs the structures whose energies are within
                window (kcal/mol) from the lowest energy
    Returns:
        None
    '''
    def __init__(self, input_path,
                 out_path, threshold, k=False, window=False):
        self.input_path = input_path
        self.out_path = out_path
        self.threshold = threshold
        self.atomic_number2symbol = {1: 'H', 
                                     5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 
                                     14: 'Si', 15: 'P', 16: 'S', 17: 'Cl',
                                     32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br',
                                     51: 'Sb', 52: 'Te', 53: 'I'}
        self.k = k
        self.window = window

    @staticmethod
    def similar(name, names):
        name2 = name.strip().split("_")[0]
        names2 = names[0].strip().split('_')[0]
        return (name2 == names2)

    @staticmethod
    def add_relative_e(list0):
        """Adding relative energies compared with lowest-energy structure
        
        Argument:
            list: a list of tuple (idx, name, energy)
            
        Return:
            list of tuple (idx, name, energy, relative_energy)
        """
        list0_ = []
        _, _, e_m = list0[0]
        for idx_name_e in list0:
            idx_i, name_i, e_i = idx_name_e
            e_relative = e_i - e_m
            list0_.append((idx_i, name_i, e_i, e_relative))
        return list0_

    def top_k(self, energies, names, mols, k=1):
        '''
        Given a group of energy_name_idxes,
        return the top-k lowest name-energies pairs with idxes as keys.
        '''
        assert(len(energies) == len(names))
        assert(len(energies) == len(mols))
        assert(len(set(names)) == 1)

        df = pd.DataFrame({"names": names, "energies": energies, "mols": mols})
        df2 = df.sort_values(by=['energies'])
        
        out_mols_ = filter_unique(list(df2["mols"]), self.threshold)
        if k < len(out_mols_):
            out_mols = out_mols_[:k]
        else:
            out_mols = out_mols_

        if len(out_mols) == 0:
            name = names[0].split("_")[0].strip()
            print(f"No structure converged for {name}, try a larger convergence threshold.")
            logging.info(f"No structure converged for {name}, try a larger convergence threshold.")
        else:
            #Adding relative energies
            ref_energy = out_mols[0].data['E_tot']
            for mol in out_mols:
                my_energy = mol.data['E_tot']
                rel_energy = float(my_energy) - float(ref_energy)
                mol.data['E_rel(eV)'] = rel_energy
        return out_mols

    def top_window(self, energies, names, mols, window=1):
        '''
        Given a group of energy_name_idxes,
        return all (idx, name, e) tuples whose energies are within
        window (Hatree) from the lowest energy. Unit table is based on: 
        http://wild.life.nctu.edu.tw/class/common/energy-unit-conv-table.html
        '''
        window = (window/ev2kcalpermol)  # convert energy window into eV unit
        assert(len(energies) == len(names))
        assert(len(energies) == len(mols))
        assert(window >= 0)
        assert(len(set(names)) == 1)


        df = pd.DataFrame({"names": names, "energies": energies, "mols": mols})
        df2 = df.sort_values(by=['energies'])

        out_mols_ = filter_unique(list(df2['mols']), self.threshold)
        out_mols = []

        if len(out_mols_) == 0:
            name = names[0].split("_")[0].strip()
            print(f"No structure converged for {name}.")
            logging.info(f"No structure converged for {name}.")
        else:
            ref_energy = float(out_mols_[0].data['E_tot'])
            for mol in out_mols_:
                my_energy = float(mol.data['E_tot'])
                rel_energy = my_energy - ref_energy
                if rel_energy <= window:
                    mol.data['E_rel(eV)'] = rel_energy
                    out_mols.append(mol)
                else:
                    break
        return out_mols

    def run(self):
        """
        When runs, lowest-energy structure will be stored in out_path folder.
        """
        print("Beggin to slelect structures that satisfy the requirements...")
        logging.info("Beggin to slelect structures that satisfy the requirements...")
        results = []

        file_type = guess_file_type(self.input_path)
        data2 = pybel.readfile(file_type, self.input_path)

        mols = [mol for mol in data2 if mol.data["Converged"].lower() == "true"]
        names = [mol.title.strip() for mol in mols]
        energies = [float(mol.data['E_tot']) for mol in mols]
        assert(len(mols) == len(names))
        assert(len(mols) == len(energies))

        #Grouping, ranking
        names2 = map(lambda x: x.strip().split("_")[0].strip(), names)
        df = pd.DataFrame({"names": names2, "energies": energies, "mols": mols})
        df2 = df.groupby("names")
        for group_name in df2.indices:
            group = df2.get_group(group_name)
            # Filter similar structures
            g_mols = list(group.loc[:, 'mols'])
            g_names = list(group.loc[:, 'names'])
            g_energies = list(group.loc[:, 'energies'])

            if self.k:
                top_results = self.top_k(g_energies, g_names,
                                            g_mols, self.k)
            elif self.window:
                top_results = self.top_window(g_energies, g_names,
                                                g_mols, self.window)
            else:
                raise ValueError(('Parameter k or window needs to be '
                                    'specified. Append "--k=1" if you'
                                    'only want one structure per SMILES'))
            results += top_results

        f = pybel.Outputfile('sdf', self.out_path)
        for mol in results:
            # Change the energy unit from eV back to Hartree
            mol.data['E_tot'] = (float(mol.data['E_tot'])/hartree2ev)
            mol.data['E_rel(kcal/mol)'] = (float(mol.data['E_rel(eV)']) * ev2kcalpermol)
            del mol.data['E_rel(eV)']
            #Remove _ in the molecule title
            t = mol.title
            t_simplified = t.split("_")[0].strip()
            mol.title = t_simplified
            f.write(mol)
        f.close()
