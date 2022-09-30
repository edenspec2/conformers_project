import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import os
import math
from enum import Enum
import help_functions 
import statistics
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

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
        os.chdir(self.molecule_path)
        self.list_of_files=[file_name for file_name in os.listdir(self.molecule_path)]
        self.origin_coordinates_df=self.get_unedited_df('coor')
        os.chdir('../')
     
    def get_file_name(self,file_type):
        return [file_name for file_name in self.list_of_files if file_type in file_name][0]



if __name__=='__main__':
    molecule_1=Molecule('try_molecule')
    # fig = plt.figure()
    x=molecule_1.get_unedited_df('coor').drop([0],axis=1)
  