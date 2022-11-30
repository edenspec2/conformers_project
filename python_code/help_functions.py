import pandas as pd
import numpy as np
import os
import glob

from enum import Enum


class FileExtensions(Enum):
    """
    Hold commonly used file extensions
    """
    SMI='.smi'
    XYZ='xyz'
    CSV='.csv'
    ZIP='.zip'
    PPT='.ppt'
    CIF='.cif'
    MOL='.mol'
    PDB='.pdb'

class XYZConstants(Enum):
    """
    Constants related to XYZ file processing
    """
    DF_COLUMNS=[['atom','x','y','z']]
    
    

def get_file_name_list(file_identifier):
    """
    The function gets a file identifier as input and returns a list of all files in the working 
    which contain the identifier in the files name
    ----------
    Parameters
    ----------
    identifier : str.
        The wanted file identifier like 'txt','info','nbo' contained in the filename
    -------
    Returns
    -------
    list
        A list of all files in the working directory with the chosen extension 
    --------
    Examples
    --------
    
    all_files_in_dir=listdir()
    print(all_files_in_dir)
        ['0_1106253-mod-mod.xyz', '0_1106253-mod.xyz', '1106253.cif', '1109098.cif', '1_1106253-mod.xyz', 'centered_0_BASCIH.xyz', 'cif_handler.py']
        
    xyz_files_in_dir=get_filename_list('.xyz')
    print(xyz_files_in_dir)
        ['0_1106253-mod-mod.xyz', '0_1106253-mod.xyz', '1_1106253-mod.xyz', 'centered_0_BASCIH.xyz']
  
    """
    return [filename for filename in os.listdir() if file_identifier in filename]

def get_df_from_csv(filename,columns=None,index=None,splitter=None):
    """
    Parameters
    ----------
    filename : str
        full file name to read.
    columns : str , optional
        list of column names for DataFrame. The default is None.
    splitter : str, optional
        input for [.split().] , for csv-',' for txt leave empty. The default is None.
    dtype : type, optional
        type of variables for dataframe. The default is None.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    with open(filename, 'r') as f:
        lines=f.readlines()
    splitted_lines=[line.split(splitter) for line in lines]
    df=pd.DataFrame(splitted_lines,columns=columns,index=index)
    return df

def convert_file_to_xyz_df(filename,splitter=','):
    df=get_df_from_csv(filename,columns=XYZConstants.DF_COLUMNS.value,splitter=splitter)
    # df['atom'].replace(main.GeneralConstants.ATOMIC_NUMBERS.value, inplace=True)
    return df

def dataframe_to_xyz(dataframe, output_name, comment_line=''):
    """

     a function that recieves a dataframe, output name, and comment line and creates a xyz type file.
     
    parameters

    ---

    dataframe: an array that can contain different classes of data, needs to be 4 colums to run.

    output_name:str, the name for the file created.

    comment_line: str, the headline of the file .
    ---

    examples:
    ---
    """
    number_of_atoms=dataframe.shape[0]
    atoms_np_array=dataframe.to_numpy()
    with open(output_name, 'w') as xyz_file:
        xyz_file.write("{}\n{}\n".format(number_of_atoms, comment_line))
        for atom_np_array in atoms_np_array:
            try:
                xyz_file.write("{:1} {:11.20} {:11.20} {:11.20}".format(*atom_np_array))
            except:
                xyz_file.write("{:1}".format(*atom_np_array))

                               
def change_filetype (filename,new_type='xyz'):
    """
    a function that recieves a file name, and a new type, and changes the type-ending of the file's name to the new one.

    parameters
    ---
    filename: str, the file we want to change

    new_type:str, the new ending we want for the file

    returns
    ---
    the same file name with a new type ending

    examples
    ---
    filename='B_THR_127_5Angs_noHOH.pdb'
    new_filename=change_filetype(filename,'xyz')
    OUTPUT:'B_THR_127_5Angs_noHOH.xyz'
    
    """
    split_result=filename.split('.')
    if '.' in new_type:
        new_filename=split_result[0]+new_type
    else:
        new_filename=split_result[0]+'.'+new_type
    return new_filename

def xyz_string_to_df(lines):

    strip_lines=[line.strip().rstrip('\n') for line in lines]
    
    
    
    
def create_molecule_directories():
    list_of_dirs=[name.split('.')[0] for name in os.listdir()]
    for dir_name in list_of_dirs:
        os.mkdir(dir_name)
    return
def delete_type_files(file_type='xyz'): ## my help function to delete xyz files
    """
    a function that gets a directory path and file type, and deletes all the files of said type.
    """
    list_of_molecules=[file for file in os.listdir() if file.endswith(file_type)]
    for molecule in list_of_molecules:
        os.remove(os.path.abspath(molecule))
        
        
def move_files_directory(file_type):#need edit
    """
    a function that moves xyz type files from one directory to another.
    help function for xyz_file_generator_library to move files to the new directory created
    A void function
    """
    list_of_dirs=[name for name in os.listdir() if os.path.isdir(os.path.abspath(name))]
    list_of_files=get_file_name_list(file_type)
    for file_name,dir_name in zip(list_of_files,list_of_dirs):   
        new_path=os.path.join(os.path.abspath(dir_name),file_name)
        os.replace(os.path.abspath(file_name),new_path)
    return
## write help functions for mol2 to dfs
