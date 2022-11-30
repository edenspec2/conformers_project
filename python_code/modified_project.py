
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import os
import math
from enum import Enum
import help_functions 
import statistics
import networkx as nx
import igraph as ig
import matplotlib.pyplot as plt

pd.options.plotting.backend
coordinates_columns=[['x','y','z']]
dipole_columns=[['dip_x','dip_y','dip_z','total']]

class GeneralConstants(Enum):
    """
    Holds constants for calculations and conversions
    1. covalent radii from Alvarez (2008) DOI: 10.1039/b801115j
    2. atomic numbers
    2. atomic weights
    """
    COVALENT_RADII= {
            'H': 0.31, 'He': 0.28, 'Li': 1.28,
            'Be': 0.96, 'B': 0.84, 'C': 0.76, 
            'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
            'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
            'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
            'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
            'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
            'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
            'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
            'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
            'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
            'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
            'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
            'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
            'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
            'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
            'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
            'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
            'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
            'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
            'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
            'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
            'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
            'Am': 1.80, 'Cm': 1.69
    }
    
    CPK_RADII={
        'H': 1.10, 'C': 1.70, 'F': 1.47,
        'S': 1.80, 'B': 1.92, 'I': 1.98, 
        'N': 1.55, 'O': 1.52, 'Co': 2.00, 
        'Br': 1.83, 'Si': 2.10,'Ni': 2.00,
        'P': 1.80, 'Cl': 1.75, 
    }
    
    BONDI_RADII={
        'C':1.50,   'H':1.00,   'S.O':1.70,  'Si':2.10,
        'C2':1.60,  'N':1.50,   'S1':1.00,   'Co':2.00,
        'C3':1.60,  'C66':1.70, 'F':1.35,    'Ni':2.00,
        'C4':1.50,  'N4':1.45,  'Cl':1.80,
        'C5/N5':1.70, 'O':1.35, 'S4':1.40,
        'C6/N6':1.70, 'O2':1.35, 'Br':1.95,
        'C7':1.70,    'P':1.40,  'I':2.15,
        'C8':1.50,    'S':1.70,  'B':1.92,
    
    }
    
    BOND_TYPE={
        
        'O.2':'O2', 'N.2':'C6/N6','S.3':'S4',
        'O.3':'O', 'N.1':'N', 'S.O2':'S1',
        'O.co2':'O', 'N.3':'C6/N6','P.3':'P',
        'C.1':'C3', 'N.ar':'C6/N6',
        'C.2':'C2', 'N.am':'C6/N6',
        "C.cat":'C3', 'N.pl3':'C6/N6',
        'C.3':'C', 'N.4':'N',
        'C.ar':'C6/N6', 'S.2':'S',
        }
    REGULAR_BOND_TYPE={
        
        'O.2':'O', 'N.2':'N','S.3':'S',
        'O.3':'O', 'N.1':'N', 'S.O2':'S',
        'O.co2':'O', 'N.3':'N','P.3':'P',
        'C.1':'C', 'N.ar':'N',
        'C.2':'C', 'N.am':'N',
        "C.cat":'C', 'N.pl3':'N',
        'C.3':'C', 'N.4':'N',
        'C.ar':'C', 'S.2':'S',
        }
    
    ATOMIC_NUMBERS ={
    '1':'H', '5':'B', '6':'C', '7':'N', '8':'O', '9':'F', '14':'Si',
             '15':'P', '16':'S', '17':'Cl', '35':'Br', '53':'I', '27':'Co', '28':'Ni'}
        

    ATOMIC_WEIGHTS = {
            'H' : 1.008,'He' : 4.003, 'Li' : 6.941, 'Be' : 9.012,
            'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,
            'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,
            'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,
            'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,
            'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,
            'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,
            'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,
            'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,
            'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,
            'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,
            'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,
            'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,
            'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,
            'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,
            'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'Gd' : 157.25,
            'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,
            'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,
            'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,
            'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,
            'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,
            'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,
            'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,
            'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247
    }


def delete_type_files(file_type='xyz'): ## my help function to delete xyz files
    """
    a function that gets a directory path and file type, and deletes all the files of said type.
    #works in the current directory
    """
    list_of_molecules=[file for file in os.listdir() if file.endswith(file_type)]
    for molecule in list_of_molecules:
        os.remove(os.path.abspath(molecule))
        

def xyz_file_generator(folder_path):###works, name in R:'xyz_file_generator'-currently the R function generates file_name without 'xyz_'.
    """
    a function that gets a directory path as folder_path, makes xyz files from all csv files in the same directory.
    This is a void function
    works in a specific molecule directory.
    example:
        path=GitHub\learning_python\project\main_python\test_dipole\molecule1
        xyz_file_generator(path)
    """
    os.chdir(folder_path)
    list_of_csv_files=help_functions.get_file_name_list('xyz_') #changed to 'xyz_' from .csv
    for csv_file_name in list_of_csv_files:
        xyz_df=help_functions.convert_file_to_xyz_df(csv_file_name)
        new_file_name=help_functions.change_filetype(csv_file_name, new_type=help_functions.FileExtensions.XYZ.value)
        help_functions.dataframe_to_xyz(xyz_df, new_file_name.replace('xyz_','txt_'))
    os.chdir('../')
    return 


def move_xyz_files_directory(current_directory,new_directory):#my help function
    """
    a function that moves xyz type files from one directory to another.
    help function for xyz_file_generator_library to move files to the new directory created
    A void function
    """
    os.chdir(current_directory)
    list_of_xyz_files=help_functions.get_file_name_list(help_functions.FileExtensions.XYZ.value)
    for xyz_file_name in list_of_xyz_files:
        current_path=os.path.join(current_directory, xyz_file_name)
        new_path=os.path.join(new_directory, xyz_file_name)
        os.replace(current_path,new_path)
    return
                       
def xyz_file_generator_library(files_directory_path, directory_name): ###not working, name in R:'xyz_file_generator_library'
    """
    a void function

    """
    path=os.path.join(files_directory_path,directory_name)
    try:
        os.mkdir(path)
        xyz_file_generator(files_directory_path)            ## edit to one function
        move_xyz_files_directory(files_directory_path,path)
    except FileExistsError:
        xyz_file_generator(files_directory_path)
        move_xyz_files_directory(files_directory_path,path)
    os.chdir('../')
    return

def change_file_name(files_directory_path,old_file_name,new_file_name):###works, name in R: 'name_changer' 
    """
    a function that gets a directory of the desired file, the name of the file to change, and changes it to a new specified name
    A void function
    """
    os.chdir(files_directory_path)
    list_of_molecules=[file for file in os.listdir(files_directory_path) ]
    for molecule in list_of_molecules:
        if (molecule==old_file_name):
            os.rename(molecule,new_file_name)
    return

def calc_angle(p1, p2): ###works, name in R: 'angle' , radians
    dot_product=np.dot(p1, p2)
    norm_x=np.linalg.norm(p1)
    norm_y=np.linalg.norm(p2)
    thetha=np.arccos(dot_product/(norm_x*norm_y))
    return thetha
                                                   

    
def calc_new_base_atoms(atom_indexes,coordinates_array):  #help function for calc_coordinates_transformation
    indexes=np.array(atom_indexes)-1
    if (len(atom_indexes)==4):
        new_origin=(coordinates_array[indexes[0]]+coordinates_array[indexes[1]])/2
    else:
        new_origin=coordinates_array[indexes[0]]
        new_y=(coordinates_array[indexes[-2]]-new_origin)/np.linalg.norm((coordinates_array[indexes[-2]]-new_origin))
        coplane=((coordinates_array[indexes[-1]]-new_origin)/np.linalg.norm((coordinates_array[indexes[-1]]-new_origin)))
    return (new_origin,new_y,coplane)


def calc_basis_vector(origin,y,coplane):#help function for calc_coordinates_transformation
    cross_y_plane=np.cross(coplane,y)
    coef_mat=np.vstack([y, coplane, cross_y_plane])
    angle_new_y_coplane=calc_angle(coplane,y)
    cop_ang_x=angle_new_y_coplane-(np.pi/2)
    result_vector=[0,np.cos(cop_ang_x),0]
    new_x=np.linalg.solve(coef_mat,result_vector)
    new_z=np.cross(new_x,y)
    new_basis=np.vstack([new_x, y, new_z])
    return new_basis

                                             
def calc_coordinates_transformation(coordinates_array,base_atoms_indexes):#origin_atom, y_direction_atom, xy_plane_atom
    """
    a function that recives coordinates_df and new base_atoms_indexes to transform the coordinates by
    and returns a dataframe with the shifted coordinates
    parameters:
    ----------
    coordinates_df: pd.dataframe
        xyz molecule dataframe
    base_atoms_indexes: list of nums
        indexes of new atoms to shift coordinates by.
 
    returns:
        transformed xyz molecule dataframe
    -------
        
    example:
    -------
    calc_coordinates_transformation(coordinates_array,[2,3,4])
    
    Output:
        atom       x       y       z
      0    H  0.3477 -0.5049 -1.3214
      1    B     0.0     0.0     0.0
      2    B    -0.0  1.5257     0.0
    """
    
    new_basis=calc_basis_vector(*calc_new_base_atoms(base_atoms_indexes,coordinates_array))
    new_origin=coordinates_array[base_atoms_indexes[0]-1] ### check with shahar cause this line make no sense for 4 atom coordinates. this line exists in R.
    transformed_coordinates=np.array([np.dot(new_basis,(row-new_origin)) for row in coordinates_array]).round(4)
    return transformed_coordinates
    

def calc_npa_charges(coordinates_array,charges_array,base_atoms_indexes,sub_atoms=None):##added option for subunits
    """
    a function that recives coordinates and npa charges, transform the coordinates
    by the new base atoms and calculates the dipole in each axis
    
    Parameters
    ---------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    charges_array: np.array
        array of npa charges
    base_atoms_indexes:list
        3/4 atom indexes for coordinates transformation
        
    optional-sub_atoms:list
        calculate npa charges from a set of sub_atoms instead of all atoms.       
    Returns:
    -------
    dip_df=calc_npa_charges(coordinates_array,charges,base_atoms_indexes,sub_atoms)
    Output:
    dip_df : pd.DataFrame
        output:            dip_x     dip_y     dip_z     total
                       0  0.097437 -0.611775  0.559625  0.834831
    """
    if sub_atoms==None:
        dip_xyz=np.vstack([(row[0]*row[1]) for row in list(zip(coordinates_array,charges_array))])
    else:
        dip_xyz=np.vstack([(row[0]*row[1]) for row in np.array(list(zip(coordinates_array,charges_array)),dtype=object)[sub_atoms,:]])
    dip_vector=np.sum(dip_xyz,axis=0)
    array_dipole=np.hstack([dip_vector,np.linalg.norm(dip_vector)])
    dip_df=pd.DataFrame(array_dipole,index=['dip_x','dip_y','dip_z','total']).T
    return dip_df

def calc_dipole_gaussian(coordinates_array,dipole_array,base_atoms_indexes):
    basis_vector=calc_basis_vector(*calc_new_base_atoms(base_atoms_indexes,coordinates_array))
    dipole_array[0,0:3]=np.matmul(basis_vector,dipole_array[0,0:3])
    return pd.DataFrame(dipole_array,columns=[['dip_x','dip_y','dip_z','total']])

def convert_atoms_indexes_to_bonds_length(coordinates_array,atoms_indexes): ##for calc_angle_between_atoms
    indexes=np.array(atoms_indexes)-1 #three atoms-angle four atoms-dihedral
    if len(indexes)==3:
        indexes=[indexes[0],indexes[1],indexes[1],indexes[2]]
    else:
        indexes=[indexes[0],indexes[1],indexes[1],indexes[2],indexes[2],indexes[3]]
    index_pair =list(zip(indexes[::2],indexes[1::2]))
    return [(coordinates_array[pair[0]]-coordinates_array[pair[1]]) for pair in index_pair]
  


def calc_angle_between_atoms(coordinates_array,atoms_indexes): #gets a list of atom indexes
    """
    a function that gets 3/4 atom indexes, and returns the angle between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    atoms_indexes- list of ints
        a list of atom indexes to calculate the angle between- [2,3,4]
   
    Returns
    -------
    angle: float
        the bond angle between the atoms
    """
    bonds_list=convert_atoms_indexes_to_bonds_length(coordinates_array,atoms_indexes)
    if len(atoms_indexes)==3:
        angle=calc_angle(*bonds_list)*(180/math.pi)
    else:
        first_cross=np.cross(bonds_list[0],bonds_list[1])
        second_cross=np.cross(bonds_list[2],bonds_list[1]) 
        angle=calc_angle(first_cross,second_cross)*(180/math.pi)
    return angle

def calc_bond_length(coordinates_array,atom_pairs): 
    """
    a function that calculates the distance between each pair of atoms.
    help function for molecule class
    
    Parameters
    ----------
    coordinates_array: np.array
        xyz coordinates 
        
    atom_pairs : iterable
        list containing atom pairs-(([2,3],[4,5]))
        
    Returns
    -------
    pairs_df : dataframe
        distance between each pair 
        
        Output:
                               0
bond length[2, 3]            1.525692
bond length[4, 5]            2.881145
    
    """
    pairs=np.array(atom_pairs)-1
    index=[('bond length')+str(pair) for pair in pairs]
    bond_list=[np.linalg.norm(coordinates_array[pair[0]]-coordinates_array[pair[1]]) for pair in pairs]
    pairs_df=pd.DataFrame(bond_list,index=index)
    return pairs_df

def organize_info_df(info_df):
    """
    this function takes the starting dataframe from info file and
    organizes it for calculations
    
    Parameters
    ----------
    info_df : dataframe
        dataframe created with help_function from csv file 'info'.

    Returns
    -------
    ordered_info_df : dataframe

    """
    info_df.set_axis(info_df[0],inplace=True)
    frequencies=info_df.loc['Frequencies'][[2,3,4]]
    ir=info_df.loc['IR'][[3,4,5]]
    frequencies_list = [item for sublist in frequencies.values.tolist() for item in sublist]#ugly
    ir_list=[item for sublist in ir.values.tolist() for item in sublist]
    ordered_info_df=pd.DataFrame([frequencies_list,ir_list], index=['Frequency[1/cm]','IR intensity'],dtype=float)
    return ordered_info_df

def vib_df_to_coordinates_array(vib_df): ##name in R: atom.vectors
    """
    a function that takes unedited vib df and organizes it to xyz coordinates
    
    Parameters
    ----------
    vib_df : dataframe
        dataframe get_unedited_df

    Returns
    -------
    vib_array : array containing xyz coordinates
    
    example:
        input:
            0   10  1   0.00   0.00   0.04   0.01  -0.01   0.30  -0.02   0.02  -0.20
            
        vib_df_to_coordinates_array(molecule_1.get_unedited_df('vib_1'))
        
        array([[ 0.  ,  0.  ,  0.04],
               [ 0.01, -0.01,  0.3 ],
               [-0.02,  0.02, -0.2 ]
        
    """
    vib_array=np.array((vib_df.drop([0,1],axis=1)).astype(float))
    return vib_array.reshape(-1, 3)

def calc_magnituede_from_coordinates_array(coordinates_array): ##help function for get_vib_info method
    return [np.linalg.norm(row) for row in coordinates_array]
      
      
def calc_max_frequency_magnitude(vib_array, info_df):##add option to return ordered_info_df-like dot.prod.info
    """
    a function that gets vibration and info dataframes and returns
    the frequecy and IR with the max magnitude for the vibration.
    splits the coordinates of vib to 3 coordinates and calculates the magnituede. takes frequencys greter than 1500
    and returns the frequency and IR corresponding to max magnitude.
    
    works in a specific molecule directory
    
    Parameters
    ----------
    vib_array: np.array
        organized vib file in array.
    info_df: np.dataframe
        organized info file in dataframe. 
        
    Returns
    -------
    dataframe
        max frequency and IR for specific vibration.
        
    Output:
                                            54
                        Frequency[1/cm]  1689.5945
                        IR intensity        6.5260
    """
    
    magnitude=calc_magnituede_from_coordinates_array(vib_array)
    df=(info_df.T)
    df['magnitude']=magnitude
    outer_finger=(df['Frequency[1/cm]'].astype(float)>1500)
    index_max=df[outer_finger]['magnitude'].idxmax()
    return info_df[index_max]

def organize_pol_info(pol_df):
    organized_df=(pol_df.iloc[4:6])
    values=[float(value.split('D')[0])*1000 for value in list(organized_df[1])]
    return pd.DataFrame(values,index=organized_df[0])

def check_pair_in_bonds(atom_pairs,bonds_df): ##help functions for gen_vib
    """
    a function that checks that the all given atom pair exists as a bond in the bonds_df
    """
    bonds_list=(bonds_df.astype(int)).values.tolist()   
    bool_check=all([pair in bonds_list for pair in atom_pairs])
    return bool_check

def atom_pairs_to_vib_file_names(atom_pairs):
    """
    a function that gets a list of atoms pairs- [[1,2],[3,4],[5,6]]

    Parameters
    ----------
    atom_pairs : list
        a list of atom pairs.

    Returns
    -------
    vib_file_names : string
        name of the vibration file.
        
        atom_pairs_to_vib_file_names([[1,3],[4,8]])
        output: 
        ['vib_1_A1_o.csv', 'vib_3_A1_o.csv', 'vib_4_A1_o.csv', 'vib_8_A1_o.csv']

    """
    atom_indexes=[atom for sublist in atom_pairs for atom in sublist] 
    vib_name=(help_functions.get_file_name_list('vib')[0]) 
    vib_file_names=[vib_name.replace(vib_name.split('_')[1],str(index)) for index in atom_indexes]
    return vib_file_names

def vib_names_to_array_pair_list(vib_file_names):
    """
    a function that gets a list of vib file names and returnst

    Parameters
    ----------
    vib_file_names : TYPE
        DESCRIPTION.

    Returns
    -------
    vib_array_pairs ([0]) : list
    containing the arrays from the vib file, each variable containing two
    arrays corresponding to the bond pairs 
        .
    vib_array_list ([1]): list
        containing all vib array coordinates.

    """
    vib_df_list=[help_functions.get_df_from_csv(name) for name in vib_file_names] 
    vib_array_list=[vib_df_to_coordinates_array(df) for df in vib_df_list]
    vib_array_pairs=list(zip(vib_array_list[::2], vib_array_list[1::2]))
    return vib_array_pairs,vib_array_list ##optional- do not split to pairs

def atom_pairs_to_coordinates_vector(coordinates_array,atom_pairs):
    atom_pairs_indexes=np.array(atom_pairs)-1
    coordinates_vector=[(coordinates_array[pair[0]]-coordinates_array[pair[1]]) for pair in atom_pairs_indexes]
    return coordinates_vector

def calc_vib_dot_product(vib_array_pairs,coordinates_vector):
    vib_dot_product=[(abs(np.dot(pair[0],coordinates_vector[i]))+abs(np.dot(pair[1],coordinates_vector[i]))) for i,pair in list(enumerate(vib_array_pairs))]
    return vib_dot_product

def calc_vib_dot_product_from_pairs(atom_pairs,coordinates):
    vib_file_names=atom_pairs_to_vib_file_names(atom_pairs)
    vib_array_pairs=vib_names_to_array_pair_list(vib_file_names)[0]
    coordinates_vector=atom_pairs_to_coordinates_vector(coordinates,atom_pairs)
    vib_dot_product=calc_vib_dot_product(vib_array_pairs,coordinates_vector)
    return vib_dot_product

def calc_max_frequency_gen_vib(vib_dot_product,info_df,threshhold=1600): ##last for gen_vib
    df_list=[pd.DataFrame(np.vstack([vib,(info_df.T)['Frequency[1/cm]']])) for vib in vib_dot_product]
    threshhold_mask=[df.T[(df.iloc[1]>threshhold)] for df in df_list]
    index_max=[df[0].idxmax() for  df in threshhold_mask]
    max_frequency_vib=[pd.DataFrame(list(df[index_max[i]]),index=[['IR','Frequency']]) for i,df in enumerate(df_list)]
    return max_frequency_vib


def vib_ring_array_list_to_vector(vib_array_list): ##help function ring_vib
    vec_sum_1_3_5=(vib_array_list[0]+vib_array_list[2]+vib_array_list[3])
    vec_sum_2_4_6=(vib_array_list[4]+vib_array_list[1]+vib_array_list[5])
    return vec_sum_1_3_5, vec_sum_2_4_6

def get_data_for_ring_vib(info_df,vib_array_list,coordinates_vector):
    prods=[np.dot(row_1,row_2) for row_1,row_2 in zip(*vib_ring_array_list_to_vector(vib_array_list))]
    calc=[abs(math.sin(calc_angle(row,coordinates_vector))) for row in vib_ring_array_list_to_vector(vib_array_list)[0]]
    data_df=pd.DataFrame(np.vstack([prods,(info_df.T)['Frequency[1/cm]'],calc]))
    return data_df

def get_filter_ring_vib_df(data_df): ##return filtered_df and check 
    filtered_df=data_df.T[(data_df.iloc[0]!=0)&(abs(data_df.iloc[0])>0.1)][data_df.iloc[1]>500].reset_index() ##randon threshhold to work, no >1500
    if (filtered_df.shape[0]==0):
        filtered_df=data_df.T[(data_df.iloc[0]!=0)][(data_df.iloc[1]>500)&(data_df.iloc[1]<700)].reset_index()
        print('Dot products are lower than 0.1 - returning the default 500 - 700 1/cm')
    return filtered_df

def calc_min_max_ring_vib(filtered_df):
    max_vib_frequency=filtered_df.iloc[filtered_df[2].idxmax()][1]
    asin_max=math.asin(filtered_df[2].max())*(180/np.pi)
    min_vib_frequency=filtered_df.iloc[filtered_df[2].idxmin()][1]
    asin_min=math.asin(filtered_df[2].min())*(180/np.pi)
    df=pd.DataFrame((max_vib_frequency,asin_max,min_vib_frequency,asin_min),index=['cross','cross_angle','para','para_angle'])
    return df.T

def get_filtered_ring_df(info_df,coordinates,ring_atom_indexes): ## last ring_vib
    coordinates_vector=atom_pairs_to_coordinates_vector(coordinates,ring_atom_indexes)[0]
    vib_file_names=atom_pairs_to_vib_file_names(ring_atom_indexes)
    vib_array_list=vib_names_to_array_pair_list(vib_file_names)[1]
    data_df=get_data_for_ring_vib(info_df,vib_array_list,coordinates_vector)
    filtered_df=get_filter_ring_vib_df(data_df)
    return filtered_df

######## sterimol functions

def direction_atoms_for_sterimol(bonds_df,base_atoms): #help function for sterinol
    """
    a function that return the base atom indexes for coordination transformation according to the bonded atoms.
    you can insert two atom indexs-[1,2] output [1,2,8] or three if the first one repeats-[1,2,1] output [1,2,3]
    """
    try:
        origin,direction=base_atoms[0],base_atoms[1]
        try :
            base_atoms[2]==origin
            if(any(bonds_df[0]==direction)):
                base_atoms[2]=int(bonds_df[(bonds_df[0]==direction)][1].iloc[1])
            else:
                base_atoms[2]=int(bonds_df[(bonds_df[1]==direction)][0].iloc[1])
        except: 
            if (any(bonds_df[0]==direction)):
                base_atoms.append(int(bonds_df[(bonds_df[0]==direction)][1].iloc[0]))
            else:
                base_atoms.append(int(bonds_df[(bonds_df[1]==direction)][0].iloc[0]))
        return base_atoms
    except IndexError:
        return print('no such direction')

def get_molecule_connections(bonds_df,source,direction):
 
    graph=ig.Graph.DataFrame(edges=bonds_df,directed=True)
    paths=graph.get_all_simple_paths(v=source,mode='all')
    with_direction=[path for path in paths if (direction in path)]
    return np.unique([path for sublist in with_direction for path in sublist] )

def get_tc_for_sterimol(degree,plane):
    cos_deg=np.cos(degree*(np.pi/180))
    sin_deg=np.sin(degree*(np.pi/180))
    rot_matrix=np.array([[cos_deg,-1*sin_deg],[sin_deg,cos_deg]])
    tc=np.vstack([np.matmul(rot_matrix,row) for row in plane]).round(3)
    
    return tc

def get_specific_bonded_atoms_df(bonds_df,atom_filter,coordinates_df):
    edited_bonds_df=bonds_df[(bonds_df.isin(atom_filter))].dropna().reset_index(drop=True)
    bonds_array=(np.array(edited_bonds_df)-1).astype(int)
    atom_bonds=np.vstack([(coordinates_df.iloc[bond]['element'].values) for bond in bonds_array]).reshape(-1,2)
    bonded_atoms_df=(pd.concat([pd.DataFrame(atom_bonds),edited_bonds_df],axis=1))
    bonded_atoms_df.columns=[['atom_1','atom_2','index_1','index_2']]
    return bonded_atoms_df

def remove_atom_bonds(atom_bonds_df,atom_remove='H'):
    delete_rows=[row[0] for row in atom_bonds_df.iterrows() if row[1][0]==atom_remove]
    atom_bonds_df.drop(delete_rows,axis=0,inplace=True)
    return atom_bonds_df


def remove_nof_bonds(bonded_atoms_df,bonds_array,coordinates_df):
    nof_indexes=[row[0] for row in bonded_atoms_df.iterrows() if (row[1][1]=='H') & (row[1][0] in ('N.1','N.2','N.3','N.4','N.am','N.ar','N.pl3','O.2','O.3','F'))]
    for index in nof_indexes:
        pair=bonds_array[index]
        coordinates=coordinates_df[['x','y','z']].astype(float)
        bond_lenght=(np.linalg.norm(coordinates.iloc[pair[0]]-coordinates.iloc[pair[1]]))
        if bond_lenght>1.6:
            bonded_atoms_df.drop(index,axis=0,inplace=True)
    
def filter_atoms_for_sterimol(bonded_atoms_df,coordinates_df):### filtering nof bonds and H bonds *****deleting H only from the left column 
    bonds_array=np.array(bonded_atoms_df[['index_1','index_2']]).astype(int)-1
    bonded_atoms_df=remove_atom_bonds(bonded_atoms_df)
    remove_nof_bonds(bonded_atoms_df, bonds_array, coordinates_df)
          
    clean_bonds= pd.concat([bonded_atoms_df['index_1'],bonded_atoms_df['index_2']],axis=1).reset_index(drop=True)
    atom_filter=np.unique([atom for sublist in clean_bonds.values.tolist() for atom in sublist])-1
    edited_coordinates_df=coordinates_df.iloc[atom_filter].reset_index(drop=True)
    return edited_coordinates_df


def get_extended_df_for_sterimol(coordinates_df,bond_type,radii):
    """
    a function that adds information to the regular coordinates_df
    Parameters
    ----------
    coordinates_df : dataframe
        
    radii : dictionary
        CPK or bondi. deafult used in the methos is 'bondi'

    Returns
    -------
    coordinates_df : dataframe
       
    input:
        atom             x             y               z
     0     H   4.175633987  -0.285908497  -0.962856209\n
     1     B   2.804633987   0.066091503  -1.306856209\n

    """
    coordinates_df['atype']=bond_type
    coordinates_df['atype'].replace(GeneralConstants.BOND_TYPE.value, inplace=True)
    coordinates_df['magnitude']=calc_magnituede_from_coordinates_array(np.array(coordinates_df[['x','z']].astype(float)))
    if radii=='bondi':
        coordinates_df['radius']=coordinates_df['atype'].replace(GeneralConstants.BONDI_RADII.value)
    elif radii=='CPK':
        coordinates_df['radius']=coordinates_df['atype'].replace(GeneralConstants.CPK_RADII.value)
    coordinates_df['Bs']=coordinates_df[['radius','magnitude']].sum(axis=1)
    coordinates_df['L']=coordinates_df[['y','radius']].sum(axis=1)
    return coordinates_df

def calc_B1(tc,avs,edited_coordinates_df,column_index):
    ind=np.where(np.abs(tc[:,column_index])==(avs.min()).round(3))[0][0]
    if tc[ind,column_index]<0:
        new_ind=np.where(tc[:,column_index]==tc[:,column_index].min())[0][0]
        bool_list=np.logical_and(tc[:,column_index]>=tc[new_ind,column_index], tc[:,column_index]<=tc[new_ind,column_index]+1)
        tc[:,column_index]=-tc[:,column_index]
    else:
        bool_list=np.logical_and(tc[:,column_index]>=tc[ind,column_index]-1, tc[:,column_index]<=tc[ind,column_index])
        
    against,against_loc=[],[]
    B1,B1_loc=[],[]
    for i in range(1,tc.shape[0]): ######whyyy
        if bool_list[i]:
            against.append(np.array(tc[i,column_index]+edited_coordinates_df['radius'].iloc[i]))
            against_loc.append(edited_coordinates_df['L'].iloc[i])
        if against:
            B1.append(max(against))
            B1_loc.append(against_loc[against.index(max(against))])
            
    return [B1,B1_loc]

def get_b1s_list(extended_df):
    b1s,b1s_loc=[],[]
    scans=90//5
    degree_list=list(range(18,108,scans))
    plane=np.array(extended_df[['x','z']].astype(float))
    for degree in degree_list:
        tc=get_tc_for_sterimol(degree, plane)
        avs=np.abs([max(tc[:,0]),min(tc[:,0]),max(tc[:,1]),min(tc[:,1])])
        if np.where(avs==avs.min())[0][0] in [0,1]:
            B1,B1_loc=calc_B1(tc,avs,extended_df,0)
        elif np.where(avs==avs.min())[0][0] in [2,3]:
            B1,B1_loc=calc_B1(tc,avs,extended_df,1)
        b1s.append(np.unique(np.vstack(B1)).max())####check
        b1s_loc.append(np.unique(np.vstack(B1_loc)).max())
               
    back_ang=degree_list[np.where(b1s==min(b1s))[0][0]]-scans   
    front_ang=degree_list[np.where(b1s==min(b1s))[0][0]]+scans                      
    
    for i in range(back_ang,front_ang+1):
        tc=get_tc_for_sterimol(i,plane)
        avs=np.abs([max(tc[:,0]),min(tc[:,0]),max(tc[:,1]),min(tc[:,1])])
        if np.where(avs==avs.min())[0][0] in [0,1]:
            B1,B1_loc=calc_B1(tc,avs,extended_df,0)
        elif np.where(avs==avs.min())[0][0] in [2,3]:
            B1,B1_loc=calc_B1(tc,avs,extended_df,1)

        b1s.append(np.unique(np.vstack(B1)).max())####check
        b1s_loc.append(np.unique(np.vstack(B1_loc)).max())
    
    return [np.array(b1s),np.array(b1s_loc)]



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
        delete_type_files('xyz')
        self.list_of_files=[file_name for file_name in os.listdir(self.molecule_path)]
        self.coordinates_df=[help_functions.convert_file_to_xyz_df(file_name) for file_name in self.list_of_files if 'xyz_' in file_name][0] #fix atomic number ,dtype
        # self.dfs=pd.DataFrame([help_functions.convert_file_to_xyz_df(file_name) for file_name in self.list_of_files],index=self.list_of_files)
        os.chdir('../')
     
    def get_file_name(self,file_type):
        return [file_name for file_name in self.list_of_files if file_type in file_name][0]
    
    def get_unedited_df(self,file_type,splitter=None): ##info df will not be ordered, 'vib','xyz' need ','
        os.chdir(self.molecule_path)
        df=help_functions.get_df_from_csv([file_name for file_name in self.list_of_files if file_type in file_name][0],splitter=splitter)
        # os.chdir('../')
        return df
    
    def get_sterimol(self,base_atoms,radii='bondi'):
        bonds_df=self.get_unedited_df('bonds',',').astype(int)
        bonds_direction=direction_atoms_for_sterimol(bonds_df,base_atoms)
        coordinates_df=self.get_coordination_transformation_df(bonds_direction)
        
        connected_from_direction=get_molecule_connections(bonds_df,bonds_direction[0],bonds_direction[1])
        ## creating df of atoms and bonds
        bonded_atoms_df=get_specific_bonded_atoms_df(bonds_df,connected_from_direction,self.coordinates_df)     
         ### filtering nof bonds and H bonds *****deleting H only from the left column       
        edited_coordinates_df=filter_atoms_for_sterimol(bonded_atoms_df,coordinates_df)
        ##adding colums
        extended_df=get_extended_df_for_sterimol(edited_coordinates_df,self.get_unedited_df('atype'),radii)
        ###calculations
        b1s,b1s_loc=get_b1s_list(extended_df)
        B1=min(b1s[b1s>=0])
        loc_B1=max(b1s_loc[np.where(b1s[b1s>=0]==min(b1s[b1s>=0]))])
        B5=max(extended_df['Bs'].values)
        L=max(extended_df['L'].values+0.4)
        loc_B5=min(extended_df['y'].iloc[np.where(extended_df['Bs'].values==B5)[0][0]])
        df=pd.DataFrame([B1,B5,L,loc_B5,loc_B1],index=[['B1','B5','L','loc_B5','loc_B1']]) ## loc_B1 produces different values
      
        return df
        
    def export_coordinates_to_xyz(self):
        os.chdir(self.molecule_path)
        output_name=help_functions.get_file_name_list('xyz')[0]
        help_functions.dataframe_to_xyz(self.coordinates_df,help_functions.change_filetype(output_name,'_tc.xyz'))
        os.chdir('../')
        
    def swap_atom_pair(self,pair_index): #swapping is permanently changed
        pairs=np.array(pair_index)-1
        temp=self.coordinates_df.iloc[pairs[0]].copy()
        self.coordinates_df.iloc[pairs[0]]=self.coordinates_df.iloc[pairs[1]]
        self.coordinates_df.iloc[pairs[1]]=temp
        return self.coordinates_df
        
    def get_coordination_transformation_df(self,base_atoms_indexes):
        coordinates_array=np.array(self.coordinates_df[['x','y','z']].astype(float))  
        transformed_coordinates= calc_coordinates_transformation(coordinates_array,base_atoms_indexes)
        array_for_df=(np.column_stack([np.array(self.coordinates_df['atom']),transformed_coordinates]))
        return pd.DataFrame(array_for_df,columns=help_functions.XYZConstants.DF_COLUMNS.value)
    
    def get_npa_df(self,base_atoms_indexes,sub_atoms=None):
        charges=np.array(self.get_unedited_df('npa').astype(float))
        coordinates_array=np.array(self.get_coordination_transformation_df(base_atoms_indexes)[['x','y','z']].astype(float))
        dip_df=calc_npa_charges(coordinates_array,charges,base_atoms_indexes,sub_atoms)
        dip_df.rename(index={0:self.get_file_name('npa')},inplace=True)
        return dip_df
    
    
    def get_dipole_gaussian_df(self,base_atoms_indexes):
        dipole_array=np.array(self.get_unedited_df('dipole',',').astype(float))
        coordinates_array=np.array(self.coordinates_df[['x','y','z']].astype(float))
        return calc_dipole_gaussian(coordinates_array,dipole_array,base_atoms_indexes)       
    
    def get_bond_angle(self,atom_indexes):##need expand to many angles
        coordinates=np.array(self.coordinates_df[['x','y','z']].astype(float))
        return [calc_angle_between_atoms(coordinates,atom_index) for atom_index in atom_indexes]
    
    
    def get_bond_length(self,atom_pairs):##input-([2,3],[4,5]) doesnt work with one input
        coordinates=np.array(self.coordinates_df[['x','y','z']].astype(float))
        bond_df= calc_bond_length(coordinates,atom_pairs)
        bond_df.rename(columns={0:self.get_file_name('xyz')},inplace=True)
        return bond_df
    
    def get_info_df(self):
        info_df=self.get_unedited_df('info')
        return organize_info_df(info_df)
    
    def get_vib_max_frequency(self,vib_file_name):
        vib_array=vib_df_to_coordinates_array(self.get_unedited_df(vib_file_name))
        info_df=self.get_info_df()
        return calc_max_frequency_magnitude(vib_array,info_df)
    
    def get_pol_info(self):
        return(organize_pol_info(self.get_unedited_df('Pol')))
    
    def get_gen_vib(self,atom_pairs):
        """

        Parameters
        ----------
        atom_pairs : molecule_1.get_gen_vib([[1,6],[3,4]])
            atom pairs must have a corresponding vib file and appear in bonds_df.

        Returns
        -------
        dataframe
            
        [                   0
         IR            0.7310
         Frequency  1655.5756,
                             0
         IR            0.35906
         Frequency  1689.59450]

        """
        os.chdir(self.molecule_path)
        bonds_df=self.get_unedited_df('bonds',',')
        if check_pair_in_bonds(atom_pairs,bonds_df)==True:
            coordinates=np.array(self.coordinates_df[['x','y','z']].astype(float))
            info_df=self.get_info_df()
            vib_dot_product=calc_vib_dot_product_from_pairs(atom_pairs,coordinates)
            return calc_max_frequency_gen_vib(vib_dot_product,info_df)
        else:
            print('Error: the following bonds do not exist-check atom numbering')
            
    def get_ring_vibs(self,ring_atom_indexes):
        """

        Parameters
        ----------
        ring_atom_indexes :working example: molecule_1.get_ring_vibs([[1,3],[1,6],[3,4]]) 
            a list of atom pairs, there must be a vib file with a corresponding number to work.
            
            # Enter a list of ring atoms with the order: primary axis - para followed by primary,
            # ortho - ortho atoms relative to primary atom and meta - meta atoms relative to primary.
            # For example - for a ring of atoms 1-6 where 4 is connected to the main group and 1 is para to it
            # (ortho will be 3 & 5 and meta will be 2 & 6) - enter the input [[1,4],[3,5],[2,6]].
        Returns
        -------
        dataframe
            cross  cross_angle      para  para_angle
      0  657.3882    81.172063  834.4249   40.674833

        """
        os.chdir(self.molecule_path)
        info_df=self.get_info_df()
        coordinates=np.array(self.coordinates_df[['x','y','z']].astype(float))
        try:
            filtered_df=get_filtered_ring_df(info_df,coordinates,ring_atom_indexes)
        except FileNotFoundError:
            return print('no vibration data file exists for this atom')
        # bool_check=filtered_df[2].duplicated().any() #check for duplicates in calc, is it needed ?? slim chances
        return calc_min_max_ring_vib(filtered_df)
    
    def check_imaginary_frequency(self):##return True if no complex frequency, called ground.state in R
        info=self.get_info_df()
        bool_imaginary=not any([isinstance(frequency, complex) for frequency in info.loc['Frequency[1/cm]']])
        return bool_imaginary
    
    
    def get_nbo_info(self,atom_pairs):
        pairs_indexes=np.array(atom_pairs)-1
        index=[('nbo charge')+str(pair) for pair in pairs_indexes]
        nbo_array= np.array((self.get_unedited_df('nbo').drop([0,1],axis=1)).astype(float))
        nbo_diff=[abs(statistics.mean(nbo_array[pair[0]])-statistics.mean(nbo_array[pair[1]])) for pair in pairs_indexes]
        return pd.DataFrame(nbo_diff, index=index)
    
    def get_molecule_comp_set(self,dipole_mode = 'gaussian', radii = 'bondi'):
        """
      molecule_1.get_molecule_comp_set()
      Ring atoms - by order -> primary axis (para first), ortho atoms and meta atoms: 1 3 1 6 3 4
      __main__:193: RuntimeWarning: invalid value encountered in double_scalars
      __main__:544: UserWarning: Boolean Series key will be reindexed to match DataFrame index.
      your atom pairs: 1 6 3 4
      Enter atoms - origin atom, y axis atom and xy plane atom: 1 2 3
      Insert atom pairs for which you wish calculate differences: 1 2 3 4
      Primary axis along: 1 2
      Distances - Atom pairs: 1 2 3 4
      Do you want to compute any angles/dihedrals? y/n: y
      Insert a list of atom triads/quartets for which you wish to have angles/dihedrals.
      For several angles/dihedrals, insert triads/quartets with a double space between them.
      Make sure there are spaces between atoms as well.
      For example - 1 2 3  1 2 3 4 will give angle(1, 2, 3) and dihedral(1, 2, 3, 4)
      Insert a list of atom triads/quartets for which you wish to have angles/dihedrals: 1 2 3  1 2 3 4
        -------
        res : TYPE
            [[  cross  cross_angle      para  para_angle
              0  657.3882    81.172063  834.4249   40.674833
              
              [                   0
               IR            0.7310
               Frequency  1655.5756,
                                   0
               IR            0.35906
               Frequency  1689.59450],
                    dip_x     dip_y     dip_z   total
              0  0.058148  1.210907  0.949174  1.5397,
              [0.2706980000000003, 0.013509999999999689],
                                        0
              B1                    2.556
              B5      [4.131529971026473]
              L                  [5.8357]
              loc_B5               1.4121
              loc_B1               3.6855,
                                xyz_molecule_1.csv
              bond length[0 1]            1.456668
              bond length[2 3]            1.227760],
             [69.7212447900024, 104.74132137986373]]

        """
        res=[]
        ring_atom_index=help_functions.split_to_pairs(input("Ring atoms - by order -> primary axis (para first), ortho atoms and meta atoms: "))
        ring_df=self.get_ring_vibs(ring_atom_index)
        gen_vib_index=help_functions.split_to_pairs(input("your atom pairs: "))
        gen_vib_df=self.get_gen_vib(gen_vib_index)
        
        dipole_index=help_functions.string_to_int_list(input('Enter atoms - origin atom, y axis atom and xy plane atom: '))
        if dipole_mode=='compute':
            dipole_df=self.get_npa_df(dipole_index)
        elif dipole_mode=='gaussian':
            dipole_df=self.get_dipole_gaussian_df(dipole_index)
        nbo_index=help_functions.split_to_pairs(input('Insert atom pairs for which you wish calculate differences: '))
        nbo_df=self.get_nbo_info(nbo_index)
        sterimol_index=help_functions.string_to_int_list(input('Primary axis along: '))
        if radii=='CPK':
            sterimol_df=self.get_sterimol(sterimol_index,radii)
        else:
           sterimol_df=self.get_sterimol(sterimol_index) 
        index_for_dist=help_functions.split_to_pairs(input('Distances - Atom pairs: '))
        atom_dist_df=self.get_bond_length(index_for_dist)
        res.append([ring_df,gen_vib_df,dipole_df,nbo_df,sterimol_df,atom_dist_df])
        angle_answer=input('Do you want to compute any angles/dihedrals? y/n: ')
        if angle_answer=='y':
            print('Insert a list of atom triads/quartets for which you wish to have angles/dihedrals.\n'
        'For several angles/dihedrals, insert triads/quartets with a double space between them.\n'
        'Make sure there are spaces between atoms as well.\n'
        'For example - 1 2 3  1 2 3 4 will give angle(1, 2, 3) and dihedral(1, 2, 3, 4)')
            angle_string=help_functions.split_for_angles(input('Insert a list of atom triads/quartets for which you wish to have angles/dihedrals: '))
            angle_df=self.get_bond_angle(angle_string)
            res.append(angle_df)
        return res
            
    
class Molecules():
    def __init__(self,molecules_dir_name):
        self.molecules_path=os.path.abspath(molecules_dir_name)
        os.chdir(self.molecules_path) 
        self.molecules=[Molecule(molecule_dir) for molecule_dir in os.listdir() if os.path.isdir(molecule_dir)] 
        self
        os.chdir('../')
        
        
if __name__=='__main__':
    # xyz_file_generator_library(r'C:\Users\edens\Documents\GitHub\learning_python\project\main_python','new_directory') #works
    path=r'C:\Users\edens\Documents\GitHub\main_python'
   
    molecules=Molecules('test_dipole')
    molecule_1=molecules.molecules[0]
    molecule_2=molecules.molecules[1]
    os.chdir(path)
    
 
    