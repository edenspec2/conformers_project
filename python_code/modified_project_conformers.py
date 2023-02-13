import pandas as pd
import numpy as np
import os
import math
import help_functions
from enum import Enum
import igraph as ig

from gaussian_handler import gauss_file_handler
import help_functions # personal imports comes after global import

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
    
    BONDI_RADII={
        'H': 1.10, 'C': 1.70, 'F': 1.47,
        'S': 1.80, 'B': 1.92, 'I': 1.98, 
        'N': 1.55, 'O': 1.52, 'Co': 2.00, 
        'Br': 1.83, 'Si': 2.10,'Ni': 2.00,
        'P': 1.80, 'Cl': 1.75, 
    }
    
    CPK_RADII={
        'C':1.50,   'H':1.00,   'S.O':1.70,  'Si':2.10,
        'C2':1.60,  'N':1.50,   'S1':1.00,   'Co':2.00,
        'C3':1.60,  'C66':1.70, 'F':1.35,    'Ni':2.00,
        'C4':1.50,  'N4':1.45,  'Cl':1.75,
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
    REGULAR_BOND_TYPE = {

        'O.2': 'O', 'N.2': 'N', 'S.3': 'S',
        'O.3': 'O', 'N.1': 'N', 'S.O2': 'S',
        'O.co2': 'O', 'N.3': 'N', 'P.3': 'P',
        'C.1': 'C', 'N.ar': 'N',
        'C.2': 'C', 'N.am': 'N',
        "C.cat": 'C', 'N.pl3': 'N',
        'C.3': 'C', 'N.4': 'N',
        'C.ar': 'C', 'S.2': 'S',
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


def adjust_indices(indices,adjustment_num=1):
    """
    adjust indices to start from 0
    """
    return np.array(indices)-adjustment_num

def calc_angle(p1, p2, degrees=False): ###works, name in R: 'angle' , radians
    dot_product=np.dot(p1, p2)
    norm_x=np.linalg.norm(p1)
    norm_y=np.linalg.norm(p2)
    thetha=np.arccos(dot_product/(norm_x*norm_y))
    if degrees:
        thetha=np.degrees(thetha)
        
    return thetha
                                                   

    
def calc_new_base_atoms( coordinates_array,atom_indices):  #help function for calc_coordinates_transformation
    if (len(atom_indices)==4):
        new_origin=(coordinates_array[atom_indices[0]]+coordinates_array[atom_indices[1]])/2
    else:
        new_origin=coordinates_array[atom_indices[0]]
    new_y=(coordinates_array[atom_indices[-2]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[-2]]-new_origin))
    coplane=((coordinates_array[atom_indices[-1]]-new_origin)/(np.linalg.norm((coordinates_array[atom_indices[-1]]-new_origin))))
    
    return (new_origin,new_y,coplane)


def calc_basis_vector(origin, y, coplane):#help function for calc_coordinates_transformation
    """
    origin: origin of the new basis
    y: y direction of the new basis
    coplane: a vector that is coplanar with the new y direction
    """
    cross_y_plane=np.cross(coplane,y)
    coef_mat=np.vstack([y, coplane, cross_y_plane])
    angle_new_y_coplane=calc_angle(coplane,y)
    cop_ang_x=angle_new_y_coplane-(np.pi/2)
    result_vector=[0,np.cos(cop_ang_x),0]

    new_x=np.linalg.solve(coef_mat,result_vector)
    new_z=np.cross(new_x,y)
    new_basis=np.vstack([new_x, y, new_z])
    return new_basis

                                             
def calc_coordinates_transformation(coordinates_array, base_atoms_indices):#origin_atom, y_direction_atom, xy_plane_atom
    """
    a function that recives coordinates_array and new base_atoms_indices to transform the coordinates by
    and returns a dataframe with the shifted coordinates
    parameters:
    ----------
    coordinates_array: np.array
        xyz molecule array
    base_atoms_indices: list of nums
        indices of new atoms to shift coordinates by.
 
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
    indices=adjust_indices(base_atoms_indices)
    new_basis=calc_basis_vector(*calc_new_base_atoms(coordinates_array,indices))
    new_origin=coordinates_array[indices[0]-1] ### check with shahar cause this line make no sense for 4 atom coordinates. this line exists in R.
    transformed_coordinates=np.array([np.dot(new_basis,(row-new_origin)) for row in coordinates_array]).round(4)
    return transformed_coordinates
    

def calc_npa_charges(coordinates_array,charges_array,base_atoms_indices,sub_atoms=None):##added option for subunits
    """
    a function that recives coordinates and npa charges, transform the coordinates
    by the new base atoms and calculates the dipole in each axis
    
    Parameters
    ---------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    charges_array: np.array
        array of npa charges
    base_atoms_indices:list
        3/4 atom indices for coordinates transformation
        
    optional-sub_atoms:list
        calculate npa charges from a set of sub_atoms instead of all atoms.       
    Returns:
    -------
    dipole_df=calc_npa_charges(coordinates_array,charges,base_atoms_indices,sub_atoms)
    Output:
    dipole_df : pd.DataFrame
        output:            dip_x     dip_y     dip_z     total
                       0  0.097437 -0.611775  0.559625  0.834831
    """
    indices=adjust_indices(base_atoms_indices)
    transformed_coordinates=calc_coordinates_transformation(coordinates_array, indices)
    atom_mask=range(charges_array) if sub_atoms==None else sub_atoms
    dipole_xyz=np.vstack([(row[0]*row[1]) for row in np.array(list(zip(transformed_coordinates,charges_array)),dtype=object)[atom_mask,:]])
    dipole_vector=np.sum(dipole_xyz,axis=0)
    array_dipole=np.hstack([dipole_vector,np.linalg.norm(dipole_vector)])
    dipole_df=pd.DataFrame(array_dipole,indices=['dipole_x','dipole_y','dipole_z','total'])
    return dipole_df

def calc_dipole_gaussian(coordinates_array, gauss_dipole_array, base_atoms_indices):
    """
    a function that recives coordinates and gaussian dipole, transform the coordinates
    by the new base atoms and calculates the dipole in each axis
    """
    indices=adjust_indices(base_atoms_indices)
    basis_vector=calc_basis_vector(*calc_new_base_atoms(coordinates_array, indices))
    gauss_dipole_array[0,0:3]=np.matmul(basis_vector,gauss_dipole_array[0,0:3])
    return pd.DataFrame(gauss_dipole_array,columns=['dipole_x','dipole_y','dipole_z','total']).T


def calc_single_bond_length(coordinates_array, atom_indices):
    """
    a function that gets 2 atom indices, and returns the distance between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates

    atom_indices- list of ints
        a list of atom indices to calculate the distance between- [2,3]

    Returns
    -------
    distance: float
        the bond distance between the atoms
    """
    indices = adjust_indices(atom_indices)
    distance = np.linalg.norm(coordinates_array[indices[0]] - coordinates_array[indices[1]])
    return distance


def calc_bonds_length(coordinates_array, atom_pairs):
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

    bond_list = [calc_single_bond_length(coordinates_array, pair) for pair in atom_pairs]
    index = [('bond_length') + str(pair) for pair in pairs]
    pairs_df = pd.DataFrame(bond_list, index=index)
    return pairs_df


def convert_atoms_indices_to_bonds_length_vector(coordinates_array,atoms_indices): ##for calc_angle_between_atoms
    """

    """
    indices=adjust_indices(atoms_indices)#three atoms-angle four atoms-dihedral
    augmented_indices=[indices[0],indices[1],indices[1],indices[2]]
    if len(indices)==4:
        augmented_indices.extend([indices[2],indices[3]])
    indices_pairs=list(zip(augmented_indices[::2],augmented_indices[1::2]))
    bond_lengths_vector=[(coordinates_array[pair[0]]-coordinates_array[pair[1]]) for pair in indices_pairs]
    return bond_lengths_vector
  


def calc_angle_between_atoms(coordinates_array,atoms_indices): #gets a list of atom indices
    """
    a function that gets 3/4 atom indices, and returns the angle between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    atoms_indices- list of ints
        a list of atom indices to calculate the angle between- [2,3,4]
   
    Returns
    -------
    angle: float
        the bond angle between the atoms
    """
    bonds_list=convert_atoms_indices_to_bonds_length_vector(coordinates_array,atoms_indices)
    print(bonds_list)
    if len(atoms_indices)==3:
        angle=calc_angle(bonds_list[0],bonds_list[1]*(-1))*(180/math.pi)
    else:
        first_cross=np.cross(bonds_list[0],bonds_list[1]*(-1))
        second_cross=np.cross(bonds_list[2]*(-1),bonds_list[1]*(-1)) 
        angle=calc_angle(first_cross,second_cross)*(180/math.pi)
    return angle

def get_angle_df(coordinates_array, atom_indices):
    """
    a function that gets a list of atom indices, and returns a dataframe of angles between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
    
    atom_indices- list of lists of ints
        a list of atom indices to calculate the angle between- [[2,3,4],[2,3,4,5]]
    """
    indices_list=['angle_{}'.format(index) if len(index)==3 else 'dihedral_{}'.format(index) for index in atom_indices]
    angle_list= [calc_angle_between_atoms(coordinates_array,index) for index in atom_indices]
    return pd.DataFrame(angle_list,index=indices_list)



def calc_magnitude_from_coordinates_array(coordinates_array): 
    """
    a function that gets a coordinates array and returns a list of magnitudes.
    help function for get_vibration_info method
    """
    magnitude=[np.linalg.norm(row) for row in coordinates_array]
    return magnitude
      
      
def calc_max_frequency_magnitude(vibration_array, info_df, threshhold=1500):##add option to return ordered_info_df-like dot.prod.info
    """
    a function that gets vibration and info dataframes and returns
    the frequecy and IR with the max magnitude for the vibration.
    splits the coordinates of vibration to 3 coordinates and calculates the magnituede. takes frequencys greter than 1500
    and returns the frequency and IR corresponding to max magnitude.
    
    works in a specific molecule directory
    
    Parameters
    ----------
    vibration_array: np.array
        organized vibration file in array.
    info_df: np.dataframe
        organized info file in dataframe. 
        
        Frequency      IR
   0      20.3253  0.0008
   1      25.3713  0.0023
   2      29.0304  0.0019
   
    Returns
    -------
    dataframe
        max frequency and IR for specific vibration.

    Output:
                                    54
            Frequency         1689.5945
            IR intensity        6.5260
    """
    
    magnitude=calc_magnitude_from_coordinates_array(vibration_array)
    df=info_df.T
    df['magnitude']=magnitude
    outer_finger=(df['Frequency'].astype(float)>threshhold)
    index_max=df[outer_finger]['magnitude'].idxmax()
    return info_df[index_max]



def check_pair_in_bonds(atom_pairs,bonds_df): ##help functions for gen_vibration
    """
    a function that checks that the all given atom pair exists as a bond in the bonds_df
    """
    bonds_list=(bonds_df.astype(int)).values.tolist()   
    bool_check=all([pair in bonds_list for pair in atom_pairs])
    return bool_check



def vibrations_dict_to_list(vibration_dict,vibration_atom_nums):
    """
    a function that gets a list of vibration file names and returns a list of vibration arrays.
    Parameters
    ----------
    vibration_file_names : TYPE
        DESCRIPTION.

    Returns
    -------
    vibration_array_pairs ([0]) : list
    containing the arrays from the vibration file, each variable containing two
    arrays corresponding to the bond pairs 
        .
    vibration_array_list ([1]): list
        containing all vibration array coordinates.

    """
    try:
        vibration_array_list=[vibration_dict['vibration_atom_{}'.format(num)] for num in vibration_atom_nums] 
    except KeyError:
        return print('there is no vibration for those atoms-pick another one')
    vibration_array_pairs=list(zip(vibration_array_list[::2], vibration_array_list[1::2]))
    return vibration_array_pairs,vibration_array_list ##optional- do not split to pairs

def atom_pairs_to_coordinates_vector(coordinates_array,atom_pairs):
    atom_pairs_indices=adjust_indices(atom_pairs)
    coordinates_vector=[(coordinates_array[pair[0]]-coordinates_array[pair[1]]) for pair in atom_pairs_indices] #### vector not distance
    return coordinates_vector

def calc_vibration_dot_product(vibration_array_pairs,coordinates_vector):
    vibration_dot_product=[(abs(np.dot(pair[0],coordinates_vector[i]))+abs(np.dot(pair[1],coordinates_vector[i]))) for i,pair in list(enumerate(vibration_array_pairs))]
    return vibration_dot_product

def calc_vibration_dot_product_from_pairs(atom_pairs,coordinates_array,vibration_dict):
    vibration_atom_nums=[item for sublist in atom_pairs for item in sublist]
    vibration_array_pairs,_=vibrations_dict_to_list(vibration_dict,vibration_atom_nums)
    coordinates_vector=atom_pairs_to_coordinates_vector(coordinates_array,atom_pairs)
    vibration_dot_product=calc_vibration_dot_product(vibration_array_pairs,coordinates_vector)
    return vibration_dot_product

def calc_max_frequency_gen_vibration(vibration_dot_product,info_df,threshhold=500): ## fix so 0 is IR and 1 is frequency
    """
    info_df:
        Frequency      IR
   0      20.3253  0.0008
   1      25.3713  0.0023
   2      29.0304  0.0019
    """
    df_list=[pd.DataFrame(np.vstack([vibration,(info_df)['Frequency']])) for vibration in vibration_dot_product]
    print(df_list[0].T)
    filtered_df_list=[df.T[(df.iloc[1]>threshhold)] for df in df_list]
    index_max=[df[0].idxmax() for  df in filtered_df_list]
    max_frequency_vibration=[pd.DataFrame(list(df[index_max[i]]),index=['IR','Frequency']) for i,df in enumerate(df_list)]
    return pd.concat(max_frequency_vibration)


def vibration_ring_array_list_to_vector(vibration_array_list): ##help function ring_vibration
    vec_sum_1_3_5=(vibration_array_list[0]+vibration_array_list[2]+vibration_array_list[3])
    vec_sum_2_4_6=(vibration_array_list[4]+vibration_array_list[1]+vibration_array_list[5])
    return vec_sum_1_3_5, vec_sum_2_4_6

def get_data_for_ring_vibration(info_df,vibration_array_list,coordinates_vector)->pd.DataFrame:

    """
    data_df:
                0          1         2    ...         135        136        137
    0 prods  -0.000900   0.011600  -0.00410  ...     0.00120     0.0000     0.0000
    1 frequency  20.325300  25.371300  29.03040  ...  3273.31690  3287.5192  3627.5123
    2 calc  0.858618   0.997548   0.22508  ...     0.22508        NaN        NaN


    """
    prods=[np.dot(row_1,row_2) for row_1,row_2 in zip(*vibration_ring_array_list_to_vector(vibration_array_list))]
    _, vibration_array_list=vibration_ring_array_list_to_vector(vibration_array_list)
    calc=[abs(math.sin(calc_angle(row,coordinates_vector))) for row in vibration_array_list]
    data_df=pd.DataFrame(np.vstack([prods,(info_df)['Frequency'],calc]))
    
    return data_df

def get_filter_ring_vibration_df(data_df): ##return filtered_df and check 
    filtered_df=data_df.T[(data_df.iloc[0]!=0)&(abs(data_df.iloc[0])>0.1)][data_df.iloc[1]>500].reset_index() ##randon threshhold to work, no >1500
    if (filtered_df.shape[0]==0):
        filtered_df=data_df.T[(data_df.iloc[0]!=0)][(data_df.iloc[1]>500)&(data_df.iloc[1]<700)].reset_index()
        print('Dot products are lower than 0.1 - returning the default 500 - 700 1/cm')
    
    return filtered_df

def calc_min_max_ring_vibration(filtered_df)->pd.DataFrame:
    max_vibration_frequency=filtered_df.iloc[filtered_df[2].idxmax()][1]
    asin_max=math.asin(filtered_df[2].max())*(180/np.pi)
    min_vibration_frequency=filtered_df.iloc[filtered_df[2].idxmin()][1]
    asin_min=math.asin(filtered_df[2].min())*(180/np.pi)
    df=pd.DataFrame((max_vibration_frequency,asin_max,min_vibration_frequency,asin_min),index=['cross','cross_angle','para','para_angle'])
    return df

def get_filtered_ring_df(info_df,coordinates_array,vibration_dict,ring_atom_indices): ## last ring_vibration
    coordinates_vector=atom_pairs_to_coordinates_vector(coordinates_array,ring_atom_indices)[0]
    vibration_atom_nums=[item for sublist in ring_atom_indices for item in sublist]
    _, vibration_array_list=vibrations_dict_to_list(vibration_dict,vibration_atom_nums)
    data_df=get_data_for_ring_vibration(info_df,vibration_array_list,coordinates_vector)
    filtered_df=get_filter_ring_vibration_df(data_df)
    return filtered_df

######## sterimol functions

def direction_atoms_for_sterimol(bonds_df,base_atoms)->list: #help function for sterinol
    """
    a function that return the base atom indices for coordination transformation according to the bonded atoms.
    you can insert two atom indicess-[1,2] output [1,2,8] or three if the first one repeats-[1,2,1] output [1,2,3]
    """
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

def get_molecule_connections(bonds_df,source,direction):
    graph=ig.Graph.DataFrame(edges=bonds_df,directed=True)
    paths=graph.get_all_simple_paths(v=source,mode='all')
    with_direction=[path for path in paths if (direction in path)]
    return np.unique([path for sublist in with_direction for path in sublist] )


def get_specific_bonded_atoms_df(bonds_df,atom_filter,coordinates_df):
    edited_bonds_df=bonds_df[(bonds_df.isin(atom_filter))].dropna().reset_index(drop=True)
    bonds_array=(np.array(edited_bonds_df)-1).astype(int)
    atom_bonds=np.vstack([(coordinates_df.iloc[bond]['atom'].values) for bond in bonds_array]).reshape(-1,2)
    bonded_atoms_df=(pd.concat([pd.DataFrame(atom_bonds),edited_bonds_df],axis=1))
    bonded_atoms_df.columns=[['atom_1','atom_2','index_1','index_2']]
    return bonded_atoms_df


def remove_atom_bonds(bonded_atoms_df,atom_remove='H'):
    atom_bonds_array=np.array(bonded_atoms_df)
    delete_rows_left=np.where(atom_bonds_array[:,0]==atom_remove)[0] #itterrow [0] is index [1] are the values
    delete_rows_right=np.where(atom_bonds_array[:,1]==atom_remove)[0]
    atoms_to_delete=np.concatenate((delete_rows_left,delete_rows_right))
    new_bonded_atoms_df=bonded_atoms_df.drop((atoms_to_delete),axis=0)
    return new_bonded_atoms_df



def remove_nof_bonds(bonded_atoms_df,bonds_array,coordinates_df):
    bonded_atoms_array=np.array(bonded_atoms_df)
    nof_indices=[]
    for idx, row in enumerate(bonded_atoms_array):
        if (row[0]=='H') & (row[1] in ('N','O','F')):
            nof_indices.append(idx)
        elif (row[1]=='H') & (row[0] in ('N','O','F')):
            nof_indices.append(idx)
    new_bonded_atoms_df=bonded_atoms_df.copy()
    for index in nof_indices:
        pair=bonds_array[index]
        coordinates_array=np.array(coordinates_df[['x','y','z']].astype(float))
        bond_lenght=calc_single_bond_length(coordinates_array, pair)
        if bond_lenght>1.6:
            new_bonded_atoms_df=new_bonded_atoms_df.drop(index,axis=0)
    return new_bonded_atoms_df

            
def filter_atoms_for_sterimol(bonded_atoms_df,coordinates_df):### change to filter both sides
    bonds_array=adjust_indices((bonded_atoms_df[['index_1','index_2']]).astype(int))
    bonded_atoms_df=remove_nof_bonds(bonded_atoms_df, bonds_array, coordinates_df)
    bonded_atoms_df=remove_atom_bonds(bonded_atoms_df) ## a bit redundent becuse remove_atom_bonds doesnt leave nof bonds anyways
    allowed_bonds_indices= pd.concat([bonded_atoms_df['index_1'],bonded_atoms_df['index_2']],axis=1).reset_index(drop=True)
    atom_filter=adjust_indices(np.unique([atom for sublist in allowed_bonds_indices.values.tolist() for atom in sublist]))
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
    coordinates_df['magnitude']=calc_magnitude_from_coordinates_array(np.array(coordinates_df[['x','z']].astype(float)))
    if radii=='bondi':
        coordinates_df['radius']=coordinates_df['atom'].replace(GeneralConstants.BONDI_RADII.value)
    elif radii=='CPK':
        coordinates_df['radius']=coordinates_df['atype'].replace(GeneralConstants.CPK_RADII.value)
    coordinates_df['B5']=coordinates_df[['radius','magnitude']].sum(axis=1)
    coordinates_df['L']=coordinates_df[['y','radius']].sum(axis=1)
    return coordinates_df


def get_transfomed_plane_for_sterimol(plane,degree):
    """
    a function that gets a plane and rotates it by a given degree
    in the case of sterimol the plane is the x,z plane.
    Parameters:
    ----------
    plane : np.array
        [x,z] plane of the molecule coordinates.
        example:
            [-0.6868 -0.4964]
    degree : float
    """
    cos_deg=np.cos(degree*(np.pi/180))
    sin_deg=np.sin(degree*(np.pi/180))
    rot_matrix=np.array([[cos_deg,-1*sin_deg],[sin_deg,cos_deg]])
    transformed_plane=np.vstack([np.matmul(rot_matrix,row) for row in plane]).round(4)
    return transformed_plane


def calc_B1(transformed_plane,avs,edited_coordinates_df,column_index):
    """
    Parameters
    ----------
    transformed_plane : np.array
        [x,z] plane of the molecule coordinates.
        example:
            [-0.6868 -0.4964]
            [-0.7384 -0.5135]
            [-0.3759 -0.271 ]
            [-1.1046 -0.8966]
            [ 0.6763  0.5885]
    avs : list
        the max & min of the [x,z] columns from the transformed_plane.
        example:[0.6763, -1.1046, 0.5885, -0.8966
                 ]
    edited_coordinates_df : TYPE
        DESCRIPTION.
    column_index : int
        0 or 1 depending- being used for transformed plane.
    """
    
    ## get the index of the min value in the column compared to the avs.min
    idx=np.where(np.isclose(np.abs(transformed_plane[:,column_index]),(avs.min()).round(4)))[0][0]
    if transformed_plane[idx,column_index]<0:
        new_idx=np.where(np.isclose(transformed_plane[:,column_index],transformed_plane[:,column_index].min()))[0][0]
        bool_list=np.logical_and(transformed_plane[:,column_index]>=transformed_plane[new_idx,column_index],
                                 transformed_plane[:,column_index]<=transformed_plane[new_idx,column_index]+1)
        
        transformed_plane[:,column_index]=-transformed_plane[:,column_index]
    else:
        bool_list=np.logical_and(transformed_plane[:,column_index]>=transformed_plane[idx,column_index]-1,
                                 transformed_plane[:,column_index]<=transformed_plane[idx,column_index])
        
    against,against_loc=[],[]
    B1,B1_loc=[],[]
    for i in range(1,transformed_plane.shape[0]): 
        if bool_list[i]:
            against.append(np.array(transformed_plane[i,column_index]+edited_coordinates_df['radius'].iloc[i]))
            against_loc.append(edited_coordinates_df['L'].iloc[i])
        if len(against)>0:
            B1.append(max(against))
            B1_loc.append(against_loc[against.index(max(against))])
        else:
            B1.append(np.abs(transformed_plane[idx,column_index]+edited_coordinates_df['radius'].iloc[idx]))
            B1_loc.append(edited_coordinates_df['radius'].iloc[idx])
            
    return [B1,B1_loc]

def b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane):
    """
    a function that gets a plane transform it and calculate the b1s for each degree.
    checks if the plane is in the x or z axis and calculates the b1s accordingly.
    Parameters:
    ----------
    extended_df : pd.DataFrame
    b1s : list
    b1s_loc : list
    degree_list : list
    plane : np.array
    """
    for degree in degree_list:
        transformed_plane=get_transfomed_plane_for_sterimol(plane, degree)
        avs=np.abs([max(transformed_plane[:,0]),min(transformed_plane[:,0]), 
                    max(transformed_plane[:,1]),min(transformed_plane[:,1])])
        
        if np.where(avs==avs.min())[0][0] in [0,1]:
                B1,B1_loc=calc_B1(transformed_plane,avs,extended_df,0)
  
        elif np.where(avs==avs.min())[0][0] in [2,3]:
            B1,B1_loc=calc_B1(transformed_plane,avs,extended_df,1)
            
        b1s.append(np.unique(np.vstack(B1)).max())####check
        b1s_loc.append(np.unique(np.vstack(B1_loc)).max())

def get_b1s_list(extended_df, scans=90//5):
    b1s,b1s_loc=[],[]
    scans=scans
    degree_list=list(range(18,108,scans))
    plane=np.array(extended_df[['x','z']].astype(float))
    b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane)
               
    back_ang=degree_list[np.where(b1s==min(b1s))[0][0]]-scans   
    front_ang=degree_list[np.where(b1s==min(b1s))[0][0]]+scans                      
    
    degree_list=range(back_ang,front_ang+1)
    b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane)

    return [np.array(b1s),np.array(b1s_loc)]

def calc_sterimol(bonded_atoms_df,extended_df):
    edited_coordinates_df=filter_atoms_for_sterimol(bonded_atoms_df,extended_df)
    b1s,b1s_loc=get_b1s_list(edited_coordinates_df)
    B1=min(b1s[b1s>=0])
    loc_B1=max(b1s_loc[np.where(b1s[b1s>=0]==min(b1s[b1s>=0]))])
    B5=max(edited_coordinates_df['B5'].values)
    L=max(edited_coordinates_df['L'].values+0.4)
    loc_B5=min(edited_coordinates_df['y'].iloc[np.where(edited_coordinates_df['B5'].values==B5)[0]])
    sterimol_df=pd.DataFrame([B1,B5,L,loc_B5,loc_B1],index=['B1','B5','L','loc_B5','loc_B1']) ## loc_B1 produces different values
    return sterimol_df



class Molecule():
    """
    """
    def __init__(self,molecule_log_filename):
        """
        

        Parameters:
        ----------
        molecule_dir_name : str
            the name of the molecule directory.

        Description:
        -------
        

        """
       
        self.molecule_name=molecule_log_filename.split('.')[0]
        self.molecule_path=os.path.abspath(self.molecule_name)
        # os.chdir(self.molecule_path)
        self.parameter_list=gauss_file_handler(molecule_log_filename)
        self.xyz_df=self.parameter_list[0]['standard_orientation_df']
        self.coordinates_array=np.array(self.xyz_df[['x','y','z']].astype(float))
        self.gauss_dipole_df=self.parameter_list[0]['dipole_df']
        self.polarizability_df=self.parameter_list[0]['pol_df']
        self.atype_df=self.parameter_list[0]['atype_df']
        self.bonds_df=self.parameter_list[0]['bonds_df']
        self.info_df=self.parameter_list[0]['info_df']
        self.vibration_dict=self.parameter_list[1]
        # self.pbmol=self.parameter_list[7]
        # os.chdir('../')
     
    
    
    
    def get_sterimol(self,base_atoms,radii='bondi'):
        
        bonds_direction=direction_atoms_for_sterimol(self.bonds_df,base_atoms)
        transformed_coordinates=self.get_coordination_transformation_df(bonds_direction)
        connected_from_direction=get_molecule_connections(self.bonds_df,bonds_direction[0],bonds_direction[1])
        extended_df=get_extended_df_for_sterimol(transformed_coordinates,self.atype_df,radii)
        bonded_atoms_df=get_specific_bonded_atoms_df(self.bonds_df,connected_from_direction,self.xyz_df)     
        sterimol_df=calc_sterimol(bonded_atoms_df,extended_df)
        return sterimol_df
        

    def swap_atom_pair(self,pair_indices): #swapping is permanently changed
        pairs=adjust_indices(pair_indices)
        xyz_df=self.xyz_df
        temp=xyz_df.iloc[pairs[0]].copy()
        xyz_df.iloc[pairs[0]]=self.coordinates_array[pairs[1]]
        xyz_df.iloc[pairs[1]]=temp
        return xyz_df
        
    def get_coordination_transformation_df(self,base_atoms_indices):
        transformed_coordinates=calc_coordinates_transformation(self.coordinates_array,base_atoms_indices)
        new_coordinates_array=(np.column_stack([np.array(self.xyz_df['atom']),transformed_coordinates]))
        new_coordinates_df=pd.DataFrame(new_coordinates_array,columns=help_functions.XYZConstants.DF_COLUMNS.value)
        return new_coordinates_df
    
    ###dont have charges yet
    # def get_npa_df(self,base_atoms_indices,sub_atoms=None):
    #     charges=np.array(self.get_unedited_df('npa').astype(float))
    #     coordinates_array=np.array(self.get_coordination_transformation_df(base_atoms_indices)[['x','y','z']].astype(float))
    #     dipole_df=calc_npa_charges(coordinates_array,charges,base_atoms_indices,sub_atoms)
        
    #     return dipole_df
    
    
    def get_dipole_gaussian_df(self,base_atoms_indices):
        return calc_dipole_gaussian(self.coordinates_array,np.array(self.gauss_dipole_df),base_atoms_indices)       
    
    def get_bond_angle(self,atom_indices):##need expand to many angles
        return get_angle_df(self.coordinates_array,atom_indices)
    
    def get_bond_length(self,atom_pairs):
        if type(atom_pairs)==list:
            bond_length=calc_single_bond_length(self.coordinates_array,atom_pairs)
        elif type(atom_pairs)==tuple:
            bond_df=calc_bonds_length(self.coordinates_array,atom_pairs)
        return bond_df
    

    ########## modified till here
    def get_vibration_max_frequency(self,vibration_atom_num):
        vibration_array=self.vibration_dict.get('vibration_atom_{}'.format(vibration_atom_num))
        return calc_max_frequency_magnitude(vibration_array,self.info_df.T)
    
    
    def get_gen_vibration(self,atom_pairs):
        """

        Parameters
        ----------
        atom_pairs : molecule_1.get_gen_vibration([[1,6],[3,4]])
            atom pairs must have a corresponding vibration file and appear in bonds_df.

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
        if check_pair_in_bonds(atom_pairs,self.bonds_df)==True:
  
            try:
                vibration_dot_product=calc_vibration_dot_product_from_pairs(atom_pairs,self.coordinates_array,self.vibration_dict)
            except TypeError:
                return print('Error: no vibration file for the molecule')
            return calc_max_frequency_gen_vibration(vibration_dot_product,self.info_df)
        else:
            print('Error: the following bonds do not exist-check atom numbering')
            
    def get_ring_vibrations(self,ring_atom_indices):
        """

        Parameters
        ----------
        ring_atom_indices :working example: molecule_1.get_ring_vibrations([[1,3],[1,6],[3,4]]) 
            a list of atom pairs, there must be a vibration file with a corresponding number to work.
            
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
        try:
            filtered_df=get_filtered_ring_df(self.info_df,self.coordinates_array,self.vibration_dict,ring_atom_indices)
        except FileNotFoundError:
            return print('no vibration data file exists for this atom')
        # bool_check=filtered_df[2].duplicated().any() #check for duplicates in calc, is it needed ?? slim chances
        return calc_min_max_ring_vibration(filtered_df)
    
    def check_imaginary_frequency(self):##return True if no complex frequency, called ground.state in R
        bool_imaginary=not any([isinstance(frequency, complex) for frequency in self.info_df['Frequency']])
        return bool_imaginary
    
    
    # def get_nbo_info(self,atom_pairs):
    #     pairs_indices=np.array(atom_pairs)-1
    #     index=[('nbo charge')+str(pair) for pair in pairs_indices]
    #     nbo_array= np.array((self.nbo_df.drop([0,1],axis=1)).astype(float))
    #     nbo_diff=[abs(statistics.mean(nbo_array[pair[0]])-statistics.mean(nbo_array[pair[1]])) for pair in pairs_indices]
    #     return pd.DataFrame(nbo_diff, index=index)
    
    def get_molecule_comp_set_hetro(self,dipole_mode = 'gaussian', radii = 'bondi'):
        """
      Ring atoms - by order -> primary axis (para first), ortho atoms and meta atoms: 19 20 19 21 20 22
      your atom pairs: 19 20 19 21
      Enter atoms - origin atom, y axis atom and xy plane atom: 1 4 5
      Primary axis along: 2 3
      Distances - Atom pairs: 2 3 2 5
      __main__:192: RuntimeWarning: invalid value encountered in double_scalars
      __main__:483: UserWarning: Boolean Series key will be reindexed to match DataFrame index.
      __main__:485: UserWarning: Boolean Series key will be reindexed to match DataFrame index.
      Dot products are lower than 0.1 - returning the default 500 - 700 1/cm
      Do you want to compute any angles/dihedrals? y/n: n
        -------
        res : dataframe
        	molecule1
        cross	657.3882
        cross_angle	81.17206344670865
        para	834.4249
        para_angle	40.67483293898943
        dipole_x	0.058148468892639825
        dipole_y	1.210906929546481
        dipole_z	0.9491739479893856
        total	1.5397
        nbo charge[0 1]	0.2706980000000003
        nbo charge[2 3]	0.013509999999999689
        B1	2.556
        B5	4.131529971026473
        L	5.8357
        loc_B5	1.4121
        loc_B1	3.6855
        bond length[0 1]	1.4566677726921815
        bond length[2 3]	1.2277597484850202
        angle_[1, 2, 3]	69.7212447900024
        dihedral_[1, 2, 3, 4]	104.74132137986373


        """
        res=[]
        ring_atom_indices=help_functions.split_to_pairs(input("Ring atoms - by order -> primary axis (para first), ortho atoms and meta atoms: "))
        gen_vibration_indices=help_functions.split_to_pairs(input("your atom pairs: "))
        dipole_indices=help_functions.string_to_int_list(input('Enter atoms - origin atom, y axis atom and xy plane atom: '))
        # nbo_indices=help_functions.split_to_pairs(input('Insert atom pairs for which you wish calculate differences: '))
        sterimol_indices=help_functions.string_to_int_list(input('Primary axis along: '))
        indices_for_dist=help_functions.split_to_pairs(input('Distances - Atom pairs: '))
        res.append(self.get_ring_vibrations(ring_atom_indices))
        res.append(self.get_gen_vibration(gen_vibration_indices))
        if dipole_mode=='compute':
            res.append(self.get_npa_df(dipole_indices))
        elif dipole_mode=='gaussian':
            res.append(self.get_dipole_gaussian_df(dipole_indices))
        # res.append(self.get_nbo_info(nbo_indices))
        if radii=='CPK':
            res.append(self.get_sterimol(sterimol_indices,radii))
        else:
           res.append(self.get_sterimol(sterimol_indices) )
        res.append(self.get_bond_length(indices_for_dist))
        
        angle_answer=input('Do you want to compute any angles/dihedrals? y/n: ')
        if angle_answer=='y':
            print('Insert a list of atom triads/quartets for which you wish to have angles/dihedrals.\n'
        'For several angles/dihedrals, insert triads/quartets with a double space between them.\n'
        'Make sure there are spaces between atoms as well.\n'
        'For example - 1 2 3  1 2 3 4 will give angle(1, 2, 3) and dihedral(1, 2, 3, 4)')
            angle_string=help_functions.split_for_angles(input('Insert a list of atom triads/quartets for which you wish to have angles/dihedrals: '))
            res.append(self.get_bond_angle(angle_string))
        res.append(self.polarizability_df)
        
        df=pd.concat(res).rename(columns={0:self.molecule_name}).T
        df.to_csv('comp_set.csv')
        print(df)
        return df
            
    
class Molecules():
    def __init__(self,molecules_dir_name):
        self.molecules_path=os.path.abspath(molecules_dir_name)
        os.chdir(self.molecules_path) 
        self.molecules=[Molecule(log_file) for log_file in os.listdir() ] 
        os.chdir('../')
        
    def get_molecules_comp_set(self,dipole_mode = 'gaussian', radii = 'bondi'):
        """
        molecules.get_molecules_comp_set()
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
        df : 
            
        """
        df=pd.DataFrame()
        ring_atom_indices=help_functions.split_to_pairs(input("Ring atoms - by order -> primary axis (para first), ortho atoms and meta atoms: "))
        gen_vibration_indices=help_functions.split_to_pairs(input("your atom pairs: "))
        dipole_indices=help_functions.string_to_int_list(input('Enter atoms - origin atom, y axis atom and xy plane atom: '))
        nbo_indices=help_functions.split_to_pairs(input('Insert atom pairs for which you wish calculate differences: '))
        sterimol_indices=help_functions.string_to_int_list(input('Primary axis along: '))
        indices_for_dist=help_functions.split_to_pairs(input('Distances - Atom pairs: '))
        angle_answer=input('Do you want to compute any angles/dihedrals? y/n: ')
        if angle_answer=='y':
            print('Insert a list of atom triads/quartets for which you wish to have angles/dihedrals.\n'
        'For several angles/dihedrals, insert triads/quartets with a double space between them.\n'
        'Make sure there are spaces between atoms as well.\n'
        'For example - 1 2 3  1 2 3 4 will give angle(1, 2, 3) and dihedral(1, 2, 3, 4)')
            angle_string=help_functions.split_for_angles(input('Insert a list of atom triads/quartets for which you wish to have angles/dihedrals: '))
            
        for molecule in self.molecules :
            res=[]
            res.append(molecule.get_ring_vibrations(ring_atom_indices))
            res.append(molecule.get_gen_vibration(gen_vibration_indices))
            if dipole_mode=='compute':
                res.append(molecule.get_npa_df(dipole_indices))
            elif dipole_mode=='gaussian':
                res.append(molecule.get_dipole_gaussian_df(dipole_indices))
            res.append(molecule.get_nbo_info(nbo_indices))
            if radii=='CPK':
                res.append(molecule.get_sterimol(sterimol_indices,radii))
            else:
               res.append(molecule.get_sterimol(sterimol_indices) )
            res.append(molecule.get_bond_length(indices_for_dist))
            if angle_answer=='y':
                res.append(molecule.get_bond_angle(angle_string))
            res.append(molecule.get_pol_info())
            df[molecule.molecule_name]=(pd.concat(res))
        os.chdir('../')
        (df.T).to_csv('comp_set.csv')
        print(df.T)
        return df.T



if __name__=='__main__':
    # xyz_file_generator_library(r'C:\Users\edens\Documents\GitHub\learning_python\project\main_python','new_directory') #works
    path=r'C:\Users\edens\Documents\GitHub\main_python\gauss_log'
    os.chdir(path)
    molecules=Molecule('AB_para_methoxy.log')
    # molecules.get_molecule_comp_set_hetro()
 
    
 
    
 
    