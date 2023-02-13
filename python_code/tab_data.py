"""
This moudle provides analysis for xyz files, such as tab data (like distances, bonds, etc.) and overlay data (like overlay graph of 2 or more xyz files)
This basis of this code was taken from Github and rewritten for efficency and production runs (since it worked only from command line)

Source github: https://github.com/radi0sus/xyzoverlay (xyz overlay)
               https://github.com/radi0sus/xyz2tab (tab data)

Kabsch algorithm: Kabsch W., 1976, A solution for the best rotation to relate two sets of vectors,
                  Acta Crystallographica, A32:922-923, doi: http://dx.doi.org/10.1107/S0567739476001873

"""

import sys                                                      #stdout

import os.path                                                  #filename split
import argparse                                                 #argument parser
import re                                                       #regex
import pandas as pd                                             #pandas tables
import numpy as np                                              #for calculations
from scipy.spatial.distance import pdist, squareform, cosine    #for the calculations of the distance matrix and angles (cosine)
##from tabulate import tabulate                                   #nice table output
from itertools import combinations, cycle                       #for r-length tuples, in sorted order, no repeated elements
import matplotlib.pyplot as plt                                 #for molecule display
from mpl_toolkits.mplot3d import Axes3D                         #for molecule display
from matplotlib.patches import FancyArrowPatch                  #for fancy arrows in xyz
from mpl_toolkits.mplot3d import proj3d                       	#for fancy arrows in xyz
from cycler import cycler                                       #generate color cycle
import io                                                       #IO for (easy) saving multi xyz
from enum import Enum

path_to_add=r'C:\Users\edens\Documents\GitHub\Crystal_structure'
# os.chdir(path_to_add)
sys.path.insert(0, path_to_add)

from tools.general_constants import *
from tools.file_handlers import *
from tools.rmsd_wrappers import *
from openbabel import pybel
from openbabel import openbabel as ob
############
#Constansts#
############

class TABConstants(Enum):
    """
    Holds constants for tab constants
    1. covalent radii from Alvarez (2008) DOI: 10.1039/b801115j
    2. atomic weights
    3. dict for numbers to subscript numbers
    """
    covalent_radii = {
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
    
    atomic_weights = {
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

    utf_sub_dict = {
        "0" : "₀",
        "1" : "₁",
        "2" : "₂",
        "3" : "₃",
        "4" : "₄",
        "5" : "₅",
        "6" : "₆",
        "7" : "₇",
        "8" : "₈",
        "9" : "₉",
    }

class PlotConstants(Enum):
    """
    Holds constants for plot preference (tab and overlay)
    1. lists for color asignment
    2. plot colors preferences
    """
    METAL_ATOMS=['Ac', 'Ag', 'Al', 'Am', 'Ar', 'As', 'At', 'Au', 'Ba', 'Be',
                 'Bi', 'Ca', 'Cd', 'Ce', 'Cf', 'Cm', 'Co', 'Cr', 'Cs', 'Cu',
                 'Db', 'Dy', 'Er', 'Es', 'Eu', 'Fe', 'Fm', 'Fr', 'Ga', 'Gd',
                 'Ge', 'Hf', 'Hg', 'Ho', 'Hs', 'In', 'Ir', 'K', 'La', 'Li',
                 'Lr', 'Lu', 'Md', 'Mg', 'Mn', 'Mo', 'Na', 'Nb', 'Nd', 'Ni',
                 'Np', 'Os', 'Pa', 'Pb', 'Pd', 'Pm', 'Po', 'Pr', 'Pt', 'Pu',
                 'Ra', 'Rb', 'Re', 'Rf', 'Rh', 'Rn', 'Ru', 'Sc', 'Sm', 'Sn',
                 'Sr', 'Ta', 'Tb', 'Tc', 'Te', 'Th', 'Ti', 'Tl', 'Tm', 'U',
                 'V', 'W', 'Y', 'Yb', 'Zn', 'Zr']
    GREEN_ATOMS=['F','Cl']
    BROWN_ATOMS=['Br']
    PURPLE_ATOMS=['P','I']
    ORANGE_ATOMS=['Si']
    OUTPUT_DPI=150 #dpi for figure
    ALPHA_ATOMS=0.55 #alpha atoms
    ALPHA_BONDS=0.55 #alpha bonds
    ATOM_SCALER=1e60 #sphere radius for atom view, change exponent
    BOND_SCALER=1e6 #cylinder radius for bond view, change exponent

class PlotMapping(Enum):
    """
    WRITE ME
    """
    color_grouping={'C': {'elements': ['C'], 'size': 100, 'color': 'black'},
                    'H': {'elements': ['H'], 'size': 100, 'color': 'tan'},
                    'N': {'elements': ['N'], 'size': 100, 'color': 'blue'},
                    'O': {'elements': ['O'], 'size': 100, 'color': 'green'},
                    'S': {'elements': ['S'], 'size': 100, 'color': 'yellow'},
                    'Metals': {'elements': PlotConstants.METAL_ATOMS.value, 'size': 100, 'color': 'red', 'alpha': 0.85},
                    'green_atoms': {'elements': PlotConstants.GREEN_ATOMS.value, 'size': 100, 'color': 'green'},
                    'brown_atoms': {'elements': PlotConstants.BROWN_ATOMS.value, 'size': 100, 'color': 'brown'},
                    'purple_atoms': {'elements': PlotConstants.PURPLE_ATOMS.value, 'size': 100, 'color': 'purple'},
                    'orange_atoms': {'elements': PlotConstants.ORANGE_ATOMS.value, 'size': 100, 'color': 'orange'},
                    'leftover_atoms': {'size': 100, 'color': 'gray'},
                    }   
##        COLOR_CYCLE=cycler(c=['b', 'r', 'g', 'c', 'm', 'y']) #default color cycle

class TabData():
    """
    A void class to help arrange final data for production run
    """
    pass

def set_tab_parser():
    """
    A function that setup a Parser for tab data
    Parser is a commpand port interface that connects user input to commands
    """
    parser = argparse.ArgumentParser(
                    prog='xyz2tab', 
                    description = "Print bond, lengths angles and more from xyz files.")

    #filename is required
    parser.add_argument("filename",
            help = "filename, xyz; e.g. mymolecule.xyz")

    #exclude atoms
    parser.add_argument('-ea','--excludeAt',
            nargs="+",
            type=str,
            help='exclude bonds and angles to specified atoms; e.g. -ea N1 or -ea N1 N2')

    #exclude elements
    parser.add_argument('-ee','--excludeEl',
            nargs="+",
            type=str,
            help='exclude bonds and angles to specified elements; e.g. -ee C or -ee C N')

    #sort by value
    parser.add_argument('-sa','--sortasc',
            default=0,
            action='store_true',
            help='sort values for bond lengths and angles ascending')

    #sort by value
    parser.add_argument('-sd','--sortdes',
            default=0, 
            action='store_true',
            help='sort values for bond lengths and angles descending')

    #sort by name
    parser.add_argument('-sae','--sortascEl',
            default=0, 
            action='store_true',
            help='ascending alphabetical sort of elements')

    #sort by name
    parser.add_argument('-sde','--sortdesEl',
            default=0, 
            action='store_true',
            help='descending alphabetical sort of elements')

    #include contacts
    parser.add_argument('-ic','--includeCon',
            nargs="+",
            type=str,
            help='include contacts, e.g. -ic C0 C11 or -ic C0 C11 N14')


    #calculate dihedral angle of selected atoms
    parser.add_argument('-d','--dihedral',
            nargs=4,
            type=str,
            help='calculate the dihedral angle of 4 atoms, e.g. -d C0 C11 C12 N13')
            
    #calculate the best plane 1
    parser.add_argument('-p1','--plane1',
            nargs='+',
            type=str,
            help='calculate the best plane through selected (A B C), a range of (A : C) or all atoms (:)\
                    and show the distances to plane 1, \
                    e.g. -p1 C0 C11 C12 N13 or -p1 C0 : N13 or -p1 :')

    #calculate the best plane 2
    parser.add_argument('-p2','--plane2',
            nargs='+',
            type=str,
            help='calculate the best plane through selected (A B C), a range of (A : C) or all atoms (:)\
                    and show the distances to plane 1 and 2 and the angle to plane 1, \
                    e.g. -p2 C21 C22 C23 N23 or -p2 C21 : N23 or -p2 :')

    #add +x% to radius
    parser.add_argument('-r','--radius',
            default=8,
            type=float,
            help='enlarge atomic radii by x %%, e.g. -r 15.2, default is 8 %%')

    #verbose
    parser.add_argument('-v','--verbose',
            default=0, action='store_true',
            help='verbose print, includes 2 more tables')

    #index
    parser.add_argument('-i','--index',
            default=0, action='store_true',
            help='the index for the first atom is 1 instead of 0')

    #plot
    parser.add_argument('-s','--show',
            default=0, action='store_true',
            help='plot xyz coordinates, bonds and planes')

    #plot with bond lengths
    parser.add_argument('-sb','--showbl',
            default=0, action='store_true',
            help='same as -s with bond lengths')

    #plot with no labels
    parser.add_argument('-sn','--shownl',
            default=0, action='store_true',
            help='same as -s with no labels')

    #show orientation
    parser.add_argument('-so','--showori',
            default=0, action='store_true',
            help='plot three arrows along the xyz-axes at 0,0,0 to show the orientation of the molecule')
    return parser

def set_overlay_parser():
    """
    A function that setup a Parser for xyz overlay 
    Parser is a commpand port interface that connects user input to commands
    """
    #argument parser
    parser = argparse.ArgumentParser(
                    prog='xyzoverlay', 
                    description = "Overlay and align two or more molecules from (multi) xyz files.")

    #filename is required
    parser.add_argument("filename", 
            nargs = '+',
            help = "filename(s), xyz; e.g. mol1.xyz or mol1.xyz mol2.xyz")

    #atoms for superposition
    parser.add_argument('-a', '--atoms',
            nargs= '+',
            type= int,
            action= 'append',
            help =  "select atoms for superposition, e.g. -a 1 2 -a 3 4")

    #select the same atoms in all xyz for superposition
    parser.add_argument('-sa', '--sameatoms',
            nargs= '+',
            type= int,
            help =  "select the same atoms in all xyz for superposition, e.g. -sa 1 2 3 ")

    #select all atoms for superposition
    parser.add_argument('-aa','--allatoms',
            default=0, 
            action='store_true',
            help='select all atoms for superposition, xyz files must have the same number of atoms')

    #plot or view, color by atom
    parser.add_argument('-vca','--viewcba',
            default=0, action='store_true',
            help='view xyz coordinates (atoms) and bonds, color by atom')

    #plot or view, color by molecule
    parser.add_argument('-vcm','--viewcbm',
            default=0, action='store_true',
            help='view xyz coordinates (atoms) and bonds, color by molecule')

    #use Matplotlib colormap
    parser.add_argument('-cm','--cmap',
            type= str,
            help='use Matplotlib colormap in view mode')

    #add +x% to radius (for view mode)
    parser.add_argument('-r','--radius',
            default=8,
            type=float,
            help='enlarge atomic radii by x %%, e.g. -r 15.2, default is 8 %%')

    #exclude elements in view mode
    parser.add_argument('-ee','--excludeEl',
            nargs="+",
            type=str,
            help='exclude specified elements in view mode; e.g. -ee C or -ee C N')

    #exclude atoms in view mode
    parser.add_argument('-ea','--excludeAt',
            nargs="+",
            type=int,
            help='exclude specified atoms in view mode; e.g. -ea 1 or -ea 1 2')

    #save single xyz
    parser.add_argument('-s','--save',
            default=0, action='store_true',
            help='save superimposed / aligned data as single xyz file(s)')

    #save multi xyz
    parser.add_argument('-st','--savetr',
            default=0, action='store_true',
            help='save superimposed / aligned data as (multi) xyz or xyz trajectory file')

    return parser

## General functions
def atom_idx_to_element(atom_idx):
    """
    A function that takes a string of atom index like (C11 or Si13, etc.),
    and returns only the element symbol from the index
    ----------
    Parameters
    ----------
    atom_idx : str.
    Atom index in the form of "element"+"index" like C11 or Si13
    -------
    Returns
    -------
    string
    The element symbol from the index 
    --------
    Examples
    --------
    atom_idx='C11'
    element_symbol=atom_idx_to_element(atom_idx)
    print(element_symbol)
            'C'
    """
    return re.sub(r'\d+', '', atom_idx)


def atom_idx_to_index(atom_idx):
    """
    A function that takes a string of atom index like (C11 or Si13, etc.),
    and returns only the number from the index
    ----------
    Parameters
    ----------
    atom_idx : str.
    Atom index in the form of "element"+"index" like C11 or Si13
    -------
    Returns
    -------
    integer
    The number from the index 
    --------
    Examples
    --------
    atom_idx='C11'
    only_idx=atom_idx_to_index(atom_idx)
    print(only_symbol)
            11
    """
    return int(re.sub(r'\D+', '', atom_idx))


def num_to_subnum(number):
    """
    A function that takes a number or string and return the subscript (utf8) version of the number
    Used to represent molecular formula
    ----------
    Parameters
    ----------
    number : str or int.
    The number you want to convert into subscript 
    -------
    Returns
    -------
    tuple
    The converted number as string in the first(and only) position in the tuple 
    --------
    Examples
    --------
    number='123'
    subnum=atom_idx_to_index(atom_idx)
    print(subnum)
            ('₁₂₃')
    """
    utf_number=''
    for letter in str(number):
            utf_letter=TABConstants.utf_sub_dict.value[letter]
            utf_number=utf_number+utf_letter
    return (utf_number)


def calc_angle(xyzarr, i, j, k):
    """
    A function that takes array of all atomic coordinates (xyzarr) and 3 indices (i, j , k)
    and calculate angle from 3 vectors / atomic coordinates : i(x,y,z); j(x,y,z); k(x,y,z)
    ----------
    Parameters
    ----------
    xyzarr : np.array.
    Array of cartesian coordinates of atoms in the molecule

    i : int.
    The index of the 1st participating atom

    j : int.
    The index of the 2nd participating atom

    k : int.
    The index of the 3rd participating atom
    -------
    Returns
    -------
    float
    The angle between the participating atoms in degrees 
    --------
    Examples
    --------
    WRITE ME
    """
    rij = xyzarr[i] - xyzarr[j]
    rkj = xyzarr[k] - xyzarr[j]
    cosine_theta = cosine(rij,rkj)
    theta = np.arccos(1-cosine_theta)
    return 	np.degrees(theta)


def calc_d_angle(xyzarr, i, j, k, l):
    """
    A function that takes array of all atomic coordinates (xyzarr) and 4 indices (i, j , k, l)
    and calculate dihedral angle from 4 vectors / atomic coordinates : i(x,y,z); j(x,y,z); k(x,y,z) ; l(x,y,z)
    ----------
    Parameters
    ----------
    xyzarr : np.array.
    Array of cartesian coordinates of atoms in the molecule

    i : int.
    The index of the 1st participating atom

    j : int.
    The index of the 2nd participating atom

    k : int.
    The index of the 3rd participating atom

    l : int.
    The index of the 4th participating atom
    -------
    Returns
    -------
    float
    The angle between the participating atoms in degrees 
    --------
    Examples
    --------
    WRITE ME
    """
    np.seterr(invalid='ignore') #no warning if division by zero
    rji = -1*(xyzarr[j] - xyzarr[i])
    rkj = xyzarr[k] - xyzarr[j]
    rlk = xyzarr[l] - xyzarr[k]
    rkj /= np.linalg.norm(rkj)
    v = rji - np.dot(rji, rkj)*rkj
    w = rlk - np.dot(rlk, rkj)*rkj
    x = np.dot(v, w)
    y = np.dot(np.cross(rkj, v), w)
    return 	np.degrees(np.arctan2(y,x))


def svd_fit(X):
    """
    calculation of the best-fit plane (Singular value decomposition
    https://gist.github.com/bdrown/a2bc1da0123b142916c2f343a20784b4
    ----------
    Parameters
    ----------
    X : ??.
    Array of cartesian coordinates of atoms in the molecule ??
    -------
    Returns
    -------
    ???
    ???
    --------
    Examples
    --------
    WRITE ME
    """
    C = np.average(X, axis=0)
    CX = X - C # Create CX vector (centroid to point) matrix
    U, S, V = np.linalg.svd(CX) # Singular value decomposition
    N = V[-1] # The last row of V matrix indicate the eigenvectors of smallest eigenvalues (singular values).
    return C, N


def kabsch(P, Q):
    """
    Kabsch algorithm
    https://github.com/charnley/rmsd
    ----------
    Parameters
    ----------
    P : ??.
    Array of cartesian coordinates of atoms in the molecule ??

    Q : ??.
    Array of cartesian coordinates of atoms in the molecule ??
    -------
    Returns
    -------
    ???
    ???
    --------
    Examples
    --------
    WRITE ME
    """
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]
    U = np.dot(V, W)
    return U # create Rotation matrix U

def align_xyz(vec1,vec2,coord):
    """
    get rotmatrix from Kabsch and transform coordinates
    ----------
    Parameters
    ----------
    vec1 : ??.
    Array of cartesian coordinates of atoms in the molecule ??

    vec2 : ??.
    Array of cartesian coordinates of atoms in the molecule ??

    coord : ??.
    Array of cartesian coordinates of atoms in the molecule ??
    -------
    Returns
    -------
    ???
    ???
    --------
    Examples
    --------
    WRITE ME        
    """
    rotmatrix = kabsch(vec1,vec2)
    return np.dot(coord,rotmatrix)

## plot functions
def set_plot_settings():
    """
    A void function that concentrates all plot viewing setting defaults
    """
    plt.rcParams['savefig.dpi'] = PlotConstants.OUTPUT_DPI.value

def init_figure_and_axis():
    """
    WRITE ME
    """
    #prepare plot
    fig = plt.figure(figsize=(10,8))
    ax = plt.axes(projection='3d')
    #otherwise molecule looks strange
    ax.set_box_aspect((1, 1, 1))
    return fig, ax

def finallize_figure_and_axis(fig, ax):
    """
    WRITE ME
    """
    #no axes
    ax.set_axis_off()
    #tight layout 
    fig.tight_layout()
    #adjust 3d drawing behavior, otherwise molecules are not correctly displayed
    set_axes_equal(ax)

def get_atom_bonds_coordinates_from_distance_df(distance_df, xyzarr, index_atoms=0):
    """
    WRITE ME
    """
    #get atom numbers from sel_dist2 data frame, e.g. C11 --> 11
    atom1_num=distance_df['atom1_idx'].apply(atom_idx_to_index).astype(int).tolist()
    atom2_num=distance_df['atom2_idx'].apply(atom_idx_to_index).astype(int).tolist()
    #distance as bond label
    bond_length=distance_df['distance_calc'].apply(lambda x: '{:.3f}'.format(x)).tolist()
    #if first atom has index 1
    if index_atoms:
            atom1_num=[idx - 1 for idx in atom1_num]
            atom2_num=[idx - 1 for idx in atom2_num]
    #atom1 and atom 2 coordinates
    atom1_coord=xyzarr[atom1_num]
    atom2_coord=xyzarr[atom2_num]
    #put atom1 and atom2 coordinats in an numpy array (for bonds display)
    bond_coord=np.array(list(zip(atom1_coord,atom2_coord)))
    return atom1_num, atom1_coord, atom2_coord, zip(bond_coord, bond_length), zip(atom1_num, atom2_num, bond_length)
        
def set_axes_equal(ax):
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.
    https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
    Void function
    ----------
    Parameters
    ----------
    ax : ??.
    Matplotlib axis

    -------
    Returns
    -------
    ???
    ???
    --------
    Examples
    --------
    WRITE ME        
    """
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    
    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    # set to 0.38 --> bigger molecules
    plot_radius = 0.35*max([x_range, y_range, z_range])
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

class Arrow3D(FancyArrowPatch):
    """
    draw fancy arrows in x y z
    https://gist.github.com/WetHat/1d6cd0f7309535311a539b42cccca89c
    ----------
    Parameters
    ----------
    vec1 : ??.
    Array of cartesian coordinates of atoms in the molecule ??

    vec2 : ??.
    Array of cartesian coordinates of atoms in the molecule ??

    coord : ??.
    Array of cartesian coordinates of atoms in the molecule ??
    -------
    Returns
    -------
    ???
    ???
    --------
    Examples
    --------
    WRITE ME        
    """
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
            super().__init__((0, 0), (0, 0), *args, **kwargs)
            self._xyz = (x, y, z)
            self._dxdydz = (dx, dy, dz)
    def draw(self, renderer):
            x1, y1, z1 = self._xyz
            dx, dy, dz = self._dxdydz
            x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)
            xs, ys, zs = proj3d.proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
            super().draw(renderer)
    def do_3d_projection(self, renderer=None):
            x1, y1, z1 = self._xyz
            dx, dy, dz = self._dxdydz
            x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)
            xs, ys, zs = proj3d.proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
            return np.min(zs) 
	
def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    """
    Add an 3d arrow to an `Axes3D` instance.
    """
    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)
	
setattr(Axes3D, 'arrow3D', _arrow3D) #end draw fancy arrows in x y z

## xyz_df functions
def get_xyz_headers(xyz_filename):
    """
    Write Me
    """
    with open(xyz_filename, 'r') as xyz_file:
            raw_headers=xyz_file.readlines()[:2]
    if raw_headers[-1]=="\n":
            raw_headers[-1]=" \n"
    return ''.join(raw_headers).rstrip('\n') 
       

def get_xyz_df_from_file(xyz_filename):
    """
    Write Me
    """
    df=load_single_xyz_file(xyz_filename)
    return df
    

def expend_xyz_df(xyz_df, index_atoms=0):
    """
    xyz_df.columns=["element", "x", "y", "z"]
    expended_xyz_df:
       atom1_idx atom2_idx element       x      y       z  weight  cov_radius
    0         C0        C0       C  -1.944  0.247  -0.071  12.011        0.76
    1         C1        C1       C  -3.293 -0.055   0.126  12.011        0.76
    2         C2        C2       C  -4.288  0.537  -0.626  12.011        0.76        
    """
    expended_xyz_df=xyz_df.copy()
    expended_xyz_df['atom1_idx']=["{}{}".format(atm, idx) for atm, idx in zip(expended_xyz_df.element, expended_xyz_df.index.array)] #element + position in xyz C --> C0
    if index_atoms: #first atom has now index 1 in atom name, e.g. C --> C1
            expended_xyz_df['atom1_idx']=["{}{}".format(atm, idx+1) for atm, idx in zip(expended_xyz_df.element, expended_xyz_df.index.array)]
    expended_xyz_df['atom2_idx']=expended_xyz_df['atom1_idx'] #atom 1 is atom 2, duplicate for later use
    #atomic weight dict
    expended_xyz_df['weight']=expended_xyz_df['element'].apply(lambda x: TABConstants.atomic_weights.value[x]) # gemmi - alt: xyz_df['weight'] = xyz_df['element'].apply(lambda x: gemmi.Element(x).weight)
    #covalent radius dict
    expended_xyz_df['cov_radius']=expended_xyz_df['element'].apply(lambda x: TABConstants.covalent_radii.value[x]) # gemmi - alt: xyz_df['cov_radius'] = xyz_df['element'].apply(lambda x: gemmi.Element(x).covalent_r)
    #reorder data frame
    expended_xyz_df=expended_xyz_df[['atom1_idx','atom2_idx','element','x','y','z','weight','cov_radius']]
    return expended_xyz_df
        
def get_xyz_stats(xyz_df, radii_ext=0.08):
    """
    info_df:

    sum_formula_list:

    """
    grouped=xyz_df.groupby('element',sort=False)['atom1_idx'].count()
    fw=xyz_df['weight'].sum() #formula weight from sum of atomic weights
    info_df=grouped.reset_index() #results of grouping into new data frame
    info_df['utf_num']=info_df['atom1_idx'].apply(lambda x: num_to_subnum(x)) #numbers to utf8 subscript number, e.g. 2 --> ₂ for sum formula
    info_df.loc[info_df['utf_num'] == '₁', 'utf_num'] = '' #replace '1' with '', e.g. C1O2 --> CO2
    info_df=info_df.sort_values(by=['element']) #sort elements alphabetically
    sum_formula_list=list("{}{}".format(element, number) for element, number in zip(info_df.element, info_df.utf_num)) #generate the sum formula out of elements and the number of each element, e.g. C 5 H 12 --> C5H12
    info_df=info_df.drop(columns=['utf_num']) #drop the utf8 subscript number column
    #calculate mass fraction of each element
    info_df['atom%']=info_df.apply(lambda x: TABConstants.atomic_weights.value[x.element]*x.atom1_idx/fw*100, axis=1) # gemmi - alt: info_df['atom%']=info_df.apply(lambda x: gemmi.Element(x.element).weight*x.atom1_idx/fw*100, axis=1)
    #atomic radii from the gemmi table
    info_df['radii']=info_df.apply(lambda x: TABConstants.covalent_radii.value[x.element], axis=1) # gemmi - alt: info_df['radii']=info_df.apply(lambda x: gemmi.Element(x.element).covalent_r, axis=1)
    #atomic radii + x% used for calculation of bonds
    info_df['radii_plus']=info_df.apply(lambda x: x.radii + x.radii*radii_ext, axis=1)
    return info_df, sum_formula_list

## dist_df_functions
def dist_matrix_to_series1(dist_matrix): #dm_to_series1
    """
    removes the upper triangle of the distance matrix and zeros
    e.g., from d(A-B) = 1.234 Å = d(B-A) =1.234 Å, d(B-A) will be removed 
    d(A-B) = 0 Å will be removed as well
       C0  C1  C2           C0  C1  C2
    C0 0.0 1.1 2.3   ==\ C0 NaN NaN NaN
    C1 1.1 0.0 1.5   ==/ C1 1.1 NaN NaN
    C2 2.3 1.5 0.0       C2 2.3 1.5 NaN
    """
    dist_matrix=dist_matrix.astype(float) # do not comment this, angle list will be incomplete
    dist_matrix.values[np.triu_indices_from(dist_matrix, k=1)] = np.nan
    dist_matrix=dist_matrix.replace(0, np.nan) #replace zeros with nan
    return dist_matrix.unstack().dropna() #return and drop all nan

def get_distance_df(xyz_df, reduce_distance_matrix=True):
    """
    dist_mat_full:
       C0  C1  C2
    C0 0.0 1.1 2.3
    C1 1.1 0.0 1.5
    C2 2.3 1.5 0.0

    dist_df:
    
    """
    dist_mat_full=pd.DataFrame(squareform(pdist(xyz_df.iloc[:,3:6],'euclid')), #iloc [:,3:6] contains xyz coordinates
                               columns=xyz_df[['atom1_idx','element','cov_radius']],
                               index=xyz_df[['atom2_idx','element','cov_radius']])
    if reduce_distance_matrix:
            dist_mat_red=dist_matrix_to_series1(dist_mat_full) #remove the upper triangle and zeros
            dist_mat_red=dist_mat_red.reset_index(level=[1]) #bring it to a "normal" form
            dist_df=pd.DataFrame(dist_mat_red.to_records()) #bring it to a "normal" form, distance matrix --> disctance data frame
            normal_form_column='index'
            dist_rename_dict={'0':'distance_calc'}
    else:
            dist_mat = dist_mat_full.unstack()
            dist_mat = dist_mat.reset_index()
            dist_df=dist_mat
            normal_form_column='level_0'
            dist_rename_dict={0:'distance_calc'}
    #bring it to a "normal" form ...
    dist_df[['atom1_idx','element1','cov_radius1']]=pd.DataFrame(dist_df[normal_form_column].tolist(), index=dist_df.index)
    dist_df[['atom2_idx','element2','cov_radius2']]=pd.DataFrame(dist_df['level_1'].tolist(), index=dist_df.index)
    dist_df.drop([normal_form_column, 'level_1'], axis=1,inplace=True)
    dist_df.rename(columns=dist_rename_dict,inplace=True)
    dist_df=dist_df[['atom1_idx','element1','cov_radius1','atom2_idx','element2','cov_radius2','distance_calc']] #reorder data frame
    return dist_df, dist_mat_full

def determine_bonds(dist_df, determine_by='radii', radii_ext=0.08, ref_bonds_mask=None, include_bonds=None):
    """
    Write Me
    """
    new_dist_df=dist_df.copy()
    if determine_by=='radii':
            new_dist_df['distance_radii']=(new_dist_df['cov_radius1']+new_dist_df['cov_radius2'])*(1+radii_ext) #column with the sum of the atomic radii from two elements /atoms + x%
            new_dist_df['is_bond']=(new_dist_df['distance_calc']<new_dist_df['distance_radii']) #distance is considered as bond if the calculated distance is smaller than the sum of the atomic radii
    if determine_by=='ref':
            new_dist_df['is_bond']=ref_bonds_mask                
    #include distances from selected atom (pairs) from args include connections sets 'is_bond' to 'True' if 'False'
    if include_bonds:
            #must have atoms 1 and 2 and 'is_bond' must be 'False'
            new_dist_df.loc[(new_dist_df.atom1_idx.isin(include_bonds) & 
                                    new_dist_df.atom2_idx.isin(include_bonds) & 
                               ~new_dist_df.is_bond), 'is_bond'] = True
            #np variant of the above, slightly slower
    ##	dist_df['is_bond'] = np.where((dist_df.atom1_idx.isin(args.includeCon)) &
    ##	(dist_df.atom2_idx.isin(args.includeCon) & (~dist_df.is_bond)), True, dist_df['is_bond'])
    return new_dist_df

def expend_distance_df(dist_df, radii_ext=0.08, include_bonds=False, bonds_by_ref=False, bonds_boolean_mask=None):
    """
    Write Me
    """
    if bonds_by_ref:
            expended_dist_df=determine_bonds(dist_df, determine_by='ref', ref_bonds_mask=bonds_boolean_mask, include_bonds=include_bonds)
    else:
            expended_dist_df=determine_bonds(dist_df, determine_by='radii', radii_ext=radii_ext, include_bonds=include_bonds)
    return expended_dist_df


def decorate_distance_df(dist_df):
    """
    Write Me
    """
    expended_dist_df=dist_df.copy()
    #fusion char, A B --> A-B
    expended_dist_df['fusion_char']='–'
    #A0 B1 - --> A0-B1
    expended_dist_df['A-B']=expended_dist_df['atom1_idx']+expended_dist_df['fusion_char']+expended_dist_df['atom2_idx']
    #A B - --> A-B
    expended_dist_df['El1-El2']=expended_dist_df['element1']+expended_dist_df['fusion_char']+expended_dist_df['element2']
    #A B --> AB
    expended_dist_df['ElEl']=expended_dist_df['El1-El2'].apply(lambda x: ''.join(sorted(x)))
    return expended_dist_df

                                                 
def exclude_distances(dist_df, exclude_atoms=False, exclude_elements=False):
    """
    Write Me
    """
    new_dist_df=dist_df.copy()
    if exclude_atoms: #exclude named atoms from input (in Atom1 and Atom2)
        new_dist_df=new_dist_df[~new_dist_df.atom1_idx.isin(exclude_atoms) & ~new_dist_df.atom2_idx.isin(exclude_atoms)]
    if exclude_elements: #exclude named elements from input (in Element1 and Element2)
        new_dist_df=new_dist_df[~new_dist_df.element1.isin(exclude_elements) & ~new_dist_df.element2.isin(exclude_elements)]
    if len(new_dist_df.index) == 0: #exit if selected bonds data frame is empty 
        print("No bonds found. Include more atoms or elements. Exit.")
        sys.exit(1)
    return new_dist_df

def sort_distances(dist_df, sort_bonds_asc=False, sort_bonds_des=False, sort_elements_asc=False, sort_elements_des=False):
    """
    Write Me
    """
    new_dist_df=dist_df.copy()
    if sort_bonds_asc: #sort bond length values ascending
        new_dist_df=new_dist_df.sort_values(by=['distance_calc'])
    if sort_bonds_des: #sort bond length values descending
        new_dist_df=new_dist_df.sort_values(by=['distance_calc'],ascending=False)
    if sort_elements_asc: #sort by elements ascending, A --> Z (not PSE like)
        new_dist_df=new_dist_df.sort_values(by=['element1','element2','distance_calc','atom1_idx','atom2_idx'])
    if sort_elements_des: #sort by elements descending, A --> Z (not PSE like)
        new_dist_df=new_dist_df.sort_values(by=['element1','element2','distance_calc','atom1_idx','atom2_idx'],ascending=False)
    return new_dist_df

#angles_df
def get_angles_df(dist_df, xyz_df, xyzarr):
    #group atom2 by atom1,e.g. C0 C1, C0 C2, C0 C3...
    group1 = dist_df.groupby('atom1_idx',sort=False)['atom2_idx']
    #make 4 empty lists
    atom1 = list()
    atom2 = list()
    atom3 = list()
    anglelist = list()
    #middle atom 'B' (for angle A-B-C) is in name of the group
    #get x,y,z coordinates (a2) of 'B' from xyz data frame
    for name, group in group1:
        a2=xyz_df.index[xyz_df['atom1_idx'] == name].tolist()
        #'A' and 'C' atoms are in the group
        #get x,y,z coordinates of 'A' (a1) and 'C' (a3) from xyz data frame 
        for s in combinations(group,2):
            #very nice itertool, e.g.:
            #a1 (central atom) binds to a2, a3, a4, a5
            #angles will be a2-a1-a3, a2-a1-a4, a2-a1-a5, a3-a1-a4, a3-a1-a5, a4-a1-a5
            #exludes double entries like a3-a1-a2, a4-a1-a2,....
            a1=xyz_df.index[xyz_df['atom1_idx'] == s[0]].tolist()
            a3=xyz_df.index[xyz_df['atom1_idx'] == s[1]].tolist()
            #calculate the angle
            angle = calc_angle(xyzarr, *a1, *a2, *a3)
            #name of atom1 ('A') --> atom1 list
            atom1.append(s[0])
            #name of atom2 ('B') --> atom2 list
            atom2.append(name)
            #name of atom3 ('C') --> atom3 list
            atom3.append(s[1])
            #calculated angle to list of angles
            anglelist.append(angle)
    #all 4 lists in angles_df data frame
    angles_df=pd.DataFrame(({'atom1_idx': atom1, 'atom2_idx': atom2, 'atom3_idx': atom3, 'angle_calc': anglelist}))
    #construct elements from atom names, e.g. C1 --> C, Fe13 --> Fe
    angles_df['element1']=angles_df['atom1_idx'].apply(atom_idx_to_element)
    angles_df['element2']=angles_df['atom2_idx'].apply(atom_idx_to_element)
    angles_df['element3']=angles_df['atom3_idx'].apply(atom_idx_to_element)
    return angles_df

def decorate_angles_df(angles_df):
    """
    Write Me
    """
    expended_angles_df=angles_df.copy()
    expended_angles_df['fusion_char'] = '–' #fuse atom names A B C by '-' --> A-B-C
    #A0 B1 C2- --> A0-B1-C2
    expended_angles_df['A-B-C'] = expended_angles_df['atom1_idx']+expended_angles_df['fusion_char']+expended_angles_df['atom2_idx']+expended_angles_df['fusion_char']+expended_angles_df['atom3_idx']
    #A B C--> A-B-C
    expended_angles_df['El1-El2-El3'] = expended_angles_df['element1']+expended_angles_df['fusion_char']+expended_angles_df['element2']+expended_angles_df['fusion_char']+expended_angles_df['element3']
    expended_angles_df['El1-El3'] = expended_angles_df['element1']+expended_angles_df['fusion_char']+angles_df['element3'] #A (B) C--> A-C, for grouping
    expended_angles_df['ElElEl'] = expended_angles_df['El1-El2-El3'].apply(lambda x: ''.join(sorted(x))) #A-B-C--> ABC, for grouping
    expended_angles_df['ElEl'] = expended_angles_df['El1-El3'].apply(lambda x: ''.join(sorted(x))) #A-C--> AC, for grouping
    expended_angles_df['ElEl_ElElEl'] = expended_angles_df['ElEl'] + expended_angles_df['ElElEl'] #AC + ABC --> ACABC, for grouping
    return expended_angles_df
                                                 
def exclude_angles(angles_df, exclude_atoms=False, exclude_elements=False):
    new_angles_df=angles_df.copy()
    if exclude_atoms: #exclude named atoms from input (in Atom1 and Atom2 and Atom3)
        new_angles_df=new_angles_df[~new_angles_df.atom1_idx.isin(exclude_atoms) & ~new_angles_df.atom2_idx.isin(exclude_atoms) & ~new_angles_df.atom3_idx.isin(exclude_atoms)] 
    if exclude_elements: #exclude named elements from input (in Element1 and Element2 and Element3)
        new_angles_df=new_angles_df[~new_angles_df.element1.isin(exclude_elements) & ~new_angles_df.element2.isin(exclude_elements) & ~new_angles_df.element3.isin(exclude_elements)]     
    if len(new_angles_df) == 0: #exit if selected angles data frame is empty
        print("No angles found. Include more atoms or elements. Exit.")
        sys.exit(1)
    return new_angles_df

def sort_angles(angles_df, sort_angles_asc=False, sort_angles_des=False, sort_elements_asc=False, sort_elements_des=False):
    new_angles_df=angles_df.copy()
    if sort_angles_asc: #sort angle values ascending
        new_angles_df=new_angles_df.sort_values(by=['angle_calc'])
    if sort_angles_des: #sort angle values descending
        new_angles_df=new_angles_df.sort_values(by=['angle_calc'],ascending=False)
    if sort_elements_asc: #sort by elements ascending, A --> Z (not PSE like)
        new_angles_df=new_angles_df.sort_values(by=['element1','element2','element3','angle_calc','atom1_idx','atom2_idx','atom3_idx'])
    if sort_elements_des: #sort by elements descending, A --> Z (not PSE like)
        new_angles_df=new_angles_df.sort_values(by=['element1','element2','element3','angle_calc','atom1_idx','atom2_idx','atom3_idx'],ascending=False)
    return new_angles_df

#Connectivity
def get_expended_dist_df(xyz_df, radii_ext=0.08, include_bonds=False, bonds_by_ref=False, bonds_boolean_mask=None):
    expended_xyz_df=expend_xyz_df(xyz_df)
    dist_df, _=get_distance_df(expended_xyz_df)
    expended_dist_df=expend_distance_df(dist_df, radii_ext, include_bonds, bonds_by_ref, bonds_boolean_mask)
    return expended_dist_df

def connect_the_dots(xyz_df, radii_ext=0.08, include_bonds=False, bonds_by_ref=False, bonds_boolean_mask=None):
    expended_dist_df=get_expended_dist_df(xyz_df, radii_ext, include_bonds, bonds_by_ref, bonds_boolean_mask)
    sel_dist=pd.DataFrame(expended_dist_df[expended_dist_df.is_bond])
    return sel_dist

def get_bond_pairs(expended_dist_df):
    indices_1=expended_dist_df['atom1_idx'].apply(atom_idx_to_index).tolist() # remove ".apply(atom_idx_to_idx)"?
    indices_2=expended_dist_df['atom2_idx'].apply(atom_idx_to_index).tolist()
    bond_pairs=zip(indices_1, indices_2)
    return bond_pairs

def get_molecule_connections(bond_pairs):
    import networkx as nx
    nx_graph=nx.from_edgelist(bond_pairs)
    molecule_atoms=list(nx.connected_components(nx_graph))
##    map_dict={z:x for x, y in (enumerate(molecule_atoms)) for z in y} #{atom_idx : number_of_the_molecule}
    return molecule_atoms

def get_aligned_molecule_connections(bond_pairs_1, bond_pairs_2):
    from networkx.algorithms import isomorphism
    import networkx as nx
    graph_1=nx.from_edgelist(bond_pairs_1)
    graph_2=nx.from_edgelist(bond_pairs_2)
    graph_matcher=isomorphism.GraphMatcher(graph_2, graph_1)
    if not graph_matcher.is_isomorphic():
        return False
    else:
        relabeled_graph_2=nx.relabel_nodes(graph_2, graph_matcher.mapping, copy=False)
        return relabeled_graph_2  

#Analyzers
class OverlayAnalyzer():
    def __init__(self, parser,xyz_filenames ,xyz_dfs, fit_mode='all', atoms=['1', '2']):
        """
        fkf
        """
        self.xyz_dfs=xyz_dfs
        set_plot_settings()
        if fit_mode=='all':
            self.args=parser.parse_args(xyz_filenames+['-aa', '-vcm', '-s'])
        if fit_mode=='atoms':
            self.args=parser.parse_args(xyz_filenames+['-a']+atoms+['-vcm', '-s'])
            self.atom_indices=self.args.atoms[0]
        self.overlay_xyzs(fit_mode)

                    
    def get_inital_data(self):
        try:
            # self.head_list=[get_xyz_headers(xyz_filename) for xyz_filename in self.args.filename]
            self.xyz_df_list=self.xyz_dfs
            
            #index starts at 1, first atom is atom 1
            for xyz_df in self.xyz_df_list:
                xyz_df.index+=1
        except IOError: #file not found
            print(f"File(s) not found. Exit.")
            sys.exit(1)

    def overlay_xyzs(self, fit_mode):
        self.get_inital_data()
        if fit_mode=='all':
            self.overlay_xyz_by_all()
        if fit_mode=='atoms': 
            self.overlay_xyz_by_atoms()
        self.plot_by_molecules()
        # self.save_altered_xyz()
            

    def overlay_xyz_by_atoms(self):
        try:
            atoms_mol1 = self.xyz_df_list[0][['x','y','z']].iloc[self.atom_indices]
            centroid0 = np.mean(self.xyz_df_list[0][['x','y','z']].iloc[self.atom_indices], axis=0)
            for xyz_df in self.xyz_df_list:
                centroid = np.mean(xyz_df[['x','y','z']].iloc[self.atom_indices], axis=0)
                xyz_df[['x','y','z']]=xyz_df[['x','y','z']].apply(lambda x: x-centroid, axis=1)
                xyz_df[['x','y','z']]=align_xyz(xyz_df[['x','y','z']].iloc[self.atom_indices],
                                                atoms_mol1, xyz_df[['x','y','z']])
        except:
            #exit if xyz files have differnet number of atoms
            print('Warning! Number of atoms must be equal in all xyz files. Exit.')
            sys.exit(1)                        

    def overlay_xyz_by_all(self):
        try:
            self.xyz_df_list[0][['x', 'y', 'z']]=self.xyz_df_list[0][['x', 'y', 'z']].apply(rmsd_center_coords) #center the molecule by substracting the centroid from the coordinates
            atoms_mol1=self.xyz_df_list[0][['x', 'y', 'z']]
            for xyz_df in self.xyz_df_list:
                xyz_df[['x', 'y', 'z']], *_=rmsd_align_xyz(xyz_df[['x', 'y', 'z']], atoms_mol1) 
        except ValueError:
            #exit if xyz files have differnet number of atoms 
            print('Warning! Number of atoms must be equal in all xyz files. Exit.')
            sys.exit(1)

    def plot_by_molecules(self):
        fig, ax=init_figure_and_axis()
        num_of_xyz = len(self.xyz_df_list)
        a1_a2_bond_list=[]
        labels_list=[]
        for xyz_df in self.xyz_df_list:
            dist_df, *_=get_distance_df(expend_xyz_df(xyz_df), reduce_distance_matrix=True) #self.args.radius
            bonded_dist_df=expend_distance_df(dist_df)
            dist_df['is_bond']=bonded_dist_df.is_bond
            sel_dist=pd.DataFrame(dist_df[(dist_df.is_bond == True)])
            sel_dist=exclude_distances(dist_df=sel_dist,
                                       exclude_atoms=self.args.excludeAt,
                                       exclude_elements=self.args.excludeEl,
                                       )
            xyzarr=xyz_df.iloc[:,3:6].to_numpy()
##                        labels=[np.concatenate(row, bonded_dist_df['atom1_idx'].iloc[index])  for index, row in enumerate(xyzarr)]
##                        labels_list.append(labels)
            labels=list(zip(xyz_df['x'], xyz_df['y'], xyz_df['z'], xyz_df.index))
##                        print(labels)
            labels_list.append(labels)
            *_, bond_mat_test = get_atom_bonds_coordinates_from_distance_df(sel_dist, xyzarr, index_atoms=True)
            a1_a2_bond_list.append(pd.DataFrame(bond_mat_test, columns=['atom_1', 'atom_2','bond_length']))              
        
        #get the total number of atoms from all xyz files 
        #minus excluded atoms for the atom size in the plot
        num_atom_xyz = sum(len(xyz_df) for xyz_df in self.xyz_df_list)                

#                
        #use Matplotlib colormap
        if self.args.cmap:                                      #WONT WORK DUE TO ENUM 
            try:
                new_colors = [plt.get_cmap(self.args.cmap)(1. * i/num_of_xyz) for i in range(num_of_xyz)]
                color_cycle = cycler('color',new_colors)
            except ValueError:
                print('Warning! Not a valid Matplotlib colormap. Applying default colormap instead.')
        
        #get the coordinates of the atoms and put them in an array
        #also apply the color cycle / colormap
        for count, (xyz_df, style) in enumerate(zip(self.xyz_df_list,cycle(cycler(c=['b', 'r', 'g', 'c', 'm', 'y'])))):  #PlotConstants.COLOR_CYCLE.value
            atoms1=a1_a2_bond_list[count]['atom_1']
            atoms2=a1_a2_bond_list[count]['atom_2']
            atom1_coord = xyz_df[['x','y','z']].loc[atoms1+1]
            atom2_coord = xyz_df[['x','y','z']].loc[atoms2+1]
            atom1_2_coord = np.array(list(zip(atom1_coord.to_numpy(),atom2_coord.to_numpy())))
            
            #scatter (atom) plot
            ax.scatter(*xyz_df[['x','y','z']].to_numpy().T,s=np.log10(PlotConstants.ATOM_SCALER.value/num_atom_xyz),alpha=PlotConstants.ALPHA_ATOMS.value,**style)
            #bonds
            for bonds in atom1_2_coord:
                ax.plot(*bonds.T,linewidth=np.log10(PlotConstants.BOND_SCALER.value/num_atom_xyz),alpha=PlotConstants.ALPHA_BONDS.value,**style)
            #labels
            for x, y, z, label in labels_list[count]:
                ax.text(x+0.12 , y+0.12, z+0.12, label, fontsize=100/len(xyz_df.index) + 8 , color='black')

        finallize_figure_and_axis(fig, ax)
        ##save plot
        plot_filename=self.args.filename[0]+'_'+self.args.filename[1]+'.png'
        plt.savefig(plot_filename)
        #show the plot
        plt.show()
        plt.clf()
        
        
    def save_altered_xyz(self):
        try:
            for file_index, filename in enumerate(self.args.filename):
                my_numpy = self.xyz_df_list[file_index].to_numpy()
                np.savetxt(os.path.splitext(filename)[0] +'-mod.xyz',
                           my_numpy, fmt='%-2s  %12.8f  %12.8f  %12.8f', delimiter='',
                           header=self.head_list[file_index], comments='')
        except IOError: #write error -> exit here
            print("Write error. Exit.")
            sys.exit(1)

class TabDataAnalyzer():
    def __init__(self, parser,origin_df, xyz_filename, get_plot=False, bonds_by_ref=False, bonds_ref=None):
        if get_plot:
            self.args=parser.parse_args([xyz_filename, '-sb'])
            set_plot_settings()
        else:
            self.args=parser.parse_args([xyz_filename])
        self.bonds_by_ref=bonds_by_ref
        self.bonds_ref=bonds_ref
        self.plot_filename=None
        self.origin_df=origin_df
        
    def get_inital_data(self):
        """
        read xyz into data frame (skip first two rows of the xyz file)
        only XMol/xyz is supported, atom(as element) x y z, e.g. C 1.58890 -1.44870 -0.47000
        """
        # try:
        #     xyz_df=get_xyz_df_from_file(self.args.filename)
        # except IOError: #file not found
        #     print(f"'{self.args.filename}'" + " not found")
        #     sys.exit(1)
        #xx.x% --> 0.xxx
        self.radii_ext=self.args.radius/100
        #expend xyz_df
        self.xyz_df=expend_xyz_df(self.origin_df, index_atoms=self.args.index)
        #get total number of each element by counting elements and grouping
        self.info_df, self.sum_formula_list=get_xyz_stats(self.xyz_df, self.radii_ext)
        

    def get_bonds_data(self):
        dist_df, self.dist_mat_full=get_distance_df(self.xyz_df, reduce_distance_matrix=True)
        bonded_dist_df=expend_distance_df(dist_df, self.radii_ext, include_bonds=self.args.includeCon, bonds_by_ref=self.bonds_by_ref ,bonds_boolean_mask=self.bonds_ref)
        self.dist_df=decorate_distance_df(bonded_dist_df)
        #all_dist data frame --> sel_dist data frame if 'is_bond' is True
        sel_dist=pd.DataFrame(self.dist_df[(self.dist_df.is_bond == True)])
        sel_dist=exclude_distances(dist_df=sel_dist,
                                   exclude_atoms=self.args.excludeAt,
                                   exclude_elements=self.args.excludeEl,
                                   )
        self.sel_dist=sort_distances(dist_df=sel_dist,
                                     sort_bonds_asc=self.args.sortasc,
                                     sort_bonds_des=self.args.sortdes,
                                     sort_elements_asc=self.args.sortascEl,
                                     sort_elements_des=self.args.sortdesEl,
                                     )

    def get_summary_bonds(self):
        #lists for printed tables
        summary_bond_table_1 = list()
        summary_bond_table_2 = list()
        summary_bond_table_3 = list()

        #group El-El and distances by ElEl, e.g. C-N 1.234 Å by CN
        grouped = self.sel_dist[['El1-El2','distance_calc']].groupby(self.sel_dist['ElEl'],sort=False)

        #verbose table El1-El2 | bond length, e.g. C-C 1.223, 1.456, 1.511
        for groups in grouped:
            summary_bond_table_1.append([groups[1].iloc[0].tolist()[0], 
                    ', '.join(groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist())])

        #short table El1-El2 | bond length, e.g. C-C 1.223 - 1.511 (for >2), or C-C 1.223 / 1.511 (for 2), C-C 1.223 (for one)
        # float to 4 decimals, e.g. 0.1234 
        for groups in grouped:
            if len(groups[1]) == 1:
                summary_bond_table_2.append([groups[1].iloc[0].tolist()[0], 
                        groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist()[0]])
            elif len(groups[1]) == 2:
                summary_bond_table_2.append([groups[1].iloc[0].tolist()[0], 
                        groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist()[0] + 
                        ' / ' + groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist()[-1]])
            else:
                summary_bond_table_2.append([groups[1].iloc[0].tolist()[0], 
                        groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist()[0] + 
                        ' - ' + groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist()[-1]])

        #grouped = sel_dist[['El1-El2','distance_calc']].groupby(sel_dist['ElEl'],sort=False)

        #generate table with statistics, | El1-El2 | Count | Mean | Median | Sam. std. dev. | Pop. std. dev. | Std. error |
        for groups in grouped:
            summary_bond_table_3.append([groups[1].iloc[0].tolist()[0], groups[1].distance_calc.count(),f'{groups[1].distance_calc.mean():.4f}', \
                    f'{groups[1].distance_calc.median():.4f}', f'{groups[1].distance_calc.std():.4f}', f'{groups[1].distance_calc.std(ddof=0):.4f}', \
                    f'{groups[1].distance_calc.sem():.4f}',f'{groups[1].distance_calc.skew():.4f}'])
        self.summary_bond_table_1=summary_bond_table_1
        self.summary_bond_table_2=summary_bond_table_2
        self.summary_bond_table_3=summary_bond_table_3
##        return summary_bond_table_1, summary_bond_table_2, summary_bond_table_3

    def get_angles_data(self):
        dist_df, *_=get_distance_df(self.xyz_df, reduce_distance_matrix=False)
        bonded_dist_df=expend_distance_df(dist_df, self.radii_ext, include_bonds=self.args.includeCon, bonds_by_ref=self.bonds_by_ref ,bonds_boolean_mask=self.bonds_ref)
        dist_df['is_bond']=bonded_dist_df.is_bond&(dist_df.distance_calc>0)
        sel_dist=pd.DataFrame(dist_df[(dist_df.is_bond == True)])
        sel_dist=exclude_distances(dist_df=sel_dist,
                                   exclude_atoms=self.args.excludeAt,
                                   exclude_elements=self.args.excludeEl,
                                   )
        sel_dist=sort_distances(dist_df=sel_dist,
                                sort_bonds_asc=self.args.sortasc,
                                sort_bonds_des=self.args.sortdes,
                                sort_elements_asc=self.args.sortascEl,
                                sort_elements_des=self.args.sortdesEl,
                                )
        self.sel_dist2=sel_dist
        self.xyzarr = self.xyz_df.iloc[:,3:6].to_numpy() #xyz array with coordinates of all atoms is needed for angle (and dihedral) calculation
        sel_angles=get_angles_df(sel_dist, self.xyz_df, self.xyzarr)
        sel_angles=decorate_angles_df(sel_angles)
        sel_angles=exclude_angles(angles_df=sel_angles,
                                  exclude_atoms=self.args.excludeAt,
                                  exclude_elements=self.args.excludeEl,
                                  )
        sel_angles=sort_angles(angles_df=sel_angles,
                               sort_angles_asc=self.args.sortasc,
                               sort_angles_des=self.args.sortdes,
                               sort_elements_asc=self.args.sortascEl,
                               sort_elements_des=self.args.sortdesEl,
                               )
        self.sel_angles=sel_angles

    def get_summary_angles(self):
        #lists for printed tables
        summary_angle_table_1 = list()
        summary_angle_table_2 = list()
        summary_angle_table_3 = list()

        #group El1-El2-El3 and angles by ElElEl, e.g. O-C-N 1.234 Å by OCOCN sorted CCNOO
        grouped = self.sel_angles[['El1-El2-El3','angle_calc']].groupby(self.sel_angles['ElEl_ElElEl'],sort=False)

        #verbose table El1-El2-El3 | angle, e.g. C-C-C 122.31, 145.61, 151.11
        for groups in grouped:
            summary_angle_table_1.append([groups[1].iloc[0].tolist()[0], 
                    ', '.join(groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist())])

        #short table El1-El2-El3 | angle, e.g.C-C-C 122.32 - 151.11 (for >2), 
        #or C-C-C 122.32 / 151.11 (for 2), C-C-C 122.32 (for one)
        for groups in grouped:
            if len(groups[1]) == 1:
                summary_angle_table_2.append([groups[1].iloc[0].tolist()[0], 
                        groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist()[0]])
            elif len(groups[1]) == 2:
                summary_angle_table_2.append([groups[1].iloc[0].tolist()[0], 
                        groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist()[0] + 
                        ' / ' + groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist()[-1]])
            else:
                summary_angle_table_2.append([groups[1].iloc[0].tolist()[0], 
                        groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist()[0] + 
                        ' - ' + groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist()[-1]])
                        
        #grouped = sel_angles[['El1-El2-El3','angle_calc']].groupby(sel_angles['ElEl_ElElEl'],sort=False)

        #generate table with statistics, | El1-El2-El3 | Count | Mean | Median | Sam. std. dev. | Pop. std. dev. | Std. error |
        for groups in grouped:
            summary_angle_table_3.append([groups[1].iloc[0].tolist()[0], groups[1].angle_calc.count(),f'{groups[1].angle_calc.mean():.2f}', \
                    f'{groups[1].angle_calc.median():.2f}', f'{groups[1].angle_calc.std():.2f}', f'{groups[1].angle_calc.std(ddof=0):.2f}', \
                    f'{groups[1].angle_calc.sem():.2f}',f'{groups[1].angle_calc.skew():.2f}'])
        self.summary_angle_table_1=summary_angle_table_1
        self.summary_angle_table_2=summary_angle_table_2
        self.summary_angle_table_3=summary_angle_table_3
##        return summary_angle_table_1, summary_angle_table_2, summary_angle_table_3

    def plot_the_molecule(self):
        atom_indices, atom1_coord, atom2_coord, bond_coord, *_ = get_atom_bonds_coordinates_from_distance_df(self.sel_dist2, self.xyzarr, index_atoms=self.args.index)
        fig, ax=init_figure_and_axis()

        #for assigning colors to different elemments
        leftover_atoms=atom_indices.copy()
        scaling_factor=50/len(self.xyzarr)
        for key, value in PlotMapping.color_grouping.value.items():
            if key!='leftover_atoms':
                atom_array=self.xyz_df.index[self.xyz_df['element'].isin(value.get('elements'))].tolist()
                ax.scatter(*self.xyzarr[atom_array].T,s=scaling_factor*value.get('size'),color=value.get('color'), alpha=value.get('alpha'))
                leftover_atoms=[item for item in leftover_atoms if item not in atom_array]

        leftover_dict=PlotMapping.color_grouping.value.get('leftover_atoms')
        ax.scatter(*self.xyzarr[leftover_atoms].T,s=scaling_factor*leftover_dict.get('size'),color=leftover_dict.get('color'))
        
        #label atoms
        #show no lables if -sn option is activated
        if not self.args.shownl:
            atom_label=self.xyz_df['atom1_idx'].tolist()
            atom_coord_name = zip(self.xyzarr,atom_label)
            for coord, label in atom_coord_name:
                ax.text(*(coord+0.12).T, label, fontsize=100/len(self.xyzarr) + 8 , color='black')
        
        #draw bonds
        for bonds, labels in bond_coord:
            ax.plot(*bonds.T, color='gray', linewidth=3.0)
            if self.args.showbl:
                #show bond labels
                ax.text(*np.average(bonds+0.06,axis=0).T,labels,fontsize=(100/len(self.xyzarr) + 8)/1.5,color='gray')
        
        #define ranges forplanes
        x_pl=np.sort(self.xyzarr[:,0])
        y_pl=np.sort(self.xyzarr[:,1])
        z_pl=np.sort(self.xyzarr[:,2])

        if self.args.plane1:
            #plane1 grid
            xx1, yy1 = np.meshgrid((x_pl[0],x_pl[-1]),(y_pl[0],y_pl[-1]))
            if args.plane2: #??
                #plane2 grid
                xx2, yy2 = np.meshgrid((x_pl[0],x_pl[-1]),(y_pl[0],y_pl[-1]))
            #plane 1 d
            d1 = -c1.dot(n1)
            
            if args.plane2:
                #plane 2 d
                d2 = -c2.dot(n2)
            #plane 1 equation
            z1 = (-n1[0] * xx1 - n1[1] * yy1 - d1) * 1. /n1[2]
            if args.plane2:
                #plane 2 equation
                z2 = (-n2[0] * xx2 - n2[1] * yy2 - d2) * 1. /n2[2]
            #plot plane 1
            surf = ax.plot_surface(xx1, yy1, z1, color='blue', alpha=0.3, label='Plane 1')
            
            if args.plane2:
                #plot plane 2
                surf = ax.plot_surface(xx2, yy2, z2, color='red', alpha=0.3, label='Plane 2')
        
        # show arrows representing the xyz-axes, starting from 0,0,0
        if self.args.showori:
            arrow_length =sum(abs(i) for i in ax.get_xlim())+sum(abs(i) for i in ax.get_ylim())+sum(abs(i) for i in ax.get_zlim())
            arrow_length = (arrow_length/3)*0.5
            if arrow_length > 3:
                arrow_length = 3
            ax.arrow3D(0,0,0, 0,0,arrow_length,
                                    mutation_scale=20,
                                    ec ='black',
                                    fc='red')
            ax.arrow3D(0,0,0, 0,arrow_length,0,
                                    mutation_scale=20,
                                    ec ='black',
                                    fc='green')
            ax.arrow3D(0,0,0, arrow_length,0,0,
                                    mutation_scale=20,
                                    ec ='black',
                                    fc='blue')
            ax.text(0, 0, arrow_length, 'z',color='red',fontsize=15)
            ax.text(0, arrow_length, 0, 'y',color='green',fontsize=15)
            ax.text(arrow_length, 0, 0, 'x',color='blue',fontsize=15)
            ax.scatter(0,0,0,s=50,color='black',alpha=0.8)
                
        finallize_figure_and_axis(fig, ax)
        plt.gca().set_zlim(z_pl[0],z_pl[-1]) #set z limits for plots, otherwise planes are sometimes very large

        #show the plot
        plt.show()
        #save the plot
        plot_filename=self.args.filename.split('.')[0]+'.png'
        plt.savefig(plot_filename)
        plt.clf()
        

    def record_data(self):
        data=TabData()
        setattr(data, 'filename', self.args.filename)
        setattr(data, 'No_atoms', self.xyz_df.shape[0])
        setattr(data, 'Sum_formula', ''.join(self.sum_formula_list))
        setattr(data, 'Excluded_atoms', re.sub(r'[^a-zA-Z0-9,]','',str(self.args.excludeAt)))
        setattr(data, 'Excluded_bonds', re.sub(r'[^a-zA-Z0-9,]','',str(self.args.excludeEl)))
        setattr(data, 'Included_contacts', re.sub(r'[^a-zA-Z0-9,]','',str(self.args.includeCon)))
        setattr(data, 'Covalent_radius', '{:.2f} %'.format(self.args.radius))
        setattr(data, 'Info_df', self.info_df)
        setattr(data, 'Bonds_boolean_mask', self.dist_df.is_bond == True)
        setattr(data, 'bonds_df_1', self.sel_dist[['atom1_idx', 'atom2_idx', 'A-B','distance_calc']])        
        setattr(data, 'angles_df_1', self.sel_angles[['atom1_idx', 'atom2_idx', 'atom3_idx', 'A-B-C','angle_calc']])
        return data

    def get_tab_data(self):
        self.get_inital_data()
        self.get_bonds_data()
        self.get_angles_data()
        if self.args.show or self.args.showbl or self.args.shownl:
                self.plot_the_molecule()
        data=self.record_data()
        return data.__dict__, self.plot_filename

if __name__=='__main__':
    pass
#     import os
#     home_path=r'C:\Users\edens\Documents\GitHub'
#     path_to_add=r'C:\Users\edens\Documents\GitHub\conformers_project\python_code\molecules\crest_conformers'
#     os.chdir(path_to_add)
# ##        os.chdir(home_path)
#     sys.path.insert(1, path_to_add)
# ##        xyz_filename_1='crest_conformers.xyz'
#     xyz_filename_1='crest_conformers_0.xyz'
#     # molecule_tab_1=TabDataAnalyzer(set_tab_parser(), xyz_filename_1, get_plot=True)
#     # data_dict, plot_filename=molecule_tab_1.get_tab_data()
# ##        print(plot_filename)
    
#     xyz_filename_2='crest_conformers_1.xyz'
# ##        molecule_tab_2=TabDataAnalyzer(set_tab_parser(), xyz_filename_2, get_plot=True) #bonds_by_ref=True, bonds_ref=data_dict.get('Bonds_boolean_mask')
# ##        data_dict, plot_filename=molecule_tab_2.get_tab_data()

    overlay_analyzer=OverlayAnalyzer(parser=set_overlay_parser(), xyz_filenames=[xyz_filename_1, xyz_filename_2], fit_mode='all',
##                                         atoms=['30', '20', '16', '0', '24', '8', '5', '33'],
                                      ) #'all'

        
