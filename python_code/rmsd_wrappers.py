"""
This moudle provides access to utlity functions for molecule alignment and rmsd calculation
rmsd between two molecules is defines as follows:
            sum((atom_i_mol_1_location-atom_i_mol_2_location)**2)
RMSD=sqrt(---------------------------------------------------------)
                                number_of_atoms
All functions are wrappers for the functions in the rmsd moudle
rmsd moudle: http://github.com/charnley/rmsd
Kabsch algorithm: Kabsch W., 1976, A solution for the best rotation to relate two sets of vectors,
                  Acta Crystallographica, A32:922-923, doi: http://dx.doi.org/10.1107/S0567739476001873
Hungarian algorithm : https://en.wikipedia.org/wiki/Hungarian_algorithm

"""
import sys
path_to_add=r'C:\Users\\itaro\OneDrive\Documents\GitHub'
sys.path.insert(0, path_to_add)

from tools.general_constants import *
from tools.file_handlers import *

class RMSDData():
    """
    Holds rmsd data for easy usage
    """
    pass

class RMSDConstants(Enum):
    """
    Holds constants related to rmsd calculation and data representation
    """
    ROUND_DIGITS=3

def rmsd_center_coords(coordinates):
    """
    The function gets a np.array with molecule coordinates, and the function translates the coordinates
    so that center of the molecule will sit on [0,0] or [0,0,0]
    ----------
    Parameters
    ----------
    coordinates : np.array.
        Array of cartesian coordinates of atoms, each item is a float
    -------
    Returns
    -------
    np.array
        Array of centered cartesian coordinates of atoms, each item is a float 
    --------
    Examples
    --------
    coordinates=np.array([[0.00000, 0.00000, 0.11779], [0.00000, 0.75545, -0.47116], [0.00000, -0.75545, -0.47116]])
    centered_coordinates=rmsd_center_coords(coordinates)
    print(centered_coordinates)
        [[ 0.          0.          0.39263333]
         [ 0.          0.75545    -0.19631667]
         [ 0.         -0.75545    -0.19631667]]
    """
    from rmsd import centroid
    new_coordinates=coordinates-centroid(coordinates)
    return new_coordinates

def rmsd_translate_molecules(coordinates_1, coordinates_2):
    """
    The function gets two coordinates np.arrays, and the function translates the coordinates
    so that center of the molecules will sit on [0,0] or [0,0,0]
    ----------
    Parameters
    ----------
    coordinates_1 : np.array.
        Array of cartesian coordinates of atoms in molecule_1, each item is a float

    coordinates_2 : np.array.
        Array of cartesian coordinates of atoms in molecule_2, each item is a float
    -------
    Returns
    -------
    np.array
        Array of centered cartesian coordinates of atoms in molecule_1, each item is a float

    np.array
        Array of centered cartesian coordinates of atoms in molecule_2, each item is a float 
    --------
    Examples
    --------
    coordinates_1=np.array([[0.00000, 0.00000, 0.11779], [0.00000, 0.75545, -0.47116], [0.00000, -0.75545, -0.47116]])
    coordinates_2=np.array([[1.00000, 1.00000, 1.11779], [1.00000, 1.75545, 0.52884], [1.00000, 0.24455, 0.52884]])
    centered_coordinates_1, centered_coordinates_2=rmsd_translate_molecules(coordinates_1, coordinates_2)
    print(centered_coordinates_1)
        [[ 0.          0.          0.39263333]
         [ 0.          0.75545    -0.19631667]
         [ 0.         -0.75545    -0.19631667]]
    print(centered_coordinates_2)
        [[ 0.00000000e+00  1.11022302e-16  3.92633333e-01]
         [ 0.00000000e+00  7.55450000e-01 -1.96316667e-01]
         [ 0.00000000e+00 -7.55450000e-01 -1.96316667e-01]]
    """
    return rmsd_center_coords(coordinates_1), rmsd_center_coords(coordinates_2)

def rmsd_rotate_molecules(coordinates_1, coordinates_2):
    """
    The function gets two coordinates np.arrays, and the function rotates the one of them
    so that the molecules will sit on each other as much as possible
    ----------
    Parameters
    ----------
    coordinates_1 : np.array.
        Array of cartesian coordinates of atoms in molecule_1, each item is a float

    coordinates_2 : np.array.
        Array of cartesian coordinates of atoms in molecule_2, each item is a float
    -------
    Returns
    -------
    np.array
        Array of rotated cartesian coordinates of atoms in molecule_1, each item is a float

    np.array
        Array of (unchanged) cartesian coordinates of atoms in molecule_2, each item is a float 
    --------
    Examples
    --------
    coordinates_1=np.array([[0.00000, 0.00000, 0.11779], [0.00000, 0.75545, -0.47116], [0.00000, -0.75545, -0.47116]])
    coordinates_2=np.array([[1.00000, 1.00000, 1.11779], [1.00000, 1.75545, 0.52884], [1.00000, 0.24455, 0.52884]])
    rotated_coordinates_1, *_=rmsd_rotate_molecules(coordinates_1, coordinates_2)
    print(rotated_coordinates_1)
        [[-0.09981108 -0.04406901 -0.04438642]
         [ 0.14099132  0.87686183  0.06269944]
         [ 0.65749734 -0.52430972  0.29239189]]
    """
    from rmsd import kabsch
    U=kabsch(coordinates_1,coordinates_2)
    rotated_coordinates_1=np.dot(coordinates_1,U)
    return rotated_coordinates_1, coordinates_2

def rmsd_align_xyz(coordinates_1, coordinates_2):
    """
    The function gets two coordinates np.arrays, and the function aligns both togather
    by centering both coordinates and then rotating molecule_1 so that RMSD between the
    two molecules is minimized using the kabsch algorithm
                sum((atom_i_mol_1_location-atom_i_mol_2_location)**2)
    RMSD=sqrt(---------------------------------------------------------)
                                    number_of_atoms
    ----------
    Parameters
    ----------
    coordinates_1 : np.array.
        Array of cartesian coordinates of atoms in molecule_1, each item is a float

    coordinates_2 : np.array.
        Array of cartesian coordinates of atoms in molecule_2, each item is a float
    -------
    Returns
    -------
    np.array
        Array of rotated+centered cartesian coordinates of atoms in molecule_1, each item is a float

    np.array
        Array of centered cartesian coordinates of atoms in molecule_2, each item is a float 
    --------
    Examples
    --------
    coordinates_1=np.array([[0.00000, 0.00000, 0.11779], [0.00000, 0.75545, -0.47116], [0.00000, -0.75545, -0.47116]])
    coordinates_2=np.array([[1.00000, 1.00000, 1.11779], [1.00000, 1.75545, 0.52884], [1.00000, 0.24455, 0.52884]])
    rotated_coordinates_1, translated_coordinates_2=rmsd_align_xyz(coordinates_1, coordinates_2)
    print(rotated_coordinates_1)
        [[ 0.00000000e+00  8.11228564e-18  3.92633333e-01]
         [ 0.00000000e+00  7.55450000e-01 -1.96316667e-01]
         [ 0.00000000e+00 -7.55450000e-01 -1.96316667e-01]]
    print(translated_coordinates_2)
        [[ 0.00000000e+00  1.11022302e-16  3.92633333e-01]
         [ 0.00000000e+00  7.55450000e-01 -1.96316667e-01]
         [ 0.00000000e+00 -7.55450000e-01 -1.96316667e-01]]
    """
    translated_coordinates_1, translated_coordinates_2=rmsd_translate_molecules(coordinates_1, coordinates_2)
    rotated_coordinates_1, *_=rmsd_rotate_molecules(translated_coordinates_1, translated_coordinates_2)
    return rotated_coordinates_1, translated_coordinates_2

def rmsd_renumbering_coordinates(symbols_1, symbols_2, coordinates_1, coordinates_2, renumber_method='hungarian'):
    """
    The function gets the atom symbols and coordinates of two molecules, and the function renumbers
    symbols_2 and coordinates_2 to match the order of symbols_1 and coordinates_1
    ----------
    Parameters
    ----------
    symbols_1 : np.array.
        Array of element symbols of atoms in molecule_1, each item is a string

    symbols_2 : np-array.
        Array of element symbols of atoms in molecule_2, each item is a string
        
    coordinates_1 : np.array.
        Array of cartesian coordinates of atoms in molecule_1, each item is a float

    coordinates_2 : np.array.
        Array of cartesian coordinates of atoms in molecule_2, each item is a float

    renumber_method : str. default 'hungarian'.
        The method used to renumber the molecules. there are two possible options:
        'hungarian' - fast and efficent, could miss sometimes
        'brute' - very slow, check all possible renumbering options
    -------
    Returns
    -------
    np.array
        Array of renumbered element symbols of atoms in molecule_2, each item is a string

    np.array
        Array of renumbered cartesian coordinates of atoms in molecule_2, each item is a float 
    --------
    Examples
    --------
    symbols_1=np.array(['C', 'C', 'H', 'C', ...])
    coordinates_1=np.array([[ 1.56     0.24705 -0.676  ],
                            [-0.466   -2.96895  0.585  ],
                            [-1.393   -0.53795  3.286  ],
                            [ 2.909   -0.05495 -0.873  ],
                            ...
                            )
    symbols_2=np.array(['C', 'C', 'C', 'H', ...])
    coordinates_2=np.array([ 0.466 -2.969 -0.585],
                            [ 1.719  3.471 -0.752],
                            [ 1.393 -0.538 -3.286],
                            [-2.909 -0.055  0.873],
                            ...
                            )
    renumber_method='hungarian'
    renumbered_symbols_2, renumbered_coordinates_2=rmsd_renumbering_coordinates(symbols_1, symbols_2, coordinates_1, coordinates_2, renumber_method)
    print(renumbered_symbols_2)
        ['C', 'C', 'H', 'C', ...]
    print(renumbered_coordinates_2)
        [[ 1.192,  0.743, -0.983],
        [ 0.466, -2.969, -0.585],
        [-1.373,  0.526,  3.269],
        [ 2.461,  1.261, -1.256],
        ...
        ]
    """
    from rmsd import reorder_hungarian, reorder_brute
    if renumber_method=='hungarian':
        new_view=reorder_hungarian(symbols_1, symbols_2, coordinates_1, coordinates_2)
    if renumber_method=='brute':
        new_view=reorder_brute(symbols_1, symbols_2, coordinates_1, coordinates_2)
    return symbols_2[new_view], coordinates_2[new_view]

def rmsd_renumber_and_align(symbols_1, symbols_2, coordinates_1, coordinates_2, renumber_method='hungarian'):
    """
    The function gets the atom symbols and coordinates of two molecules, and the function renumbers (so that the order matches)
    and aligns the two molecules so that RMSD between them is minimized using the kabsch algorithm
                sum((atom_i_mol_1_location-atom_i_mol_2_location)**2)
    RMSD=sqrt(---------------------------------------------------------)
                                    number_of_atoms
    ----------
    Parameters
    ----------
    symbols_1 : np.array.
        Array of element symbols of atoms in molecule_1, each item is a string
        
    symbols_2 : np-array.
        Array of element symbols of atoms in molecule_2, each item is a string
        
    coordinates_1 : np.array.
        Array of cartesian coordinates of atoms in molecule_1, each item is a float

    coordinates_2 : np.array.
        Array of cartesian coordinates of atoms in molecule_2, each item is a float

    renumber_method : str. default 'hungarian'.
        The method used to renumber the molecules. there are two possible options:
        'hungarian' - fast and efficent, could miss sometimes
        'brute' - very slow, check all possible renumbering options
    -------
    Returns
    -------
    np.array
        Array of the altered element symbols of atoms in molecule_1, each item is a string

    np.array
        Array of the altered element symbols of atoms in molecule_2, each item is a string

    np.array
        Array of the altered cartesian coordinates of atoms in molecule_1, each item is a float

    np.array
        Array of the altered cartesian coordinates of atoms in molecule_2, each item is a float 
    --------
    Examples
    --------
    symbols_1=np.array(['C', 'C', 'H', 'C', ...])
    coordinates_1=np.array([[ 1.56     0.24705 -0.676  ],
                            [-0.466   -2.96895  0.585  ],
                            [-1.393   -0.53795  3.286  ],
                            [ 2.909   -0.05495 -0.873  ],
                            ...
                            )
    symbols_2=np.array(['C', 'C', 'C', 'H', ...])
    coordinates_2=np.array([ 0.466 -2.969 -0.585],
                            [ 1.719  3.471 -0.752],
                            [ 1.393 -0.538 -3.286],
                            [-2.909 -0.055  0.873],
                            ...
                            )
    renumber_method='hungarian'
    new_symbols_1, new_symbols_2, new_coordinates_1, new_coordinates_2=rmsd_renumber_and_align(symbols_1, symbols_2, coordinates_1, coordinates_2, renumber_method)
    print(new_symbols_1)
        ['C', 'C', 'H', 'C', ...]
    print(new_symbols_2)
        ['C', 'C', 'H', 'C', ...]
    print(new_coordinates_1)
        [[ 1.50555389,  0.20641808, -0.80144169],
        [ 0.16420583, -2.27987675,  2.03696921],
        [-1.43785849,  0.91156516,  3.18264781],
        [ 2.89737211,  0.15179577, -0.89979019],
        ...
        ]
    print(new_coordinates_2)
        [[ 1.192,  0.743, -0.983],
        [ 0.466, -2.969, -0.585],
        [-1.373,  0.526,  3.269],
        [ 2.461,  1.261, -1.256],
        ...
        ]
    """
    renumbered_symbols_2, renumbered_coordinates_2=rmsd_renumbering_coordinates(symbols_1, symbols_2, coordinates_1, coordinates_2, renumber_method)
    aligned_coordinates_1, aligned_coordinates_2=rmsd_align_xyz(coordinates_1, renumbered_coordinates_2)
    return symbols_1, renumbered_symbols_2, aligned_coordinates_1, aligned_coordinates_2

def rmsd_get_RMSD_score(coordinates_1, coordinates_2):
    """
    The function gets two coordinates np.arrays, and the function calculates the RMSD between
    the two coordinates
                sum((atom_i_mol_1_location-atom_i_mol_2_location)**2)
    RMSD=sqrt(---------------------------------------------------------)
                                    number_of_atoms
    ----------
    Parameters
    ----------
    coordinates_1 : np.array.
        Array of cartesian coordinates of atoms in molecule_1, each item is a float

    coordinates_2 : np.array.
        Array of cartesian coordinates of atoms in molecule_2, each item is a float
    -------
    Returns
    -------
    float
        RMSD score between the two molecules 
    --------
    Examples
    --------
    coordinates_1=np.array([[0.00000, 0.00000, 0.11779], [0.00000, 0.75545, -0.47116], [0.00000, -0.75545, -0.47116]])
    coordinates_2=np.array([[1.00000, 1.00000, 1.11779], [1.00000, 1.75545, 0.52884], [1.00000, 0.24455, 0.52884]])
    rmsd_score=rmsd_get_RMSD_score(coordinates_1, coordinates_2)
    print(rmsd_score)
        1.732
    """
    from rmsd import rmsd
    return round(rmsd(coordinates_1,coordinates_2), RMSDConstants.ROUND_DIGITS.value)

def rmsd_score_two_molecules(symbols_1, symbols_2, coordinates_1, coordinates_2, return_altered_stats=False, renumber_method='hungarian'):
    """
    The function gets the atom symbols and coordinates of two molecules, and the function renumbers (so that the order matches)
    and aligns the two molecules so that RMSD between them is minimized using the kabsch algorithm and returns the RMSD score.
                sum((atom_i_mol_1_location-atom_i_mol_2_location)**2)
    RMSD=sqrt(---------------------------------------------------------)
                                    number_of_atoms
    ----------
    Parameters
    ----------
    symbols_1 : np.array.
        Array of element symbols of atoms in molecule_1, each item is a string

    symbols_2 : np-array.
        Array of element symbols of atoms in molecule_2, each item is a string
        
    coordinates_1 : np.array.
        Array of cartesian coordinates of atoms in molecule_1, each item is a float

    coordinates_2 : np.array.
        Array of cartesian coordinates of atoms in molecule_2, each item is a float

    return_altered_stats : boolean. default False.
        Determine the output of the function.
        False - Only returns the RMSD score
        True - Returns the RMSD score, symbols of both molecules and coordinates of both molecules

    renumber_method : str. default 'hungarian'.
        The method used to renumber the molecules. there are two possible options:
        'hungarian' - fast and efficent, could miss sometimes
        'brute' - very slow, check all possible renumbering options
    -------
    Returns
    -------
    float
        The RMSD score of the two molecules
        
    np.array (if return_altered_stats==True)
        Array of the altered element symbols of atoms in molecule_1, each item is a string

    np.array (if return_altered_stats==True)
        Array of the altered element symbols of atoms in molecule_2, each item is a string

    np.array (if return_altered_stats==True)
        Array of the altered cartesian coordinates of atoms in molecule_1, each item is a float

    np.array (if return_altered_stats==True)
        Array of the altered cartesian coordinates of atoms in molecule_2, each item is a float 
    --------
    Examples
    --------
    symbols_1=np.array(['C', 'C', 'H', 'C', ...])
    coordinates_1=np.array([[ 1.56     0.24705 -0.676  ],
                            [-0.466   -2.96895  0.585  ],
                            [-1.393   -0.53795  3.286  ],
                            [ 2.909   -0.05495 -0.873  ],
                            ...
                            )
    symbols_2=np.array(['C', 'C', 'C', 'H', ...])
    coordinates_2=np.array([ 0.466 -2.969 -0.585],
                            [ 1.719  3.471 -0.752],
                            [ 1.393 -0.538 -3.286],
                            [-2.909 -0.055  0.873],
                            ...
                            )
    renumber_method='hungarian'
    return_altered_stats=False
    rmsd_score=rmsd_score_two_molecules(symbols_1, symbols_2, coordinates_1, coordinates_2, return_altered_stats, renumber_method)
    print(rmsd_score)
        1.925
    return_altered_stats=True
    rmsd_score, new_symbols_1, new_symbols_2, new_coordinates_1, new_coordinates_2=rmsd_score_two_molecules(symbols_1, symbols_2, coordinates_1, coordinates_2, return_altered_stats, renumber_method)
    print(new_symbols_1)
        ['C', 'C', 'H', 'C', ...]
    print(new_symbols_2)
        ['C', 'C', 'H', 'C', ...]
    print(new_coordinates_1)
        [[ 1.50555389,  0.20641808, -0.80144169],
        [ 0.16420583, -2.27987675,  2.03696921],
        [-1.43785849,  0.91156516,  3.18264781],
        [ 2.89737211,  0.15179577, -0.89979019],
        ...
        ]
    print(new_coordinates_2)
        [[ 1.192,  0.743, -0.983],
        [ 0.466, -2.969, -0.585],
        [-1.373,  0.526,  3.269],
        [ 2.461,  1.261, -1.256],
        ...
        ]
    """
    aligned_symbols_1, aligned_symbols_2, aligned_coordinates_1, aligned_coordinates_2=rmsd_renumber_and_align(symbols_1, symbols_2, coordinates_1, coordinates_2, renumber_method)
    rmsd_score=rmsd_get_RMSD_score(aligned_coordinates_1, aligned_coordinates_2)
    if return_altered_stats:
        return rmsd_score, aligned_symbols_1, aligned_symbols_2, aligned_coordinates_1, aligned_coordinates_2
    else:
        return rmsd_score

def rmsd_score_two_molecules_from_xyz(xyz_filename1, xyz_filename2, return_altered_stats=False, renumber_method='hungarian'):
    """
    The function gets two xyz filenames, and the function renumbers (so that the order matches)
    and aligns the two molecules so that RMSD between them is minimized using the kabsch algorithm and returns the RMSD score.
                sum((atom_i_mol_1_location-atom_i_mol_2_location)**2)
    RMSD=sqrt(---------------------------------------------------------)
                                    number_of_atoms
    ----------
    Parameters
    ----------
    xyz_filename1 : str.
        XYZ filename of molecule_1

    xyz_filename2 : str.
        XYZ filename of molecule_2

    return_altered_stats : boolean. default False.
        Determine the output of the function.
        False - Only returns the RMSD score
        True - Returns the RMSD score, symbols of both molecules and coordinates of both molecules

    renumber_method : str. default 'hungarian'.
        The method used to renumber the molecules. there are two possible options:
        'hungarian' - fast and efficent, could miss sometimes
        'brute' - very slow, check all possible renumbering options
    -------
    Returns
    -------
    float
        The RMSD score of the two molecules
        
    np.array (if return_altered_stats==True)
        Array of the altered element symbols of atoms in molecule_1, each item is a string

    np.array (if return_altered_stats==True)
        Array of the altered element symbols of atoms in molecule_2, each item is a string

    np.array (if return_altered_stats==True)
        Array of the altered cartesian coordinates of atoms in molecule_1, each item is a float

    np.array (if return_altered_stats==True)
        Array of the altered cartesian coordinates of atoms in molecule_2, each item is a float 
    --------
    Examples
    --------
    --0_1106253.xyz--
    40
    comment_line
    C     1.56000000    0.24705000   -0.67600000
    C    -0.46600000   -2.96895000    0.58500000
    H    -1.39300000   -0.53795000    3.28600000
    C     2.90900000   -0.05495000   -0.87300000
    --end_of_file--
    --1_1106253.xyz--
    40
    comment_line
    H	 -3.139	 0.721	 1.585	
    C	 -0.544	 0.341	 1.604	
    C	 1.719	 -3.471	 -0.752	
    H	 -1.373	 -0.526	 3.269
    --end_of_file--
    xyz_filename1='0_1106253.xyz'
    xyz_filename2='1_1106253.xyz' 
    renumber_method='hungarian'
    return_altered_stats=False
    rmsd_score=rmsd_score_two_molecules_from_xyz(xyz_filename1, xyz_filename2, return_altered_stats, renumber_method)
    print(rmsd_score)
        1.925
    return_altered_stats=True
    rmsd_score, new_symbols_1, new_symbols_2, new_coordinates_1, new_coordinates_2=rmsd_score_two_molecules_from_xyz(xyz_filename1, xyz_filename2, return_altered_stats, renumber_method)
    print(new_symbols_1)
        ['C', 'C', 'H', 'C', ...]
    print(new_symbols_2)
        ['C', 'C', 'H', 'C', ...]
    print(new_coordinates_1)
        [[ 1.50555389,  0.20641808, -0.80144169],
        [ 0.16420583, -2.27987675,  2.03696921],
        [-1.43785849,  0.91156516,  3.18264781],
        [ 2.89737211,  0.15179577, -0.89979019],
        ...
        ]
    print(new_coordinates_2)
        [[ 1.192,  0.743, -0.983],
        [ 0.466, -2.969, -0.585],
        [-1.373,  0.526,  3.269],
        [ 2.461,  1.261, -1.256],
        ...
        ]
    """
    symbols_1, coordinates_1=load_single_xyz_file(xyz_filename1)
    symbols_2, coordinates_2=load_single_xyz_file(xyz_filename2)
    return rmsd_score_two_molecules(symbols_1, symbols_2, coordinates_1, coordinates_2, return_altered_stats, renumber_method)

def rmsd_score_two_aligned_molecules_from_df(xyz_df_1, xyz_df_2):
    coordinates_1=xyz_df_1[['x', 'y', 'z']].to_numpy()
    coordinates_2=xyz_df_2[['x', 'y', 'z']].to_numpy()
    return rmsd_get_RMSD_score(coordinates_1, coordinates_2)

def rmsd_score_two_aligned_molecules_from_xyz(xyz_filename_1, xyz_filename_2):
    _, coordinates_1=load_single_xyz_file(xyz_filename_1)
    _, coordinates_2=load_single_xyz_file(xyz_filename_2)
    return rmsd_get_RMSD_score(coordinates_1, coordinates_2)

if __name__=='__main__':
    xyz_filename1='0_1106253-mod.xyz'
    xyz_filename2='1_1106253.xyz' 
    rmsd_score=rmsd_score_two_molecules_from_xyz(xyz_filename1, xyz_filename2, return_altered_stats=False)
    rmsd_score, aligned_symbols_1, aligned_symbols_2, aligned_coordinates_1, aligned_coordinates_2=rmsd_score_two_molecules_from_xyz(xyz_filename1, xyz_filename2, return_altered_stats=True)
    print(rmsd_score)
