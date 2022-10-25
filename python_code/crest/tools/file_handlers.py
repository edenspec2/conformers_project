"""This moudle provides access to utlity functions for file processing"""
import numpy as np
import pandas as pd

from enum import Enum

class FileExtensions(Enum):
    """
    Hold commonly used file extensions
    """
    SMI='.smi'
    XYZ='.xyz'
    CSV='.csv'
    ZIP='.zip'
    PPT='.ppt'
    CIF='.cif'
    MOL='.mol'
    PDB='.pdb'

def get_filename_list(file_extension):
    """
    The function gets a file extension as input and returns a list of all files in the working directory
    ----------
    Parameters
    ----------
    file_extension : str.
        The wanted file extension like '.csv' or '.ppt'
    -------
    Returns
    -------
    list
        A list of all files in the working directory with the chosen extension 
    --------
    Examples
    --------
    from os import listdir

    all_files_in_dir=listdir()
    print(all_files_in_dir)
        ['0_1106253-mod-mod.xyz', '0_1106253-mod.xyz', '1106253.cif', '1109098.cif', '1_1106253-mod.xyz', 'centered_0_BASCIH.xyz', 'cif_handler.py']
        
    xyz_files_in_dir=get_filename_list('.xyz')
    print(xyz_files_in_dir)
        ['0_1106253-mod-mod.xyz', '0_1106253-mod.xyz', '1_1106253-mod.xyz', 'centered_0_BASCIH.xyz']
    
    cif_files_in_dir=get_filename_list('.cif')
    print(cif_files_in_dir)
        ['1106253.cif', '1109098.cif']    
    """
    from os import listdir
    return [filename for filename in listdir() if filename.endswith(file_extension)]    

def get_file_lines(filename, encoding=None):
    """
    The function gets a file name/file location and return the raw text lines of the file as a list
    ----------
    Parameters
    ----------
    filename : str.
        The wanted file name/location like '0_1106253-mod.xyz'

    encoding : str. default None
        Encoding format like 'utf-8'. Argument passed to the native python function 'open'
    -------
    Returns
    -------
    list
        A list of the raw lines of the file
    --------
    Examples
    --------
    filename='0_1106253-mod.xyz'
    file_lines=get_file_lines(filename)
    print(file_lines)
        ['40\n', 'comment_line\n', 'C  1.56000000 0.24705000 -0.67600000\n' ...]
    """
    with open(filename, 'r', encoding=encoding) as my_file:
        lines=my_file.readlines()
    return lines

def get_file_striped_lines(filename, encoding=None):
    """
    The function gets a file name/file location and return the stripped text lines of the file as a list
    ----------
    Parameters
    ----------
    filename : str.
        The wanted file name/location like '0_1106253-mod.xyz'

    encoding : str. default None
        Encoding format like 'utf-8'. Argument passed to the native python function 'open'
    -------
    Returns
    -------
    list
        A list of the striped lines of the file
    --------
    Examples
    --------
    filename='0_1106253-mod.xyz'
    file_lines=get_file_lines(filename)
    print(file_lines)
        ['40\n', 'comment_line\n', 'C  1.56000000 0.24705000 -0.67600000\n' ...]

    striped_file_lines=get_file_lines(filename)
    print(striped_file_lines)
        ['40', 'comment_line', 'C  1.56000000 0.24705000 -0.67600000' ...]
    """
    lines=get_file_lines(filename, encoding)
    strip_lines=[line.strip().rstrip('\n') for line in lines]
    return strip_lines


def change_filename(old_filename, new_filename):
    """
    The function gets a file name and new file name and changes it. A wrapper for the os.rename function
    ----------
    Parameters
    ----------
    old_filename : str.
        The file name you want to change like '0_1106253-mod.xyz'

    new_filename : str.
        The wanted file name you want in the end '1106253.xyz'
    -------
    Returns
    -------
    None
    --------
    Examples
    --------
    from os import listdir

    all_files_in_dir=listdir()
    print(all_files_in_dir)
        ['0_1106253-mod-mod.xyz', '0_1106253-mod.xyz', '1106253.cif', '1109098.cif', '1_1106253-mod.xyz', 'centered_0_BASCIH.xyz', 'cif_handler.py']
    
    old_filename='0_1106253-mod.xyz'
    new_filename='1106253.xyz'
    change_filename(old_filename, new_filename)
    
    all_files_in_dir=listdir()
    print(all_files_in_dir)
        ['0_1106253-mod-mod.xyz', '1106253.xyz', '1106253.cif', '1109098.cif', '1_1106253-mod.xyz', 'centered_0_BASCIH.xyz', 'cif_handler.py']
    """
    from os import rename
    rename(old_filename, new_filename)

def exchange_file_extension(filename, new_extension):
    """
    The function gets a file name and new file extension and changes the file name extension.
    ----------
    Parameters
    ----------
    filename : str.
        The file name you want to change like '0_1106253-mod.xyz'

    new_extension : str.
        The new file extension you want in the end '.cif' or 'csv'
    -------
    Returns
    -------
    String
        the new filename with the new file extension
    --------
    Examples
    --------
    filename='0_1106253-mod.xyz'
    new_extension='mol'
    new_filename=exchange_filename(old_filename, new_filename)
    print(new_filename)
        '0_1106253-mod.mol'
    """
    if not new_extension.startswith('.'):
        new_extension='.'+new_extension
    return filename.split('.')[0]+new_extension

def save_to_new_zipfile(zip_filename, filenames_to_zip):
    """
    The function gets a zip file name and list of files to zip, and the function creates a zip that contains all the files.
    ----------
    Parameters
    ----------
    zip_filename : str.
        The zip file name you want to create like 'my_zipfile.zip'

    filenames_to_zip : iterable
        An iterable containing file names/locations that you want to zip like ['1106253.cif', '1109098.cif']
    -------
    Returns
    -------
    None
    --------
    Examples
    --------
    zip_filename='my_zipfile.zip'
    filenames_to_zip=['1106253.cif', '1109098.cif'] 
    save_to_new_zipfile(zip_filename, filenames_to_zip)
    """
    from zipfile import ZipFile, ZIP_DEFLATED
    with ZipFile(zip_filename, 'w', ZIP_DEFLATED) as zip_file:
        for filename in filenames_to_zip:
            zip_file.write(filename)   

def get_file_from_zipfile(zip_filename, filename_to_unzip):
    """
    The function gets a zip file name and a name of file inside the zip file, and the function extract the requested file,
    and return it's path
    ----------
    Parameters
    ----------
    zip_filename : str.
        The zip file name you want to extract from like 'my_zipfile.zip'

    filename_to_unzip : str.
        A name of file inside the zipfile that you want to extract like '1106253.cif'
    -------
    Returns
    -------
    String
        the path to the file that was extracted from the zipfile
    --------
    Examples
    --------
    zip_filename='my_zipfile.zip'
    filename_to_unzip='1106253.cif' 
    filepath=get_file_from_zipfile(zip_filename, filename_to_unzip)
    print(filepath)
        'C:\...\1106253.cif'
    """
    from zipfile import ZipFile, ZIP_DEFLATED
    with ZipFile(zip_filename, 'r', ZIP_DEFLATED) as zip_file:
        output_filepath=zip_file.extract(filename_to_unzip)
    return output_filepath

def save_single_xyz_file(symbols, coordinates, output_filename, comment_line='comment_line'):
    """
    The function gets a the atom symbols and coordinates, and the function creates a xyz file from it.
    ----------
    Parameters
    ----------
    symbols : iterable.
        Array of atom types, each item is a string

    coordinates : iterable.
        Array of cartesian coordinates of all atoms, each item is a float

    output_filename : str.
        The name given to the output xyz file
        
    comment_line : str. default 'comment_line'
        A line recorded into the second line of the xyz file
    -------
    Returns
    -------
    None
    --------
    Examples
    --------
    symbols=['O', 'H', 'H']
    coordinates=[[0.00000, 0.00000, 0.11779], [0.00000, 0.75545, -0.47116], [0.00000, -0.75545, -0.47116]]
    output_filename='water.xyz'
    comment_line='Water molecule'
    save_single_xyz_file(symbols, coordinates, output_filename, comment_line)
    --water.xyz--
    3
    Water molecule
    O          0.00000        0.00000        0.11779
    H          0.00000        0.75545       -0.47116
    H          0.00000       -0.75545       -0.47116
    --end_of_file--
    """
    with open(output_filename, 'w') as my_file:
        my_file.write('{}\n'.format(len(symbols)))
        my_file.write(comment_line+'\n')
        for loop_index, (x, y, z) in enumerate(coordinates):
            my_file.write('{}\t {:.5f}\t {:.5f}\t {:.5f}\t\n'.format(symbols[loop_index], x, y, z))

def load_single_xyz_file(xyz_filename):
    """
    The function gets a xyz filename, and the function return the atoms symbols and atom coordinates.
    ----------
    Parameters
    ----------
    xyz_filename : str.
        The xyz filename you want to load
    -------
    Returns
    -------
    symbols : iterable.
        Array of atom types, each item is a string

    coordinates : iterable.
        Array of cartesian coordinates of all atoms, each item is a float
    --------
    Examples
    --------
    --water.xyz--
    3
    Water molecule
    O          0.00000        0.00000        0.11779
    H          0.00000        0.75545       -0.47116
    H          0.00000       -0.75545       -0.47116
    --end_of_file--
    xyz_filename='water.xyz'
    symbols, coordinates=load_single_xyz_file(xyz_filename)
    print(symbols)
        ['O', 'H', 'H']
    print(coordinates)
        [[0.00000, 0.00000, 0.11779], [0.00000, 0.75545, -0.47116], [0.00000, -0.75545, -0.47116]]
    """
    xyz_file=np.genfromtxt(fname=xyz_filename, skip_header=2, dtype='unicode')
    symbols=xyz_file[:,0]
    coordinates=(xyz_file[:,1:])
    coordinates=coordinates.astype(float)
    return symbols, coordinates    

def get_xyz_df_from_file(xyz_filename):
    """
    The function gets a xyz filename, and the function return dataframe with xyz data.
    ----------
    Parameters
    ----------
    xyz_filename : str.
        The xyz filename you want to load
    -------
    Returns
    -------
    DataFrame
        A dataframe with raw data of the xyz file. the columns are: ["element", "x", "y", "z"]
    --------
    Examples
    --------
    --water.xyz--
    3
    Water molecule
    O          0.00000        0.00000        0.11779
    H          0.00000        0.75545       -0.47116
    H          0.00000       -0.75545       -0.47116
    --end_of_file--
    xyz_filename='water.xyz'
    xyz_df=get_xyz_df_from_file(xyz_filename)
    print(xyz_df)
          element    x        y        z
        0       O  0.0  0.00000  0.11779
        1       H  0.0  0.75545 -0.47116
        2       H  0.0 -0.75545 -0.47116
    """
    xyz_df=pd.read_csv(xyz_filename,
                       delim_whitespace=True,
                       skiprows=2,
                       names=["element", "x", "y", "z"],
                       error_bad_lines='skip')
    return xyz_df

def extract_single_xyz_file_from_ensamble(ensamble_xyz_file, file_index, output_filename):
    """
    The function gets a xyz ensamble filename, that contains more that one molecule, and saves a xyz file of
    a single molecule according to the give file index.
    The function assumes that every molecule in the xyz file has the same number of atoms (like xyz file of conformer/trj.)
    ----------
    Parameters
    ----------
    ensamble_xyz_filename : str.
        The ensamble xyz filename you want to load

    file_index : int.
        The number of the molecule you want to isolate from the ensamble

    output_filename : str.
        The name given to the output xyz file
    -------
    Returns
    -------
    None
    --------
    Examples
    --------
    --crest_conformers.xyz--
    24
            -43.36823318
     O         -0.4828638664       -1.9500973380       -1.1282222019
     C         -1.6027687881       -1.5929633302       -0.3599335790
     C         -1.4640162543       -0.3073244166        0.4802374845
    ..
    ..
    ..
    24
            -43.36609350
     O         -0.5149278744       -1.9326809371       -1.1484212156
     C         -1.6122653533       -1.6085025372       -0.3360185428
     C         -1.4658513491       -0.3132080922        0.4879324683
    ..
    ..
    ..
    24
            -43.36566800
     O         -3.1347625091       -1.3858833466       -0.5243128010
     C         -1.8246055953       -1.7947731283       -0.2681858912
     C         -0.9659366197       -0.6584341364        0.3045489942
    --end_of_file--
    ensamble_xyz_file='crest_conformers.xyz'
    file_index=2
    output_filename='third_conformer.xyz'
    extract_single_xyz_file_from_ensamble(ensamble_xyz_file, file_index, output_filename)
    --third_conformer.xyz--
    24
            -43.36566800
     O         -3.1347625091       -1.3858833466       -0.5243128010
     C         -1.8246055953       -1.7947731283       -0.2681858912
     C         -0.9659366197       -0.6584341364        0.3045489942
    --end_of_file--    
    """
    file_lines=get_file_lines(ensamble_xyz_file) #CrestConstants.CREST_OUTPUT_NAME.value
    number_of_single_molecule_lines=int(file_lines[0].strip())+2
    with open(output_filename, 'w') as f:
        f.writelines(file_lines[file_index*number_of_single_molecule_lines:(file_index+1)*number_of_single_molecule_lines])


##def get_indiv_xyz_from_cif(diffpy_structure):
##    number_of_molecules=get_number_of_molecules_cif_labels(diffpy_structure.label)
##    total_size=len(diffpy_structure.label)
##    molecule_atoms, molecule_locations=[], []
##    for molecule_number in range(number_of_molecules):
##        molecule_indices=slice(molecule_number, total_size, number_of_molecules)
##        cif_atoms=diffpy_structure.label[molecule_indices]
##        molecule_atoms.append(np.array([re.sub(r'[^A-Z]', '', atom) for atom in cif_atoms]))
##        molecule_locations.append(diffpy_structure.xyz_cartn[molecule_indices])
##    return molecule_atoms, molecule_locations # lists of the approiate items
###           molecule_atoms=[mol_1_atoms, mol_2_atoms, ...]

if __name__=='__main__':
    print('script loaded successfully')
