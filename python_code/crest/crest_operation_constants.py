from tools.general_constants import *
from tools.file_handlers import *
from tools.linux_operations import *

class CrestOutputFilenames(Enum):
    ENERGIES='crest.energies'
    CONFORMERS='crest_conformers.xyz'
    ROTAMERS='crest_rotamers.xyz'

class CrestConstants(Enum):
    MOLECULE_XTB_OUTPUT='_xtb_output.zip'
    CREST_OUTPUT_NAME='crest_conformers.xyz'
    SINGLE_CONFORMER_XYZ='_single_conformer.xyz'
    FINAL_OUTPUT_NAME='_crest_conformer_results.csv'

if __name__=='__main__':
    print('script loaded successfully')
