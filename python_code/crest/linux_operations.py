import re
import os
from enum import Enum

class LinuxCommands(Enum):
    OBABEL='/gpfs0/gaus/users/itamarwa/transfer_to_itamar/build/bin/obabel '
    OBABEL_XYZ_SETTINGS_1=' -O '
    OBABEL_XYZ_SETTINGS_2=' --gen3d'
    XTB_INPUT_PREFIX='/gpfs0/gaus/projects/xtb-6.4.1/bin/xtb --input ' # xtb_di.inp 
    XTB_PREFIX='/gpfs0/gaus/projects/xtb-6.4.1/bin/xtb '
    XTB_SUFIX=' --ohess --dipole --pop'
    CREST_INPUT_PREFIX='/gpfs0/gaus/projects/crest --input '
    CREST_PREFIX='/gpfs0/gaus/projects/crest '
    CREST_SUFIX=' --xnam /gpfs0/gaus/projects/xtb-6.4.1/bin/xtb'
    COPY_COMMAND='cp ' # location_1 -> location_2
    MOVE_COMMAND='mv ' # location_1 -> location_2

def run_crest_calculation(xyz_filename, detailed_input=None):
    if detailed_input:
        os.system(LinuxCommands.CREST_INPUT_PREFIX.value+detailed_input+' '+xyz_filename+LinuxCommands.CREST_SUFIX.value)
    else:
        os.system(LinuxCommands.CREST_PREFIX.value+xyz_filename+LinuxCommands.CREST_SUFIX.value)

def run_xtb_calcuation(xyz_filename, detailed_input=None):
    if detailed_input:
        os.system(LinuxCommands.XTB_INPUT_PREFIX.value+detailed_input+' '+xyz_filename+LinuxCommands.XTB_SUFIX.value)
    else:
        os.system(LinuxCommands.XTB_PREFIX.value+xyz_filename+LinuxCommands.XTB_SUFIX.value)

def grep_py(expression, text_filename):
    with open(text_filename, 'r') as file_one:
        lines=[line for line in file_one if re.search(expression, line)]
    return lines

if __name__=='__main__':
    print('script loaded successfully')
