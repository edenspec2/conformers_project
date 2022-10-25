from enum import Enum

def is_number(n):
    try:
        float(n)
        return True
    except ValueError:
        return  False

class LinuxCommands(Enum):
    OBABEL='/gpfs0/gaus/users/itamarwa/transfer_to_itamar/build/bin/obabel '
    OBABEL_XYZ_SETTINGS_1=' -O '
    OBABEL_XYZ_SETTINGS_2=' --gen3d'
    XTB_INPUT_PREFIX='/gpfs0/gaus/projects/xtb-6.4.1/bin/xtb --input ' # xtb_di.inp 
    XTB_PREFIX='/gpfs0/gaus/projects/xtb-6.4.1/bin/xtb '
    XTB_SUFIX=' --ohess --dipole --pop'
    CREST_PREFIX='/gpfs0/gaus/projects/crest '
    CREST_SUFIX=' --xnam /gpfs0/gaus/projects/xtb-6.4.1/bin/xtb'
    COPY_COMMAND='cp ' # location_1 -> location_2
    MOVE_COMMAND='mv ' # location_1 -> location_2

if __name__=='__main__':
    print('script loaded successfully')
