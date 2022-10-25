import sys
path_to_add=r'C:\Users\\itaro\OneDrive\Documents\GitHub'
sys.path.insert(0, path_to_add)

##from Crystal_structure.tools.general_functions import *
from tools.general_constants import *
from tools.file_handlers import *
from tools.linux_operations import *

class XTBOutputFilenames(Enum):
    CHARGES='charges'
    PROPERTIES='properties.out'
    WBO='wbo'
    OPT_XYZ='xtbopt.xyz'
##                      G98='g98.out'
##                      HESSIAN='hessian'
##                      HESSIAN_OUT='hessian.out'
##                      VIB='vibspectrum'
##                      JSON='xtbout.json'
##                      OPTIMIZED='xtbtopo.mol'
##                      SCREEN='xtbscreen.xyz'
##                      MO='xtbmoinfo'
##                      ESP1='xtb_esp_profile.dat'
##                      ESP2='xtb_esp.dat'
##                      ESP3='xtb_esp.cosmo'
##                      MOLDEN='molden.input'
##                      COORD='lmocent.coord'
##                      COORD2='coordprot.0'
##                      HESS_XYZ='xtbhess.xyz'
##                      OPT_LOG='xtbopt.log'
##                      XTBRESTART='xtbrestart'
##                      VIB_MODES='vib_normal_modes'

class XTBFeaturesNames(Enum):
    ENERGY_LEVELS=('HOMO-2', 'HOMO-1', 'HOMO', 'LUMO', 'LUMO+1', 'LUMO+2')
    DISPERSION_POLARIZABILITY=('C6', 'C8', 'polarizability')
    DIPOLE_Q=('dipole_q_x', 'dipole_q_y', 'dipole_q_z', 'dipole_q_tot')
    DIPOLE=('dipole_x', 'dipole_y', 'dipole_z', 'dipole_tot')
    QUAD_Q=('quad_q_xx', 'quad_q_xy', 'quad_q_yy', 'quad_q_xz', 'quad_q_yz', 'quad_q_zz')
    QUAD_Q_DIP=('quad_q_vib_xx', 'quad_q_vib_xy', 'quad_q_vib_yy', 'quad_q_vib_xz', 'quad_q_vib_yz', 'quad_q_vib_zz')
    QUAD_FULL=('quad_full_xx', 'quad_full_xy', 'quad_full_yy', 'quad_full_xz', 'quad_full_yz', 'quad_full_zz')
    INERTIA_MOMENTS=('inertia_x', 'inertia_y', 'inertia_z')
    BONDS_STATS=('avg_distance', 'max_distance', 'min_distance')
    THERMO_TYPES=('VIB', 'ROT', 'INT', 'TR', 'TOT')
    THERMO_NAMES=('enthalpy', 'heat_capacity', 'entropy')
    ZERO_ENERGY=('zero_energy')

class XTBConstants(Enum): 
    MOLECULE_XTB_OUTPUT='xtb_output.zip'
    FEATURES_XTB_SUFIX='_features.pickle'
    FEATURES_CSV_SUFIX='_features.csv'

if __name__=='__main__':
    print('script loaded successfully')
