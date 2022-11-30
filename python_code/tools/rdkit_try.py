import sys
path_to_add=r'C:\Users\\itaro\OneDrive\Documents\GitHub'
sys.path.insert(0, path_to_add)

from Crystal_structure.tools.rmsd_wrappers import *

coordinates_1=np.array([[0.00000, 0.00000, 0.11779], [0.00000, 0.75545, -0.47116], [0.00000, -0.75545, -0.47116]])
coordinates_2=np.array([[1.00000, 1.00000, 1.11779], [1.00000, 1.75545, 0.52884], [1.00000, 0.24455, 0.52884]])
rmsd_score=rmsd_get_RMSD_score(coordinates_1, coordinates_2)
print(rmsd_score)
##print(translated_coordinates_2)

