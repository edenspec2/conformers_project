# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 14:45:24 2022

@author: edens
"""

import os

if __name__ == '__main__':
    os.chdir(r'/gpfs0/gaus/users/edenspec/test')
    xyz_files=[filename for filename in os.listdir() if 'xyz' in filename]
   
    for xyz_file in xyz_files:
        os.system('/gpfs0/gaus/users/edenspec/working_obabel/build/bin/obabel {} -O {} --conformer --nconf 30--writeconformers'.format(*xyz_file,('obconformers_'+xyz_file)))
        os.system('/gpfs0/gaus/projects/crest --input {} -gfn2 -g h2o -T 4 --xnam /gpfs0/gaus/projects/xtb-6.4.1/bin/xtb'.format(xyz_file))
