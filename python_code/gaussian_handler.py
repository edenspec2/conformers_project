# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 12:01:35 2022

@author: edens
"""
from openbabel import pybel
import re
import pandas as pd
import numpy as np
import os
from enum import Enum

# import pymatgen
# from pymatgen.io.gaussian  import GaussianOutput

ob_log_handler = pybel.ob.OBMessageHandler()
ob_log_handler.SetOutputLevel(0)

os.chdir(r'C:\Users\edens\Documents\GitHub\main_python\gauss_log')

class ReExpressions(Enum):

    FLOAT= r'-?\d\.\d*'
    FLOATS_ONLY= "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?"
    BONDS= r'R\(\d+.\d+'
    FREQUENCY= r'\-{19}'

class FileFlags(Enum):
    DIPOLE_START='Dipole moment'
    DIPOLE_END='Quadrupole moment'
    MAIN_CUTOFF= r'\-{69}'
    STANDARD_ORIENTATION_START='Standard orientation:'
    POL_START='iso'
    POL_END='xx'


class Names(Enum):
    
    DIPOLE_COLUMNS=['dip_x','dip_y','dip_z','total']
    STANDARD_ORIENTATION_COLUMNS=['atom','x','y','z']
    DF_LIST=['standard_orientation_df', 'dipole_df', 'pol_df', 'atype_df', 'bonds_df', 'info_df']


class GeneralConstants(Enum):    
    ATOMIC_NUMBERS ={
     '1':'H', '5':'B', '6':'C', '7':'N', '8':'O', '9':'F', '14':'Si',
              '15':'P', '16':'S', '17':'Cl', '35':'Br', '53':'I', '27':'Co', '28':'Ni'}
     

    



def search_phrase_in_text(text_lines, key_phrase) :    
    search_result=re.compile(key_phrase).search(text_lines)
    return search_result

def extract_lines_from_text(text_lines, re_expression):
    selected_lines=re.findall(re_expression, text_lines)
    # if strip:
    #     selected_lines=selected_lines.strip()
    return selected_lines

    
def process_gaussian_dipole(log_file_lines):
    dipole_start=search_phrase_in_text(log_file_lines, key_phrase=FileFlags.DIPOLE_START.value)
    dipole_end=search_phrase_in_text(log_file_lines, key_phrase=FileFlags.DIPOLE_END.value)
    selected_lines=extract_lines_from_text(log_file_lines[dipole_start.start():dipole_end.start()], re_expression=ReExpressions.FLOAT.value)
    dipole_df=pd.DataFrame(selected_lines, index=Names.DIPOLE_COLUMNS.value,dtype=float).T
    return dipole_df

def gauss_first_split(log_file_lines):
    first_split=re.split(FileFlags.STANDARD_ORIENTATION_START.value,log_file_lines)[-1]
    gauss_data=re.split(FileFlags.MAIN_CUTOFF.value,first_split)
    return gauss_data

def process_gaussian_standard_orientation_df(log_file_lines):
    standard_orientation=np.array(extract_lines_from_text(log_file_lines, ReExpressions.FLOATS_ONLY.value )).reshape(-1,6)
    standard_orientation=np.delete(standard_orientation,(0,2),1)
    standard_orientation_df=pd.DataFrame(standard_orientation,columns=Names.STANDARD_ORIENTATION_COLUMNS.value)
    standard_orientation_df['atom'].replace(GeneralConstants.ATOMIC_NUMBERS.value, inplace=True)
    return standard_orientation_df

def process_gaussian_pol(log_file_lines):
    pol_start=search_phrase_in_text(log_file_lines, key_phrase=FileFlags.POL_START.value)
    pol_end=search_phrase_in_text(log_file_lines, key_phrase=FileFlags.POL_END.value)
    pol=extract_lines_from_text(log_file_lines[pol_start.start():pol_end.start()], re_expression=ReExpressions.FLOAT.value)
    pol_df=pd.DataFrame([float(pol[0])*1000,float(pol[3])*1000],index=['iso','aniso'],dtype=float)
    return pol_df

def process_gaussian_atype(standard_orientation_df):
    array_string = np.array2string(np.array(standard_orientation_df))
    array_string=(re.sub('[\[\]\']',' ',array_string))
    array_string=('{}\n{}\n'.format(standard_orientation_df.shape[0],''))+array_string
    pbmol=pybel.readstring('xyz',array_string)
    mol2_string=pbmol.write('mol2')
    t=np.array(re.split('\s+',(re.split('@',mol2_string)[2].split('ATOM')[1])))
    mol2_array=np.delete(t, (0,-1),0).reshape(-1,9)
    atype_df=pd.DataFrame(mol2_array[:,5])
    return atype_df

def process_gaussian_bonds(log_file_lines):
    bonds=extract_lines_from_text(log_file_lines, re_expression=ReExpressions.BONDS.value)
    bonds_df=pd.DataFrame([re.sub(r'R\(','',line).split(',') for line in bonds]).astype(int)
    return bonds_df

def process_gaussian_frequency_string(log_file_lines):
    frequency_string=re.split('Frequencies', log_file_lines)
    frequency_string.pop(0) #remove the first empty string
    frequency_string[-1]=re.split(ReExpressions.FREQUENCY.value, frequency_string[-1])[0]
    return frequency_string

def process_gaussian_vibs(frequency_string):
    vibs_list=[]
    for data in frequency_string:
        match=extract_lines_from_text(data, re_expression=ReExpressions.FLOATS_ONLY.value)
        del match[0:12] 
        match=np.array(match)
        try:
            vibs_list.append(match.reshape(-1,11))
        except ValueError:
            match=np.delete(match,[-1,-2,-3],0)
            vibs_list.append(np.array(match).reshape(-1,11))
    vibs=np.vstack(vibs_list).astype(float)
    vibs=vibs[vibs[:,0].argsort()]
    ordered_vibs_list=[vibs[i:i + len(vibs_list)] for i in range(0, len(vibs), len(vibs_list))]
    return vib_array_list_to_dict(ordered_vibs_list)

def process_gaussian_info(frequency_string):
    info_df=pd.DataFrame()
    ir,frequency=[],[]
    for data in frequency_string:
        match=extract_lines_from_text(data, re_expression=ReExpressions.FLOATS_ONLY.value)
        frequency.append((match[0:3]))
        ir.append(match[6:9])
    info_df['Frequency']=[float(item) for sublist in frequency for item in sublist]
    info_df['IR']=[float(item) for sublist in ir for item in sublist]    
    return info_df

def vib_array_list_to_dict(array_list, key_prefix='vibration_atom_'):
    my_dict={}
    for array in array_list:
        new_array=np.delete(array,[0,1],axis=1).reshape(-1,3)
        my_dict[key_prefix+'{}'.format(int(array[0,0]))]=new_array
    return my_dict

def df_list_to_dict(df_list):
    my_dict={}
    for name,df in zip(Names.DF_LIST.value,df_list):
        my_dict[name]=df
    return my_dict

def export_list_to_csv(filename, list_to_export):
    #make directory and change path to it
    path=os.path.abspath(filename)
    os.mkdir(path)
    os.chdir(path)
    #export list to csv
    for file in list_to_export:
        if type(file) == dict:
            for value in file:
                # print(type(file[value]))
                if type(file[value])==np.ndarray:
                    np.savetxt((str(value)+'.csv'), file[value],delimiter=",")  
                elif type(file[value])==pd.DataFrame:
                    file[value].to_csv(str(value)+'.csv')
    return


def gauss_file_handler(gauss_filename, export=False):
    with open(os.path.abspath(gauss_filename)) as f: #open files take memory - after loading file and the lines - close the file
        log_file_lines=f.read()
        #close file
        f.close()    
    dipole_df=process_gaussian_dipole(log_file_lines) 
    pol_df=process_gaussian_pol(log_file_lines)
    gauss_data=gauss_first_split(log_file_lines)
    standard_orientation_df=process_gaussian_standard_orientation_df(gauss_data[2])
    atype_df=process_gaussian_atype(standard_orientation_df)
    bonds_df=process_gaussian_bonds(gauss_data[4])
    frequency_string=process_gaussian_frequency_string(gauss_data[3])
    vibs_dict=process_gaussian_vibs(frequency_string)
    info_df=process_gaussian_info(frequency_string)
    filename=gauss_filename.split('.')[0]
    df_list=[standard_orientation_df,dipole_df,pol_df,atype_df,bonds_df,info_df]
    df_dict=df_list_to_dict(df_list)
    dict_list=[df_dict,vibs_dict]
    if export:
        export_list_to_csv(filename,dict_list)

    return dict_list

        
        
if __name__ == '__main__':
  
    # log_object= GaussianOutput(r'C:\Users\edens\Documents\GitHub\main_python\gauss_log\AB_para_methoxy.log')
    parameters=gauss_file_handler(r'AB_para_methoxy.log',export=False)

        
