#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 15:25:44 2023

@author: c1519699
"""

import numpy as np 
import matplotlib.pyplot as plt 
import numpy.matlib
import math
import vessels
import os
import shutil 
import Boas_replication_ESV_UPDATED as ESV


# Run all DIMAC input networks 

# Loop through DIMAC files 
path = '/home/c1519699/phd_year1/vessel_text_files/'
directory_path = '/home/c1519699/phd_year1/vessel_text_files/DIMAC_networks_varying_segments/'

subject_number = 2                                                                                                                                                                                                                                                

#text_files = np.arange(1, 10, 1)
text_files = np.array([1])
#text_files = np.arange(1, 10, 1)


filename_updated = 'DIMAC_network_'

path_output = '/home/c1519699/phd_year1/output_values/ESV_outputs/DIMAC_temp/'

path_output_updated = '/home/c1519699/phd_year1/output_values/ESV_outputs/DIMAC/varying_segments/subject_'+str(subject_number)+'/'

# For iterations 
n = 0

for current_file in sorted(os.listdir(directory_path)):    
    file_number = text_files[n]
    if current_file.endswith(str(file_number)+'.txt'):
                
        # Copy contents of chosen text file into generic input text file named in Boas_replication_ESV_UPDATED.py
        # (i) Create input text file 
        input_file = path + 'network_input.txt'
        
        # (ii) Copy contents from chosen DIMAC text file 
        shutil.copy(directory_path+current_file, input_file)
        
        print('Working on subject: ', subject_number)
        print('Working on file: ', current_file)
        
        # (iii) Run Boas_replication_ESV_UPDATED.py 
        ESV.main(file_number)
        
        # Move files to correct directory 
        ICA_flow = path_output+'ICA_flow_DIMAC_'+str(file_number)+'.txt'
        ICA_pressure = path_output+'ICA_pressure_DIMAC_'+str(file_number)+'.txt'

        MCA_flow = path_output+'MCA_flow_DIMAC_'+str(file_number)+'.txt'
        MCA_pressure = path_output+'MCA_pressure_DIMAC_'+str(file_number)+'.txt'

        level2_flow = path_output+'level2_flow_DIMAC_'+str(file_number)+'.txt'
        level2_pressure = path_output+'level2_pressure_DIMAC_'+str(file_number)+'.txt'

        cap_flow = path_output+'cap_flow_DIMAC_'+str(file_number)+'.txt'
        cap_pressure = path_output+'cap_pressure_DIMAC_'+str(file_number)+'.txt'
       
        shutil.move(ICA_flow, path_output_updated+'ICA_flow_'+filename_updated+str(file_number)+'.txt')
        shutil.move(ICA_pressure, path_output_updated+'ICA_pressure_'+filename_updated+str(file_number)+'.txt')
        
        shutil.move(MCA_flow, path_output_updated+'MCA_flow_'+filename_updated+str(file_number)+'.txt')        
        shutil.move(MCA_pressure, path_output_updated+'MCA_pressure_'+filename_updated+str(file_number)+'.txt')        
        
        shutil.move(level2_flow, path_output_updated+'level2_flow_'+filename_updated+str(file_number)+'.txt')
        shutil.move(level2_pressure, path_output_updated+'level2_pressure_'+filename_updated+str(file_number)+'.txt')
        
        shutil.move(cap_flow, path_output_updated+'cap_flow_'+filename_updated+str(file_number)+'.txt')
        shutil.move(cap_pressure, path_output_updated+'cap_pressure_'+filename_updated+str(file_number)+'.txt')        
        
        # Set n for next iteration 
        n+=1
        
        # Delete input file for next iteration
        os.remove(path+'network_input.txt')
        
        
# If looping through subjects then add another command in loop to pick correct input pressure wave for subject