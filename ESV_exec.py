# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 13:51:50 2022

@author: kajal
"""

import numpy as np 
import matplotlib.pyplot as plt 
import os
import shutil
import Boas_replication_ESV as ESV

# Loop through files 
directory_path = r'C:/Users/kajal/Documents/phd_year1/vessel_code/vessel_text_files/'

# Path to vessel text files required for simulation
simulation_directory = 'realistic_network/whole_network/'

#current_level = 13
#simulation_directory = 'level_compliance/level_' + str(current_level) + '/'

filename_updated = 'realistic_network' #last part of simulation_directory path - will be used in output filenames 

path = directory_path + simulation_directory 

# Relates to names of input text files 
#compliances = np.linspace(2,3,11)
#compliances = np.array([100000])
compliances = np.array([10])

# Original path for pressure and flow values 
path_outputs = r'C:/Users/kajal/Documents/phd_year1/vessel_code/output_values/ESV_outputs/'
path_output_pressure = path_outputs+'pressure/'
path_output_flow = path_outputs+'flow/' 
path_output_resistance = path_outputs+'resistance/'
path_output_volume = path_outputs+'volume/'

# Change this depending on simulation directory 
path_outputs_updated = path_outputs+simulation_directory #check this is where output files will be moved to 

# For iterations in loop 
n = 0

for current_file in sorted(os.listdir(path)): 
    file_number = compliances[n]
    if current_file.endswith(str(file_number)+'.txt'):
        print(current_file)
        n+=1
        
        # Copy contents of chosen text file into generic input text file 
        # i) Create input text file 
        input_file = directory_path + 'network_input.txt' #generic input file name 
        
        # ii) Copy contents from chosen network text file
        shutil.copy(path+current_file, input_file)
               
        print('Working on:', current_file)
        
        # iii) Run Boas_replication_ESV.py 
        ESV.main(file_number)
                    
        # iv) Rename pressure and flow arrays 
        # Move to correct directory 
        if (0):
            pressure_file= 'pressure_compliance_'+str(file_number)+'.txt'
            flow_file = 'flow_compliance_'+str(file_number)+'.txt' 
            resistance_file = 'resistance_compliance_'+str(file_number)+'.txt'
            volume_file = 'volume_compliance_'+str(file_number)+'.txt'
            shutil.move(path_output_pressure+pressure_file, path_outputs_updated+'pressure_1/pressure_'+filename_updated+'_'+str(file_number)+'.txt')
            shutil.move(path_output_flow+flow_file, path_outputs_updated+'flow_1/flow_'+filename_updated+'_'+str(file_number)+'.txt')
#            shutil.move(path_output_resistance+resistance_file, path_outputs_updated+'resistance/resistance_'+filename_updated+'_'+str(file_number)+'.txt')
#            shutil.move(path_output_volume+volume_file, path_outputs_updated+'volume/volume_'+filename_updated+'_'+str(file_number)+'.txt')
        
        if (1):
            pressure_file_a = 'pressure_compliance_'+str(file_number)+'_path_a.txt'
            pressure_file_b = 'pressure_compliance_'+str(file_number)+'_path_b.txt'
            flow_file_a = 'flow_compliance_'+str(file_number)+'_path_a.txt'
            flow_file_b = 'flow_compliance_'+str(file_number)+'_path_b.txt'
            volume_file_a = 'volume_compliance_'+str(file_number)+'_path_a.txt'
            volume_file_b = 'volume_compliance_'+str(file_number)+'_path_b.txt'
            resistance_file_a = 'resistance_compliance_'+str(file_number)+'_path_a.txt'
            resistance_file_b = 'resistance_compliance_'+str(file_number)+'_path_b.txt'
            shutil.move(path_output_pressure+pressure_file_a, path_outputs_updated+'pressure/pressure_'+filename_updated+'_'+str(file_number)+'_path_a.txt')
            shutil.move(path_output_pressure+pressure_file_b, path_outputs_updated+'pressure/pressure_'+filename_updated+'_'+str(file_number)+'_path_b.txt')
            shutil.move(path_output_flow+flow_file_a, path_outputs_updated+'flow/flow_'+filename_updated+'_'+str(file_number)+'_path_a.txt')
            shutil.move(path_output_flow+flow_file_b, path_outputs_updated+'flow/flow_'+filename_updated+'_'+str(file_number)+'_path_b.txt')
            shutil.move(path_output_volume+volume_file_a, path_outputs_updated+'volume/volume_'+filename_updated+'_'+str(file_number)+'_path_a.txt')
            shutil.move(path_output_volume+volume_file_b, path_outputs_updated+'volume/volume_'+filename_updated+'_'+str(file_number)+'_path_b.txt')
            shutil.move(path_output_resistance+resistance_file_a, path_outputs_updated+'resistance/resistance_'+filename_updated+'_'+str(file_number)+'_path_a.txt')
            shutil.move(path_output_resistance+resistance_file_b, path_outputs_updated+'resistance/resistance_'+filename_updated+'_'+str(file_number)+'_path_b.txt')

        
        # Set n for next iteration
        n+=1

        # v) Delete input file for next iteration 
        os.remove(directory_path+'network_input.txt')
    