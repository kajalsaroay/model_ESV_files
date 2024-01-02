# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 15:20:54 2022

@author: kajal
"""

import numpy as np 
import matplotlib.pyplot as plt 
import math 
import os

# Change whole network compliance text files
# Change path
path = r'C:/Users/kajal/Documents/phd_year1/vessel_code/vessel_text_files/realistic_network/test/'
string_code = 'mmm' #string to replace in file 

#level = 'level_11'
#string_code = 'v11'
#path = r'C:/Users/kajal/Documents/phd_year1/vessel_code/vessel_text_files/level_compliance/'+level+'/'

#path = '/home/c1519699/phd_year1/vessel_text_files/art_ven_compliance/' #use full path

start_val = 4
end_val = 9
num_values = 6

compliance_values = np.linspace(start_val, end_val, num_values)
print(compliance_values)

#Change filename for generic text file 
file = path+'realistic_network_test_generic.txt'
file_number = 13

#file = path+'realistic_network_generic_wnc.txt' 

for c in range(0, len(compliance_values)):
    # Create new input file - copy generic text file and save with new suffix (file_number)
    
    #new_file = open(path+'Boas_network_'+level+'_compliance_' + str(compliance_values[c]) + '.txt', 'w')
    
    new_file = open(path+'realistic_network_test_'+ str(int(compliance_values[c])) + '.txt', 'w') 
    
    with open(file) as file_object:
        lines = file_object.readlines()
        
        for line in range(0, len(lines)):
            # Copy first two lines from text file 
            if line <2:
                new_file.write(lines[line])
            
            else: 
                # find 'ccc' in each line, replace with compliance value
                new_line = lines[line].replace(string_code, str(int(compliance_values[c])))
                # Replace lines in text file with new lines 
                new_file.write(new_line)
                
            print(line)
            
    new_file.close()
            
    file_number = file_number + 1