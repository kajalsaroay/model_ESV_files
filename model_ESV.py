#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 20 11:35:56 2022

@author: c1519699
"""

import numpy as np 
import matplotlib.pyplot as plt 
import numpy.matlib
import math
import vessels
import os
from scipy import signal
 
""" Same as Boas_replication_vessel_classes (not including a sink segment) but with the addition of ESV. """ 

# Define pressures (mmHg)
minimum_pressure = vessels.minimum_pressure #change in vessels.py 
PIC = vessels.PIC 
print('Min pressure is: ' + str(minimum_pressure))
print('PIC is: ' + str(PIC))

# Proportion of total resistance in the sink segment
proportion_R = 0.5 #final segment has 50% of total resistance in the vessel

# %% 
# 1) Load in data listing type of object within each level 

segments = []

# Text file holds information about branch pattern, branches in the network and the type of object in each level
#filename = 'vessel_text_files/whole_network_compliance/whole_network/Boas_network_3.txt' 
#filename = 'vessel_text_files/within_level_compliance/within_level/Boas_within_level_8.txt'
#filename = 'vessel_text_files/realistic_network/test/realistic_network_test_10.txt'
#filename = 'vessel_text_files/whole_network_compliance/whole_network/Boas_network_100000.txt'

filename = 'vessel_text_files/realistic_network/dynamic_simulations/realistic_network_capillary.txt'

#filename = 'vessel_text_files/realistic_network/dynamic_simulations/realistic_network_baseline_2_capillary_10x_segments.txt'

#filename = 'vessel_text_files/realistic_network/changing_segments/realistic_network_baseline_2_capillary_10x_segments.txt'

with open(filename) as file_object: 
    lines = file_object.readlines()
    branch_pattern = eval(lines[0])
    branches = eval(lines[1])
    for line in lines[2:]:
        segments.append(eval(line))
        
total_vessels = np.sum(branch_pattern)
cumulative_total_vessels = np.cumsum(branch_pattern) 
cumulative_total_vessels = np.insert(cumulative_total_vessels, 0, 0)

# Create list of levels, includes correct number of branch objects in each level 
levels = []
n = 0
for i in range(0, len(branches)):
    branch_in_level = np.array([])
    for j in range(0, branches[i]):
        branch_in_level = np.append(branch_in_level, segments[n])
        n += 1 
    levels.append(branch_in_level)


# Addition of sink segment
num_segments = len(levels)
num_segments_sink = 0
total_segments = num_segments + num_segments_sink

required_time_steps = len(levels) * levels[0][0].num_time_steps

# %% 
# 2) Create labels for every vessel in network, save in 'branch_labels' list 
# 3) Find path of flow for every vessel in network, save in matrix 'vessels_in_path' 

# A) diverging or converging
branch_labels_A = np.zeros([np.max(branch_pattern), len(levels)])

# Create diverging(0)/converging(1) array 
y = []

for i in range(0, len(levels)):
    if i < ((len(levels))/2):
        y.append(0)
    else: 
        y.append(1)

for i in range(0, len(levels)):
    for j in range(0, branch_pattern[i]):
        branch_labels_A[j, i] = y[i]
         
# B) branch pattern
branch_labels_B = np.zeros([np.max(branch_pattern), len(levels)])
for i in range(0, len(levels)):
    for j in range(0, branch_pattern[i]):
#    if i <= branch_pattern[i]:
        branch_labels_B[j, i] = branch_pattern[i]

# C) number within branch pattern 
branch_labels_C = np.zeros([np.max(branch_pattern), len(levels)])
for i in range(0, len(levels)):
    for j in range(0, branch_pattern[i]):
        branch_labels_C[j, i] = j
         
# Generate label for every branch 
branch_labels = []
for i in range(0, len(levels)):
    for j in range(0, branch_pattern[i]):
        branch_labels.append(np.array((branch_labels_A[j, i], branch_labels_B[j, i], branch_labels_C[j, i])))

# Store branch labels from each level in arrays 
branch_labels_levels = []
for n in range(0, len(branch_pattern)):
    branch_labels_levels.append(branch_labels[cumulative_total_vessels[n]:cumulative_total_vessels[n+1]])


# Obtain vessels in the path of blood flow for every vessel in network 
# Create binary matrix to store path of blood flow for each vessel 
# 1: vessel is in path of blood flow, 0: vessel is not in path, 999: no vessels in level 
vessels_in_path = np.zeros([len(branch_labels), np.max(branch_pattern), len(levels)])
for i in range(0, len(branch_labels)):
    for j in range(0, len(branch_pattern)):
        for k in range(0, max(branch_pattern)):
            if k < branch_pattern[j]:
                vessels_in_path[i, k, j] = 0
            else:
                vessels_in_path[i, k, j] = 999
           
for i in range(0, len(branch_labels)):
    current_vessel = branch_labels[i]
    
    # Set current vessel as 1 in binary matrix 
    m = current_vessel[2] # vessel number within level, corresponds to binary matrix row
    m = int(m)
    
    a = current_vessel[0] # diverging or converging
    a = int(a)
    b = current_vessel[1] # branch pattern value, corresponds to binary matrix column 
    b = int(b)
    if a == 0: # if diverging, obtain index for corresonding value in first half of branch pattern 
        n = branch_pattern.index(b) 
    else: 
        n = branch_pattern.index(b, int((len(branch_pattern)/2)))
 
    vessels_in_path[i, m, n] = 1
        
    # Loop through levels (l) and vessels (v) to find vessels in path across whole network 
    for l in range(0, len(levels)-1):
        for v in range(0, np.max(branch_pattern)):
            # find branch label for the next vessel in the path 
            if vessels_in_path[i, v, l] == 1: 
                
                if l < len(levels)/2:
                    a = 0
                else: 
                    a = 1
                b = branch_pattern[l]
                c = v
                
                current_vessel = np.array([a, b, c])
                
                # Find vessels in blood flow path in next level 
                if current_vessel[0] == 0 and current_vessel[1] != np.max(branch_pattern): # if diverging and not max number of vessels in level
                    a_next = 0
                    b_next = current_vessel[1] * 2
                    
                    c1_next = current_vessel[2] * 2
                    c2_next = (current_vessel[2]*2) + 1    
                    
                    next_vessel = []
                    next_vessel_1 = np.array([a_next, b_next, c1_next])
                    next_vessel_2 = np.array([a_next, b_next, c2_next])
                    next_vessel.append(next_vessel_1)
                    next_vessel.append(next_vessel_2)
                    
                    num_vessels = len(next_vessel)
                    for k in range(0, num_vessels):
                        m = next_vessel[k][2]
                        m = int(m)
                        
                        a = next_vessel[k][0]
                        a = int(a)
                        b = next_vessel[k][1]
                        b = int(b)
                        if a == 0:
                            n = branch_pattern.index(b)
                        else: 
                            n = branch_pattern.index(b, int((len(branch_pattern)/2)))
                        
                        vessels_in_path[i, m, n] = 1
                    
                elif current_vessel[0] == 1: #if converging
                    a_next = 1
                    b_next = (current_vessel[1]/2)
                    c_next = math.floor(current_vessel[2]/2)
                    
                    next_vessel = np.array([a_next, b_next, c_next])
                    
                    num_vessels = 1
                    for k in range(0, num_vessels):
                        m = next_vessel[2]
                        m = int(m)
                        
                        a = next_vessel[0]
                        a = int(a)
                        b = next_vessel[1]
                        b = int(b)
                        if a == 0:
                            n = branch_pattern.index(b)
                        else: 
                            n = branch_pattern.index(b, int((len(branch_pattern)/2)))
                        
                        vessels_in_path[i, m, n] = 1
                    
                else:  # if in level with max number of vessels 
                    a_next = 1
                    b_next = np.max(branch_pattern)
                    c_next = current_vessel[2]
                    
                    next_vessel = np.array([a_next, b_next, c_next])
                    
                    num_vessels = 1
                    for k in range(0, num_vessels):
                        m = next_vessel[2]
                        m = int(m)
                        
                        a = next_vessel[0]
                        a = int(a)
                        b = next_vessel[1]
                        b = int(b)
                        if a == 0:
                            n = branch_pattern.index(b)
                        else: 
                            n = branch_pattern.index(b, int((len(branch_pattern)/2)))
                                               
                        vessels_in_path[i, m, n] = 1
                 

# Previous vessel matrices  
# For every vessel, find the vessel(s) in the previous level (output pressure[i] = input_pressure[i+1]) 
previous_vessel_matrix = np.zeros([len(branch_labels), np.max(branch_pattern), len(levels)])
for i in range(0, len(branch_labels)):
    for j in range(0, len(branch_pattern)):
        for k in range(0, max(branch_pattern)):
            if k < branch_pattern[j]:
                previous_vessel_matrix[i, k, j] = 0
            else:
                previous_vessel_matrix[i, k, j] = 999
                

for i in range(1, len(branch_labels)):
    current_vessel = branch_labels[i]
    a = current_vessel[0]
    b = current_vessel[1]
    c = current_vessel[2]
    
    if a == 0:
        a_prev = 0
        b_prev = b/2
        c_prev = math.floor(c/2)
        
        previous_vessel = np.array([a_prev, b_prev, c_prev])
        num_vessels = 1
        for k in range(0, num_vessels):
            m = previous_vessel[2] #row number
            m = int(m)
                            
            a = previous_vessel[0]
            a = int(a)
            b = previous_vessel[1]
            b = int(b)
            if a == 0:
                n = branch_pattern.index(b) #column number
            else: 
                n = branch_pattern.index(b, int((len(branch_pattern)/2)))
    
            previous_vessel_matrix[i, m, n] = 1
        
    elif a == 1 and b == np.max(branch_pattern):
        a_prev = 0
        b_prev = b
        c_prev = c
        
        previous_vessel = np.array([a_prev, b_prev, c_prev])

        num_vessels = 1
        for k in range(0, num_vessels):
            m = previous_vessel[2] #row number
            m = int(m)
                            
            a = previous_vessel[0]
            a = int(a)
            b = previous_vessel[1]
            b = int(b)
            if a == 0:
                n = branch_pattern.index(b) #column number
            else: 
                n = branch_pattern.index(b, int((len(branch_pattern)/2)))
    
            previous_vessel_matrix[i, m, n] = 1      
            
    else: 
        a_prev = 1
        b_prev = b*2
        c1_prev = c*2
        c2_prev = (c*2) + 1 
        
        previous_vessel = []
        previous_vessel_1 = np.array([a_prev, b_prev, c1_prev])
        previous_vessel_2 = np.array([a_prev, b_prev, c2_prev])
        previous_vessel.append(previous_vessel_1)
        previous_vessel.append(previous_vessel_2)
        
        num_vessels = len(previous_vessel)
        for k in range(0, num_vessels):
            m = previous_vessel[k][2]
            m = int(m)
            
            a = previous_vessel[k][0]
            a = int(a)
            b = previous_vessel[k][1]
            b = int(b)
            if a == 0:
                n = branch_pattern.index(b)
            else: 
                n = branch_pattern.index(b, int((len(branch_pattern)/2)))
            
            previous_vessel_matrix[i, m, n] = 1


# Mappings 
def levels_to_matrix(level, branch, vessel):
    """ Return the position in binary matrix corresponding to the vessel being indexed with levels list notation """
    row = (2*branch) + vessel 
    column = level
    return [row, column]

def matrix_to_levels(row, column):
    """ Return the index notation for levels index for corresponding vessel in binary matrix """
    level = column #level in network
    branch = math.floor(row/2) #branch within level
    vessel = row % 2 #vessel within branch (0 or 1 only)
    return [level, branch, vessel]

#def ESV_to_levels(ESV_number, num_time_steps):
#    """ Returns level and position within branched network corresponding to ESV segment number. 
#    Assumes all objects in all levels have the same number of segments/time steps """
#    level = math.floor(ESV_number/num_time_steps)
#    position = ESV_number % num_time_steps 
#    return [level, position]
    
# Functions     
def parallel_resistance(resistance_0, resistance_1):
    """ Return total resistance for two resistors in parallel """ 
    total_resistance = (resistance_0 * resistance_1)/(resistance_0 + resistance_1)
    return total_resistance 


# %% INITIAL VALUES                       
# 4) Create matrices to store initial values for resistance, pressure and volume across network                       
                        
# For indexing purposes in initial value matrices 
#y = np.array([0, 0, 1, 0, 1, 0]) # branch within level 
#z = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]) # vessel within branch (in order of loop)

# Obtain 0 and 1 values representing the branch number within the level 
branches_in_level = []
for i in range(0, len(branch_pattern)):
    if branch_pattern[i] == 1: 
        branches_in_level.append(int(0))
    else: 
        branches_in_level.append(int(branch_pattern[i]/2))
branches_in_level = branches_in_level[1:-1]
        
y = np.array([], dtype = int) # branch within level
for i in range(0, len(branches_in_level)):
    j = branches_in_level[i]
    for k in range(0, j):
        y = np.append(y, k)

z = np.array([], dtype = int) # vessel within branch
for i in range(0, len(branches_in_level)):
    j = branches_in_level[i]
    for k in range(0, j):
        z = np.append(z, np.array([0, 1]))
        

# Mapping - vessel number in branch labels list to position in binary matrix      
vessel_number_matrix = np.zeros([np.max(branch_pattern), len(levels)])
for i in range(0, len(branch_pattern)):
    for j in range(0, max(branch_pattern)):
        if j < branch_pattern[i]:
            vessel_number_matrix[j, i] = 0
        else:
            vessel_number_matrix[j, i] = 999
    
number = np.array(range(0, len(branch_labels)))   
n = 0
for i in range(0, len(branch_pattern)):
    for j in range(0, max(branch_pattern)):    
        if vessel_number_matrix[j, i] == 0:
            vessel_number_matrix[j, i] = number[n] 
            n = n+1

# Remaining_resistance_vessels
# All vessels in levels following vessel of interest
# Use this in remaining resistance equation for comparison with steady-state model 
remaining_resistance_vessels = np.zeros([len(branch_labels), np.max(branch_pattern), len(levels)])
for i in range(0, len(branch_labels)):
    for j in range(0, len(branch_pattern)):
        for k in range(0, max(branch_pattern)):
            if k < branch_pattern[j]:
                remaining_resistance_vessels[i, k, j] = 0
            else:
                remaining_resistance_vessels[i, k, j] = 999
  
for i in range(0, len(branch_labels)):
    m = np.where(vessel_number_matrix == i)[0]
    n = np.where(vessel_number_matrix == i)[1]         
    for j in range(0, len(levels)):
        for k in range(0, max(branch_pattern)):
            # find row and column position of i 
            if j > n[0]:
                if remaining_resistance_vessels[i, k, j] == 0:
                    remaining_resistance_vessels[i, k, j] = 1


# Create initial value matrices          
# Initial resistance, pressure and volume values for every vessel in network 
resistances = np.zeros([np.max(branch_pattern), len(levels)])
pressures = np.zeros([np.max(branch_pattern), len(levels)]) #input pressures 
volumes = np.zeros([np.max(branch_pattern), len(levels)])
output_pressures = np.zeros([np.max(branch_pattern), len(levels)])
diameters = np.zeros([np.max(branch_pattern), len(levels)])
#viscosities = np.zeros([np.max(branch_pattern), len(levels)])
for i in range(0, len(levels)): 
    for j in range(0, max(branch_pattern)):
        if j < branch_pattern[i]:
            resistances[j, i] = 0
            pressures[j, i] = 0
            volumes[j, i] = 0
            output_pressures[j, i] = 0
            diameters[j, i] = 0
            #viscosities[j, i] = 0
        else:
            resistances[j, i] = 999
            pressures[j, i] = 999
            volumes[j, i] = 999
            output_pressures[j, i] = 999 
            diameters[j, i] = 999
            #viscosities[j, i] = 999

for l in range(0, len(levels)):
    for v in range(0, np.max(branch_pattern)):
        if resistances[v, l] == 0: 
            x, y, z = matrix_to_levels(v, l)
            resistances[v, l] = levels[x][y].resistance[z]
            pressures[v, l] = levels[x][y].input_pressure[z]
            volumes[v, l] = levels[x][y].volume[z]
            output_pressures[v, l] = levels[x][y].output_pressure[z]
            diameters[v, l] = levels[x][y].diameter[z]
            #viscosities[v, l] = levels[x][y].viscosity[z]

# Calculate a0() for every object 
for level in range(0, len(levels)):
    for branch in range(0, len(levels[level])):
        levels[level][branch].calculate_a0() 

# %% MAIN LOOP 
# 5) Simulate blood flow in network across time 

#pressure_wave = np.load("pressure_wave_2.npy")
        
##sine wave between 25 and 60mmHg 
#t_step = 0.01
#t_end = 10. 
#t = np.arange(0, t_end, t_step)
#pressure_wave= (17.5 * np.sin(t-45.5)) + 42.5

        
## Steady-state pressure input 
#max_pressure = 60. #maximum pressure
#pw_len = 1000 
#pressure_wave = np.repeat(max_pressure, pw_len)
#print('Max pressure is: ', str(max_pressure)) 
#print('Simulation type: steady state')

# Pulsatile pressure input
pressure_wave = np.loadtxt('example_dimac_ts.txt')
pressure_wave = pressure_wave[0:200]
#pressure_wave = pressure_wave[0:90]
#
pressure_smoothed = signal.savgol_filter(pressure_wave, 21, 2)
pressure_wave = pressure_smoothed
#
#if (1):
#    pressure_wave = pressure_wave[19:-1] #first beat 
#    
#    
# Scaling pressure wave between 80 and 120 mmHg 
min_val = 80.
max_val = 120. 

pw_sd = np.array([])
pw_scaled = np.array([])
for i in range(0, len(pressure_wave)):
    sd = (pressure_wave[i]-min(pressure_wave))/(max(pressure_wave)-min(pressure_wave))
    pw_sd = np.append(pw_sd, sd)
    
    scaled_val = (sd * (max_val-min_val)) + min_val
    pw_scaled = np.append(pw_scaled, scaled_val)

pressure_wave = pw_scaled

##pressure_wave = pressure_wave[15:]
#
## add in baseline pressure between repeats
#if (0):
#    #baseline_pressure = np.repeat(min_val, 50)
#    baseline_pressure = np.repeat(81., 50)
#    pressure_wave = np.insert(pressure_wave, 0, baseline_pressure)
#
# Make pressure wave longer - more beats
pressure_wave_test = np.array([])
inc_factor = 5 # factor to increase number of beats (3) by 
for i in range(0, inc_factor):
    pressure_wave_test = np.concatenate((pressure_wave_test, pressure_wave))

pressure_wave = pressure_wave_test

# Add in Gaussian noise
if (0):
    noise = np.random.normal(0, 1, len(pressure_wave))
    pressure_wave = pressure_wave + noise 
    
# Repeating pressure wave points 
if (0):
    num_repeats = 10 #double, triple or quadruple 
    repeats = np.repeat(num_repeats, len(pressure_wave))
    pressure_wave = np.repeat(pressure_wave, repeats)
    
if (0):
    #baseline_pressure = np.repeat(min_val, 50)
    baseline_pressure = np.repeat(80., 200)
    pressure_wave = np.insert(pressure_wave, 0, baseline_pressure)
    
    
#    
## Pressure wave with random baseline lengths in between 
#if (1):
#    baseline_lengths = np.array([50, 10, 35, 30, 20, 15])
#    pressure_wave = np.array([])
#    # pressure_wave = np.repeat(80., 500)
#    for b in range(0, len(baseline_lengths)):
#        # Pulsatile pressure input
#        pressure_wave_1 = np.loadtxt('example_dimac_ts.txt')
#        pressure_wave_1 = pressure_wave_1[0:80]
#        
#        pressure_smoothed = signal.savgol_filter(pressure_wave_1, 21, 2)
#        pressure_wave_1 = pressure_smoothed
#        
#        
#        pressure_wave_1 = pressure_wave_1[19:-1] #first beat 
#            
#            
#        # Scaling pressure wave between 80 and 120 mmHg 
#        min_val = 80.
#        max_val = 120. 
#        
#        pw_sd = np.array([])
#        pw_scaled = np.array([])
#        for i in range(0, len(pressure_wave_1)):
#            sd = (pressure_wave_1[i]-min(pressure_wave_1))/(max(pressure_wave_1)-min(pressure_wave_1))
#            pw_sd = np.append(pw_sd, sd)
#            
#            scaled_val = (sd * (max_val-min_val)) + min_val
#            pw_scaled = np.append(pw_scaled, scaled_val)
#        
#        pressure_wave_1 = pw_scaled
#        
#        
#        pressure_wave_1 = np.append(pressure_wave_1, np.repeat(80., baseline_lengths[b]))   
#        
#        pressure_wave = np.append(pressure_wave, pressure_wave_1)
#
#
#    # Repeating pressure wave points 
#    if (1):
#        num_repeats = 10 #double, triple or quadruple 
#        repeats = np.repeat(num_repeats, len(pressure_wave))
#        pressure_wave = np.repeat(pressure_wave, repeats)
#    
#    first_baseline = np.repeat(80., 500)
#    pressure_wave = np.insert(pressure_wave, 0, first_baseline)

print('Simulation type: dynamic')
print('Pressure wave max: ', max_val)
print('Pressure wave min: ', min_val)
###############################################################

time_per_segment = 1 

time = np.array([])
for n in range(0, len(pressure_wave)+1):
    time_step = time_per_segment*n
    time = np.append(time, time_step)

# Create matrices to store values of resistance, pressure, volume and flow in every vessel in network across time 
resistance_matrix = np.zeros([len(time), np.max(branch_pattern), len(levels)])
pressure_matrix = np.zeros([len(time), np.max(branch_pattern), len(levels)])
volume_matrix = np.zeros([len(time), np.max(branch_pattern), len(levels)])
flow_matrix = np.zeros([len(time), np.max(branch_pattern), len(levels)])
output_pressure_matrix =  np.zeros([len(time), np.max(branch_pattern), len(levels)])
diameter_matrix =  np.zeros([len(time), np.max(branch_pattern), len(levels)])
remaining_resistance_matrix = np.zeros([len(time), np.max(branch_pattern), len(levels)])
#viscosity_matrix = np.zeros([len(time), np.max(branch_pattern), len(levels)])
for h in range(0, len(time)):
    for i in range(0, len(branch_pattern)):
        for j in range(0, max(branch_pattern)):
            if j < branch_pattern[i]:
                resistance_matrix[h, j, i] = 0
                pressure_matrix[h, j, i] = 0
                volume_matrix[h, j, i] = 0
                flow_matrix[h, j, i] = 0
                output_pressure_matrix[h, j, i] = 0
                diameter_matrix[h, j, i] = 0
                remaining_resistance_matrix[h, j, i] = 0
                #viscosity_matrix[h, j, i] = 0
            else:
                resistance_matrix[h, j, i] = 999
                pressure_matrix[h, j, i] = 999
                volume_matrix[h, j, i] = 999
                flow_matrix[h, j, i] = 999
                output_pressure_matrix[h, j, i] = 999
                diameter_matrix[h, j, i] = 999
                remaining_resistance_matrix[h, j, i] = 999
                #viscosity_matrix[h, j, i] = 999
                

# Set initial values in matrices 
resistance_matrix[0, :, :] = resistances
pressure_matrix[0, :, :] = pressures
volume_matrix[0, :, :] = volumes 
output_pressure_matrix[0, :, :] = output_pressures
diameter_matrix[0, :, :] = diameters 
#viscosity_matrix[0, :, :] = viscosities

# Set up Equivalent Single Vessel
num_segments = np.array([])
for i in range(0, len(levels)):
    num_segments = np.append(num_segments, levels[i][0].num_time_steps)
total_ESV_segments = int(np.sum(num_segments))
print('ESV has', total_ESV_segments, 'segments')
ESV_pressures = np.repeat(minimum_pressure, total_ESV_segments)
ESV_pressures = np.append(ESV_pressures, minimum_pressure) #add minimum pressure for whole pressure distribution

ESV_cumulative = np.cumsum(num_segments)

ESV_index_ranges = []
# ESV segment values
for i in range(0, len(levels)):
    if i == 0:
        index_min = 0
    else:
        index_min = int(ESV_cumulative[i-1]) 
    index_max = int(ESV_cumulative[i]-1)
    ESV_index_ranges.append(np.array([index_min, index_max]))
    
def ESV_to_levels(segment, index_ranges):  
    for r in range(0, len(index_ranges)):
        index_min = index_ranges[r][0]
        index_max = index_ranges[r][1]
        if segment in range(index_min, index_max+1):
            level_index = r 
    return level_index

def ESV_to_segment_position(segment, num_segments):
    segment_position = segment % num_segments
    return segment_position 

# Create segment flow matrix 
segment_flows = np.zeros([max(branch_pattern), len(levels)], dtype = object)
for i in range(0, len(levels)): 
    for j in range(0, max(branch_pattern)):
        if j < branch_pattern[i]:
            segment_flows[j, i] = 0
        else:
            segment_flows[j, i] = 999
        
for l in range(0, len(levels)):
    for v in range(0, np.max(branch_pattern)):
        if segment_flows[v, l] == 0: 
            x, y, z = matrix_to_levels(v, l)
            segment_flows[v, l] = levels[x][y].resistances[z]

segment_flow_matrix = np.zeros([len(time), np.max(branch_pattern), len(levels)], dtype = object)
for h in range(0, len(time)):
    for i in range(0, len(branch_pattern)):
        for j in range(0, max(branch_pattern)):
            if j < branch_pattern[i]:
                segment_flow_matrix[h, j, i] = 0
            else:
                segment_flow_matrix[h, j, i] = 999

segment_flow_matrix[0, :, :] = segment_flows

# ESV segment flow matrix - holds values from distribution calculation for main network NOT ESV calculation
ESV_flow_matrix = np.zeros([len(time), max(branch_pattern), total_ESV_segments])

# Hold values from ESV calculations 
ESV_pressure_matrix = np.zeros([len(pressure_wave), total_ESV_segments+1])
ESV_resistance_matrix = np.zeros([len(pressure_wave), total_ESV_segments])
ESV_flows = np.zeros([len(pressure_wave), total_ESV_segments])

resistance_ratio_matrix = np.zeros([len(time), max(branch_pattern), total_ESV_segments])
flow_ratio_matrix = np.zeros([len(time), max(branch_pattern), len(levels)])

ESV_last_segment_index = np.array([], dtype=int)
for i in range(0, len(ESV_index_ranges)):
    ind = int(ESV_index_ranges[i][-1])
    ESV_last_segment_index = np.append(ESV_last_segment_index, ind)
    
ESV_first_segment_index = np.array([], dtype=int)
for i in range(0, len(ESV_index_ranges)):
    ind = int(ESV_index_ranges[i][0])
    ESV_first_segment_index = np.append(ESV_first_segment_index, ind)


# Start of main loop     
for t in range(1, len(time)):   
    print('Time =', t, 'seconds')
    
    # 1) Change pressures
    # Loop through levels and branches (objects)
    # Take pressure input into each object as output from previous vessel at previous time
    # Update pressure arrays for objects    
    for level in range(0, len(levels)):
        for branch in range(0, len(levels[level])):
            if level == 0:
                levels[level][branch].change_pressure(pressure_wave[t-1])
            else: 
                new_pressure = np.array([])
                for vessel in range(0, levels[level][branch].number_of_vessels):
                    x = level
                    y = branch
                    z = vessel

                    # Find correct input pressure for each branch
                    # This is output pressure from previous vessel at previous time point (t-1) 
                    
                    m, n = levels_to_matrix(x, y, z) #find position in binary matrix 
                    
                    index = vessel_number_matrix[m, n]
                    index = int(index)
                                        
                    position = np.where(previous_vessel_matrix[index, :, :] == 1)
                    row_position = position[0][0]
                    column_position = position[1][0]
                    new_pressure = np.append(new_pressure, output_pressure_matrix[t-1, row_position, column_position])
                
                levels[level][branch].change_pressure(new_pressure)
    
    # Update volume, resistance, diameter and (input) pressure matrices for current time point
    for l in range(0, len(levels)):
        for v in range(0, np.max(branch_pattern)):
            if resistance_matrix[t, v, l] == 0: 
                x, y, z = matrix_to_levels(v, l)
                
                resistance_matrix[t, v, l] = levels[x][y].resistance[z] #total resistance in vessel
                volume_matrix[t, v, l] = levels[x][y].volume[z] #total volume of vessel
                diameter_matrix[t, v, l] = levels[x][y].diameter[z]
                pressure_matrix[t, v, l] = levels[x][y].input_pressure[z]
                #viscosity_matrix[t, v, l] = levels[x][y].viscosity[z]
            
    # 2) ESV calculation 
    # i) Create array of resistances for every segment in network 
    total_resistance_values = np.array([])
    for level in range(0, len(levels)):
        level_resistances = [] #arrays of segment resistance values for all vessels 
                    
        for branch in range(0, len(levels[level])):
            for vessel in range(0, levels[level][branch].number_of_vessels): 
                segment_resistances = levels[level][branch].resistances[vessel] #array of resistances for all segments in a vessel 
                level_resistances.append(segment_resistances) #holds arrays of vessel segment resistances 

        # Loop through segments 
        for segment in range(0, levels[level][branch].num_time_steps):
            resistance_values = np.array([])
#            for vessel in range(0, levels[level][branch].number_of_vessels):
            for vessel in range(0, len(level_resistances)):
                resistance_values = np.append(resistance_values, level_resistances[vessel][segment])
        
                # ii) Find total resistance value for each segment of ESV 
                while len(resistance_values)>1:
                    # Parallel calculation
                    # loop through pairs of values in resistance_values array until one value remaining
                    n0 = 0
                    n1 = 1
                    values = []
                    for j in range(0, len(resistance_values)//2):
                        values.append(parallel_resistance(resistance_values[n0], resistance_values[n1]))
                        n0 =+ 1
                        n1 =+ 1
                    resistance_values = values   
                    
            total_resistance_values = np.append(total_resistance_values, resistance_values[0])
            
    # total_resistance_values holds resistance values for all segments of ESV 
    ESV_resistances = np.copy(total_resistance_values)
    ESV_resistance_matrix[t-1, :] = ESV_resistances
    
    # iii) Calculate pressure distribution across ESV 
    # For every segment, calculate pressure drop 
    ESV_pressures = ESV_pressures[0:-1]
    ESV_pressures = np.insert(ESV_pressures, 0, pressure_wave[t-1])
    
    ESV_pressure_matrix[t-1, :] = ESV_pressures
    
    scaling_factor = np.ones([len(ESV_pressures)-1])
    for n in range(0, len(ESV_resistances)):
        total_remaining_resistance = np.sum(ESV_resistances[n:])
        pressure_drop = ((ESV_pressures[n]-minimum_pressure)/total_remaining_resistance) * ESV_resistances[n]
        perc_pressure_drop = pressure_drop/ESV_pressures[n]
        scaling_factor[n] = 1 - perc_pressure_drop

    ESV_output_pressures = np.copy(ESV_pressures)
    for n in range(0, len(ESV_pressures)-1):
        ESV_output_pressures[n] = (ESV_pressures[n] * scaling_factor[n])
    
    ESV_pressure_gradient = np.array([])
    for i in range(0, len(ESV_pressures)-1):
#        pressure_grad = ESV_pressures[i] - ESV_pressures[i+1]
        pressure_grad = ESV_pressures[i] - ESV_output_pressures[i]
        ESV_pressure_gradient = np.append(ESV_pressure_gradient, pressure_grad)
        
    # iv) Calculate flow through each segment in the ESV 
    ESV_flow = np.array([])
    for s in range(0, total_ESV_segments):
        flow = ESV_pressure_gradient[s]/ESV_resistances[s]
        ESV_flow = np.append(ESV_flow, flow)
        
    ESV_flows[t-1, :] = ESV_flow
        
    ESV_pressures = ESV_output_pressures #update ESV pressures ready for next iteration (t=t+1)
    
    # v) Distribute the flow into vessels in main network
    # If segment falls within first half of network: 
        # Find remaining resistance for vessels in path following vessel containing current segment
        # Add remaining resistance in current vessel from current segment onwards
    # Else: 
        # Find flow coming into segments (obtain from previous vessel or previous segment in current vessel)
    
    # Loop through segments in ESV - distribute flow ESV segments to corresponding segments in branched network 
    for segment in range(0, total_ESV_segments): 
    #for segment in range(0, 70):
        flow = ESV_flow[segment]
        
        segment_level = ESV_to_levels(segment, ESV_index_ranges) #level corresponding to current ESV segment
        segment_position = ESV_to_segment_position(segment, levels[segment_level][0].num_time_steps) #position within vessel corresponding to current ESV segment
        
        # Find branch, vessels in path matrix 
        vessels_in_level = branch_pattern[segment_level] #number of vessels in current level - therefore number of segments to distribute flow across 
        
        level = segment_level 

####### Updated to include changes with calculating ratio of resistance / distribution of flow ########################################
        
        # i) If in a single vessel 
        if level == 0 or level == len(levels)-1:
            ESV_segment_flow_values = np.array([])
            segment_flow = flow
            
            ratio_resistance = np.array([1])
            while len(ratio_resistance) < max(branch_pattern):
                ratio_resistance = np.append(ratio_resistance, 999)
            resistance_ratio_matrix[t, :, segment] = ratio_resistance
            
            ESV_segment_flow_values = np.append(ESV_segment_flow_values, segment_flow)
        
            while len(ESV_segment_flow_values) < max(branch_pattern):
                ESV_segment_flow_values = np.append(ESV_segment_flow_values, 999)
            
            ESV_flow_matrix[t, :, segment] = ESV_segment_flow_values
            
            
            if segment == ESV_last_segment_index[level]:
                flow_ratios = np.array([1])
                
                while len(flow_ratios) < max(branch_pattern):
                    flow_ratios = np.append(flow_ratios, 999)
                    
                flow_ratio_matrix[t, :, level] = flow_ratios

                
#            # Update segment_flow_matrix
#            level_index = level
#            for vessel_index in range(0, branch_pattern[level_index]):
#                flow_array = np.array([])
#                for i in range(ESV_index_ranges[level_index][0], ESV_index_ranges[level_index][1]+1):
#                    flow_array = np.append(flow_array, ESV_flow_matrix[t, vessel_index, i])
#                    
#                segment_flow_matrix[t, vessel_index, level_index] = flow_array 
#
#            
#            # Update flow attribute 
#            l = level 
#            for b in range(0, len(levels[l])):
#                for v in range(0, levels[l][b].number_of_vessels):
#                    r, c = levels_to_matrix(l, b, v)
#                        
#                    levels[l][b].flows[v] = segment_flow_matrix[t, r, c]

            
        # ii) Calculate remaining resistances - if in diverging branch
        elif level < int(len(levels)/2):
            remaining_resistance_ratios = np.array([]) #store remaining resistance in path ratio for every segment 
            ESV_segment_flow_values_level = np.array([])
            for branch in range(0, len(levels[level])):
                remaining_resistances = np.array([])
                total_remaining_resistance = np.array([]) #store resistance values for vessels a and b in branch in this array
                vessel_indices = np.array([])
                for vessel in range(0, levels[level][branch].number_of_vessels):
                    # Find remaining resistance for each vessel in branch 
                    # Find correct vessels_in_path matrix for current vessel 
                    m, n = levels_to_matrix(level, branch, vessel)
                    
                    index = int(vessel_number_matrix[m, n])
                    vessel_indices = np.append(vessel_indices, index)
                    remaining_resistance_vessels_temp = np.copy(vessels_in_path[index, :, :])
                    remaining_resistance_vessels_temp[m, n] = 0 #set current vessel as 0 in path, not included in calculation 
                    resistance_t = np.copy(resistance_matrix[t, :, :]) #copy of resistance matrix at current time point
                    resistance_in_path = resistance_t * remaining_resistance_vessels_temp
                    for column in range(0, len(levels)):
                        for row in range(0, np.max(branch_pattern)):
                            if resistance_in_path[row, column] == 999*999:
                                resistance_in_path[row, column] = 0 #if there is no vessel in position, set as 0
                
        
                    # Calculate parallel resistance in each remaining level 
                    level_total_resistance = np.array([])
                    for l in range(0, len(levels)-level): #to not include last vessel 
                        if l > level:
                            resistance_values = resistance_in_path[:, l] #Pick column for current level
                            # Remove values = 0 
                            resistance_values = resistance_values[resistance_values != 0]
                            while len(resistance_values)>1:
                                # Parallel calculation
                                # loop through pairs of values in resistance_values array until one value remaining
                                n0 = 0
                                n1 = 1
                                values = []
                                for j in range(0, len(resistance_values)//2):
                                    values.append(parallel_resistance(resistance_values[n0], resistance_values[n1]))
                                    n0 =+ 1
                                    n1 =+ 1
                                resistance_values = values   
    #                            print(l)
    #                            print(resistance_values)
                                
                            level_total_resistance = np.append(level_total_resistance, resistance_values[0])
        
        
                    # Calculate remaining resistance in current vessel from current segment onwards
                    current_vessel_resistance = np.sum(levels[level][branch].resistances[vessel][segment_position:])
                    
                    total_remaining_resistance = current_vessel_resistance + np.sum(level_total_resistance)
                    
                    remaining_resistances = np.append(remaining_resistances, total_remaining_resistance) #holds values for resistance in path for each vessel   
                    
                
                ESV_segment_flow_values = np.array([])
                ratios = np.array([])
                for vessel_in_branch in range(0, len(remaining_resistances)): #per vessel in branch 
                    ratio_resistance = remaining_resistances[vessel_in_branch]/np.sum(remaining_resistances)
                    ratios = np.append(ratios, ratio_resistance)
                 
                    # Ratio of flow for previous vessels
                    # Find previous ratio for all levels before current level 
                    combined_prev_vessels_ratio = np.array([]) # holds the ratio of flow for previous vessels 

                    x = level
                    y = branch
                    z = vessel_in_branch
                    m, n = levels_to_matrix(x, y, z) #find position of vessel in binary matrix 
                    index = vessel_number_matrix[m, n]
                    index = int(index)
                    
                    current_vessel_index = index
                    for i in range(level, 0, -1):        
                        position = np.where(previous_vessel_matrix[current_vessel_index, :, :] == 1)
                        vessel_position= position[0][0]
                        level_position = position[1][0]
                    
                        # find ratio of flow into previous vessel 
                        ratio_flow = flow_ratio_matrix[t, vessel_position, level_position]
                        combined_prev_vessels_ratio = np.append(combined_prev_vessels_ratio, ratio_flow)
                        
                        current_vessel_index = int(vessel_number_matrix[vessel_position, level_position])


                    # Multiply values in combined_prev_vessels_ratio together - calculation per vessel
                    flow_ratios_multiplied = np.prod(combined_prev_vessels_ratio)
                    
                    segment_flow = (1-ratio_resistance) * flow_ratios_multiplied * flow #distribution of ESV segment flow
                        
                    ESV_segment_flow_values = np.append(ESV_segment_flow_values, segment_flow)
                    
                remaining_resistance_ratios = np.append(remaining_resistance_ratios, ratios) # ratios for corresponding segments in main network 
                ESV_segment_flow_values_level = np.append(ESV_segment_flow_values_level, ESV_segment_flow_values)
                
                       
            if segment == ESV_last_segment_index[level]:
                flow_ratios = 1-remaining_resistance_ratios

            while len(flow_ratios) < max(branch_pattern): 
                flow_ratios = np.append(flow_ratios, 999)
                
            flow_ratio_matrix[t, :, level] = flow_ratios

            # Add ratios to ratio_resistance_matrix 
            while len(remaining_resistance_ratios) < max(branch_pattern): 
                remaining_resistance_ratios = np.append(remaining_resistance_ratios, 999)
            
            resistance_ratio_matrix[t, :, segment] = remaining_resistance_ratios       
            
            while len(ESV_segment_flow_values_level) < max(branch_pattern):
                ESV_segment_flow_values_level = np.append(ESV_segment_flow_values_level, 999)
            
            ESV_flow_matrix[t, :, segment] = ESV_segment_flow_values_level 

                   
        # iii) Converging branch 
        else:
            
            ESV_segment_flow_values = np.array([])
            # No flow in current iteration 
            if flow == 0: # ESV segment flow = 0 
                for branch in range(0, len(levels[level])):
                    for vessel in range(0, levels[level][branch].number_of_vessels):
                        ESV_segment_flow_values = np.append(ESV_segment_flow_values, 0)
                        
                while len(ESV_segment_flow_values) < max(branch_pattern):
                    ESV_segment_flow_values = np.append(ESV_segment_flow_values, 999)
 

            else: 
                ESV_segment_flow_values_level = np.array([])
                for branch in range(0, len(levels[level])):
                    ESV_segment_flow_values = np.array([])
                    for vessel_in_branch in range(0, len(remaining_resistances)): #per vessel in branch 
                    
                        # Ratio of flow for previous vessels
                        # Find previous ratio for all levels before current level 
                        combined_prev_vessels_ratio = np.array([]) # holds the ratio of flow for previous vessels 
    
                        x = level
                        y = branch
                        z = vessel_in_branch
                        m, n = levels_to_matrix(x, y, z) #find position of vessel in binary matrix 
                        
                        # Position in corresponding level on diverging side (symmetrical level)
                        corresponding_level = len(levels) - level - 1
                        
                        index = vessel_number_matrix[m, corresponding_level]
                        index = int(index)
                        ratio_flow = flow_ratio_matrix[t, m, corresponding_level]
                        combined_prev_vessels_ratio = np.append(combined_prev_vessels_ratio, ratio_flow)
                        
                        current_vessel_index = index
                        for i in range(corresponding_level, 0, -1):     
                            position = np.where(previous_vessel_matrix[current_vessel_index, :, :] == 1)
                            vessel_position = position[0][0]
                            level_position = position[1][0]

                            # find ratio of flow from previous vessels 
                            ratio_flow = flow_ratio_matrix[t, vessel_position, level_position]
                            combined_prev_vessels_ratio = np.append(combined_prev_vessels_ratio, ratio_flow)
                            
                            current_vessel_index = int(vessel_number_matrix[vessel_position, level_position])
                   
                        # Multiply values in combined_prev_vessels_ratio together - calculation per vessel
                        flow_ratios_multiplied = np.prod(combined_prev_vessels_ratio)
                        
                        segment_flow =  flow_ratios_multiplied * flow #distribution of ESV segment flow
                            
                        ESV_segment_flow_values = np.append(ESV_segment_flow_values, segment_flow)
                        
                    ESV_segment_flow_values_level = np.append(ESV_segment_flow_values_level, ESV_segment_flow_values)


#                if segment == ESV_last_segment_index[level]:
#                    flow_ratios = 1-remaining_resistance_ratios
#    
#                while len(flow_ratios) < max(branch_pattern): 
#                    flow_ratios = np.append(flow_ratios, 999)
#                    
#                flow_ratio_matrix[t, :, level] = flow_ratios
#    
#                # Add ratios to ratio_resistance_matrix 
#                while len(remaining_resistance_ratios) < max(branch_pattern): 
#                    remaining_resistance_ratios = np.append(remaining_resistance_ratios, 999)
#                
#                resistance_ratio_matrix[t, :, segment] = remaining_resistance_ratios       
                
                while len(ESV_segment_flow_values_level) < max(branch_pattern):
                    ESV_segment_flow_values_level = np.append(ESV_segment_flow_values_level, 999)
                
                ESV_flow_matrix[t, :, segment] = ESV_segment_flow_values_level 


    # Update segment_flow_matrix for every level                     
    # 3) Save in object's flow array attribute 
    # Loop through ESV index ranges 
    for level_index in range(0, len(ESV_index_ranges)):
        for vessel_index in range(0, branch_pattern[level_index]):
            flow_array = np.array([])
            for i in range(ESV_index_ranges[level_index][0], ESV_index_ranges[level_index][1]+1):
                flow_array = np.append(flow_array, ESV_flow_matrix[t, vessel_index, i])
                
            segment_flow_matrix[t, vessel_index, level_index] = flow_array 
            
    # Save segment flow values as object attribute 
    for l in range(0, len(levels)):
        for b in range(0, len(levels[l])):
            for v in range(0, levels[l][b].number_of_vessels):
                r, c = levels_to_matrix(l, b, v)
                
                levels[l][b].flows[v] = segment_flow_matrix[t, r, c]

 
    # Update flow attribute value for every object
    for l in range(0, len(levels)):
        for b in range(0, len(levels[l])):               
            levels[l][b].update_flows() # update self.flow_values and self.flow_matrix for every object at every time step

    # Update flow matrix 
    for l in range(0, len(levels)):
        for v in range(0, np.max(branch_pattern)):
            if flow_matrix[t, v, l] == 0: 
                x, y, z = matrix_to_levels(v, l)
                
                flow_matrix[t, v, l] = levels[x][y].flow[z] #middle segment value of flow in vessel


    # 4) Calculate pressures 
    # Update output pressure matrix             
    for level in range(0, len(levels)):
        for branch in range(0, len(levels[level])):
            levels[level][branch].calculate_pressures()
            
            for vessel in range(0, levels[level][branch].number_of_vessels):
                position = levels_to_matrix(level, branch, vessel)
                row = position[0]
                column = position[1]
                output_pressure_matrix[t, row, column] = levels[level][branch].output_pressure[vessel]


#print('The pressure distribution is:', np.append(pressure_matrix[150, 0, :], minimum_pressure))


# %% 
# Save outputs 
                
if (1):
    # Last time point
    pressure_path_a = pressure_matrix[-1, 0, :]
    flow_path_a = flow_matrix[-1, 0, :]
    volume_path_a = volume_matrix[-1, 0, :]
    
    path_b_index = [i-1 for i in branch_pattern]
    pressure_path_b = np.array([])
    flow_path_b = np.array([])
    volume_path_b = np.array([])
    for n in range(0, len(path_b_index)):
        pressure_path_b = np.append(pressure_path_b, pressure_matrix[-1, path_b_index[n], n])
        flow_path_b = np.append(flow_path_b, flow_matrix[-1, path_b_index[n], n])
        volume_path_b = np.append(volume_path_b, volume_matrix[-1, path_b_index[n], n])

    # Flow for every time point in ICA, MCA and capillaries (and corresponding levels for path a)
    flow_ICA_ts_a = np.array([])
    flow_MCA_ts_a = np.array([])
    flow_cap_ts_a = np.array([])
    flow_level2_ts_a = np.array([])
    for n in range(0, len(time)):
        flow_ICA_ts_a = np.append(flow_ICA_ts_a, flow_matrix[n, 0, 0])
        flow_MCA_ts_a = np.append(flow_MCA_ts_a, flow_matrix[n, 0, 1])
        flow_cap_ts_a = np.append(flow_cap_ts_a, flow_matrix[n, 0, 6])
        flow_level2_ts_a = np.append(flow_level2_ts_a, flow_matrix[n, 0, 2])

    flow_ICA_ts_b = np.array([])
    flow_MCA_ts_b = np.array([])
    flow_cap_ts_b = np.array([])
    flow_level2_ts_b = np.array([])
    for n in range(0, len(time)):
        flow_ICA_ts_b = np.append(flow_ICA_ts_b, flow_matrix[n, path_b_index[0], 0])
        flow_MCA_ts_b = np.append(flow_MCA_ts_b, flow_matrix[n, path_b_index[1], 1])
        flow_cap_ts_b = np.append(flow_cap_ts_b, flow_matrix[n, path_b_index[6], 6])
        flow_level2_ts_b = np.append(flow_level2_ts_b, flow_matrix[n, path_b_index[2], 2])
        
    # save path b AND corresponding path a values 
    # (time points, path)
    # Column 0: path a, column 1: path b
    flow_ICA = np.zeros([len(time), 2])
    flow_MCA = np.zeros([len(time), 2])
    flow_cap = np.zeros([len(time), 2])
    flow_level2 = np.zeros([len(time), 2])
    
    flow_ICA[:, 0] = flow_ICA_ts_a
    flow_MCA[:, 0] = flow_MCA_ts_a
    flow_cap[:, 0] = flow_cap_ts_a
    flow_level2[:, 0] = flow_level2_ts_a
        
    flow_ICA[:, 1] = flow_ICA_ts_b
    flow_MCA[:, 1] = flow_MCA_ts_b
    flow_cap[:, 1] = flow_cap_ts_b
    flow_level2[:, 1] = flow_level2_ts_b

if (0):    
    # Volume for every time point for chosen levels 
    path_b_index = [i-1 for i in branch_pattern]
    
    volume_1_ts_a = np.array([])
    volume_6_ts_a = np.array([])
    volume_8_ts_a = np.array([])
    volume_11_ts_a = np.array([])
    for n in range(0, len(time)):
        volume_1_ts_a = np.append(volume_1_ts_a, volume_matrix[n, 0, 1])
        volume_6_ts_a = np.append(volume_6_ts_a, volume_matrix[n, 0, 6])
        volume_8_ts_a = np.append(volume_8_ts_a, volume_matrix[n, 0, 8])
        volume_11_ts_a = np.append(volume_11_ts_a, volume_matrix[n, 0, 11])
    
    volume_1_ts_b = np.array([])
    volume_6_ts_b = np.array([])
    volume_8_ts_b = np.array([])
    volume_11_ts_b = np.array([])
    for n in range(0, len(time)):
        volume_1_ts_b = np.append(volume_1_ts_b, volume_matrix[n, path_b_index[1], 1])
        volume_6_ts_b = np.append(volume_6_ts_b, volume_matrix[n, path_b_index[6], 6])
        volume_8_ts_b = np.append(volume_8_ts_b, volume_matrix[n, path_b_index[8], 8])
        volume_11_ts_b = np.append(volume_11_ts_b, volume_matrix[n, path_b_index[11], 11])

    # baseline volume for each level 
    volume_baseline_a = volume_matrix[0, 0, :]
    volume_baseline_b = np.array([])
    for n in range(0, len(path_b_index)):
        volume_baseline_b = np.append(volume_baseline_b, volume_matrix[0, path_b_index[n], n])
    
    
if (0):
    np.savetxt('output_values/ESV_outputs/UPDATED_ESV_model/realistic_network_ss/pressure_compliance_9_test_100_b.txt', pressure_path_b)
    np.savetxt('output_values/ESV_outputs/UPDATED_ESV_model/realistic_network_ss/pressure_compliance_9_test_100_a.txt', pressure_path_a)

