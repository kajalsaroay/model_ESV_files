#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 11:29:11 2023

@author: c1519699
"""

import numpy as np 
import matplotlib.pyplot as plt 
import numpy.matlib
import math
import vessels
import os
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy import signal

# Calculate lengths and diameters of vessels in every level of realistic network 
# Values taken from Boas (2008) and Piechnik (2008)
# 3 compartments: 1) arteries, 2) microvessels, 3) veins 

branch_pattern = np.array([1, 2, 4, 8, 16, 32, 64, 64, 32, 16, 8, 4, 2, 1]) #number of vessels in each level 

# Define functions 
def calc_flow(velocity, diameter):
    """ Calculate flow from blood velocity and vessel diameter """
    f = (np.pi*velocity*((diameter)**2))/4 
    return f 

def calc_resistance(p1, p2, flow):
    """ Calculate resistance using change in pressure and flow """
    R = (p1-p2)/flow
    return R

def calc_length(resistance, viscosity, diameter):
    """ Calculate length using Poiseuille equation - known resistance, viscosity and diameter """
    l = (((resistance*133.322)*np.pi*((diameter)**4))/(128*(viscosity*1e-3))) #Convert R to Pas/mm^3 and viscosity to Pas  
    return l

def calc_diameter(resistance, viscosity, length):
    """ Calculate diameter using Poiseuille equation - known resistance, viscosity and length """
    d = ((128*(viscosity*1e-3)*length)/(np.pi*(resistance*133.322)))**(1/4) 
    return d

#######################################################################
# To estimate viscosity from diameter (interpolation) - use for microvessels only 
# Diameter and viscosity values from VAN model (Boas et al., 2008)
initial_diameter = np.array([30.5, 24.4, 19.5, 15.6, 12.5, 10., 8, 12., 15., 18.7, 23.4, 29.3, 36.6])
initial_viscosity = np.array([2.49, 2.34, 2.25, 2.20, 2.16, 2.12, 2.10, 2.15, 2.18, 2.22, 2.32, 2.51, 2.70]) #Convert to different units (x10**-3)
estimated_viscosity = interp1d(initial_diameter, initial_viscosity, fill_value="extrapolate") #in cP

if (0):
    plt.figure()
    plt.plot(initial_diameter, initial_viscosity, 'x')
    plt.xlabel('Diameter (\u03BCm)')
    plt.ylabel('Viscosity (cP)')
    plt.title('Viscosity against diameter for microvessels')
    
########################################################################


########################################################################
# To estimate input pressure for diameters in microvessel compartment (obtain desired pressure distribution)
# Diameter and pressure values from VAN model (Boas et al., 2008)
#Boas_pressure_dist = np.array([60., 58.32811981, 56.41019377, 54.14977327, 51.45178739, 48.23886899, 44.38947205, 32.75281805, 30.87014588,
#29.3063378 , 27.98774958, 26.8637201 , 25.87428611])
#Boas_diameters = np.array([-30.5, -24.4, -19.5, -15.6, -12.5, -10., 8., 12., 15., 18.7, 23.4, 29.3, 36.6]) * 1e-3


#Boas_pressure_dist = np.array([60., 58.32811981, 56.41019377, 54.14977327, 51.45178739, 48.23886899, 44.38947205, 38.57115663, 32.75281805, 30.87014588,
#29.3063378 , 27.98774958, 26.8637201 , 25.87428611])
#Boas_diameters = np.array([-30.5, -24.4, -19.5, -15.6, -12.5, -10., -8., 8., 12., 15., 18.7, 23.4, 29.3, 36.6]) * 1e-3

#plt.figure()
#plt.plot(Boas_diameters, Boas_pressure_dist)
#spl = UnivariateSpline(Boas_diameters, Boas_pressure_dist)

###########################################################################

# Create storage matrices 
realistic_diameter_matrix = np.zeros([np.max(branch_pattern), len(branch_pattern)])
realistic_length_matrix = np.zeros([np.max(branch_pattern), len(branch_pattern)])
realistic_flow_matrix = np.zeros([np.max(branch_pattern), len(branch_pattern)])
realistic_resistance_matrix = np.zeros([np.max(branch_pattern), len(branch_pattern)])
realistic_viscosity_matrix = np.zeros([np.max(branch_pattern), len(branch_pattern)])

for i in range(0, len(branch_pattern)): 
    for j in range(0, max(branch_pattern)):
        if j < branch_pattern[i]:
            realistic_diameter_matrix[j, i] = 0
            realistic_length_matrix[j, i] = 0
            realistic_flow_matrix[j, i] = 0
            realistic_resistance_matrix[j, i] = 0
            realistic_viscosity_matrix[j, i] = 0
            
        else:
            realistic_diameter_matrix[j, i] = 999
            realistic_length_matrix[j, i] = 999
            realistic_flow_matrix[j, i] = 999
            realistic_resistance_matrix[j, i] = 999
            realistic_viscosity_matrix[j, i] = 999
            
# %% 
# Compartment 1: Arteries 
# levels 0-2 

num_vessels_art = np.array([1, 2, 4]) #number of vessels in each level

# Choose ICA and MCA diameters and velocities
d_ICA = 4 #mm
v_ICA = 1250 #mm/s

d_MCA = 1.92 #mm
v_MCA = 270 #mm/s

viscosity_art = np.array([2.5, 2.5, 2.5]) #Assume same viscosity value for arterial levels

start_pressure_art = 120.
end_pressure_art = 60.
pressure_dist_art = np.linspace(start_pressure_art, end_pressure_art, 4) #Assume linear drop in pressure 

# Calculate flow in levels 0 and 1
f_ICA = calc_flow(v_ICA, d_ICA)
f_MCA = calc_flow(v_MCA, d_MCA)
f_else = f_ICA - f_MCA 

# Level 0 (ICA)
pressure_drop_level_0 = pressure_dist_art[0]-pressure_dist_art[1]
resistance_ICA = pressure_drop_level_0/f_ICA 

length_level_0 = calc_length(resistance_ICA, viscosity_art[0], d_ICA)

# Level 1 (MCA and everything else)
pressure_drop_level_1 = pressure_dist_art[1]-pressure_dist_art[2]
resistance_MCA = pressure_drop_level_1/f_MCA

length_level_1 = calc_length(resistance_MCA, viscosity_art[1], d_MCA)

resistance_else = pressure_drop_level_1/f_else
d_else = calc_diameter(resistance_else, viscosity_art[1], length_level_1)

# Level 2 (distribute flow from previous level)
f_branch_0 = np.array([f_else/2, f_else/2]) #flow in each vessel in branch 
f_branch_1 = np.array([f_MCA/2, f_MCA/2])

pressure_drop_level_2 = pressure_dist_art[2]-pressure_dist_art[3]
resistance_branch_0 = np.array([pressure_drop_level_2/f_branch_0[0], pressure_drop_level_2/f_branch_0[1]])
resistance_branch_1 = np.array([pressure_drop_level_2/f_branch_1[0], pressure_drop_level_2/f_branch_1[1]])

# Calculate diameters 
length_level_2 = length_level_1 #set lengths in level 2 to be the same as level 1

d_branch_0 = calc_diameter(resistance_branch_0, viscosity_art[2], length_level_2)
d_branch_1 = calc_diameter(resistance_branch_1, viscosity_art[2], length_level_2)

# Store values in matrices 
realistic_diameter_matrix[0, 0] = d_ICA
realistic_diameter_matrix[0, 1] = d_else
realistic_diameter_matrix[1, 1] = d_MCA
realistic_diameter_matrix[0, 2] = d_branch_0[0]
realistic_diameter_matrix[1, 2] = d_branch_0[1]
realistic_diameter_matrix[2, 2] = d_branch_1[0]
realistic_diameter_matrix[3, 2] = d_branch_1[1]

realistic_length_matrix[0, 0] = length_level_0
realistic_length_matrix[0, 1] = length_level_1
realistic_length_matrix[1, 1] = length_level_1
realistic_length_matrix[0, 2] = length_level_2
realistic_length_matrix[1, 2] = length_level_2
realistic_length_matrix[2, 2] = length_level_2
realistic_length_matrix[3, 2] = length_level_2

realistic_flow_matrix[0, 0] = f_ICA
realistic_flow_matrix[0, 1] = f_else
realistic_flow_matrix[1, 1] = f_MCA 
realistic_flow_matrix[0, 2] = f_branch_0[0]
realistic_flow_matrix[1, 2] = f_branch_0[1]
realistic_flow_matrix[2, 2] = f_branch_1[0]
realistic_flow_matrix[3, 2] = f_branch_1[1]

realistic_resistance_matrix[0, 0] = resistance_ICA
realistic_resistance_matrix[0, 1] = resistance_else
realistic_resistance_matrix[1, 1] = resistance_MCA
realistic_resistance_matrix[0, 2] = resistance_branch_0[0]
realistic_resistance_matrix[1, 2] = resistance_branch_0[1]
realistic_resistance_matrix[2, 2] = resistance_branch_1[0]
realistic_resistance_matrix[3, 2] = resistance_branch_1[1]

realistic_viscosity_matrix[0, 0] = viscosity_art[0]
realistic_viscosity_matrix[0, 1] = viscosity_art[1]
realistic_viscosity_matrix[1, 1] = viscosity_art[1]
realistic_viscosity_matrix[0, 2] = viscosity_art[2]
realistic_viscosity_matrix[1, 2] = viscosity_art[2]
realistic_viscosity_matrix[2, 2] = viscosity_art[2]
realistic_viscosity_matrix[3, 2] = viscosity_art[2]

##################################################################
# Calculate compliance for arteries 
P_IC = vessels.PIC
P_0 = pressure_dist_art[-1]

P_i_ICA = pressure_dist_art[0]
P_i_MCA = 100.6 # Change for input value from ESV only using ICA compliance 

def calc_arterial_compliance(ICP, P0, Pi, AC):
    """ Convert arterial compliance (AC) value (defined as % change in volume per mmHg pressure change) to Boas equivalent compliance value. 
        ICP is intracranial pressure, P0 is baseline pressure and PI is input pressure. """
    x = ((Pi-ICP)/(P0-ICP))
    base_AC = 1+ (AC/100)
    beta = math.log(x, base_AC)
    return beta

beta_ICA = np.round(calc_arterial_compliance(P_IC, P_0, P_i_ICA, 0.5), 1)
beta_MCA = np.round(calc_arterial_compliance(P_IC, P_0, P_i_MCA, 0.5), 1)

# %%
# Compartment 2: Microvessels 
# Levels 3-10 
# Distribution of blood flow in level 1 (MCA and everything else) means that the top half of the network will differ to the bottom half
# To keep vessel lengths the same within a level, calculate the diameters of vessels in top half of network 

num_vessels_micro = np.array([8, 16, 32, 64, 64, 32, 16, 8])

diameter_micro = np.array([30., 20., 10., 5.6, 5.6, 15., 30., 45.]) * 1e-3 #mm 
viscosity_micro = estimated_viscosity(diameter_micro*1000) #cP #estimated from Boas values #convert diameter to micrometres

#pressure_dist_micro = np.array([60., 56.4, 51.5, 44.4, 32.8, 29.3, 26.9, 25.9, 25.]) #input pressures, subset of Boas VAN input pressures 

#pressure_set = np.array([60., 56.4, 51.5, 44.4, 38.6, 32.8, 29.3, 26.9, 25.]) #ensures pressure drop across capillaries is the same as Boas network 
#pressure_dist_micro = pressure_set 

## Use UnivariateSpline to obtain input pressures for microvessels - finding middle value between levels 6 and 8
#if (0):
#    diameter_micro_spl = np.array([-30., -20., -10., 5.6, 15., 30., 36.6]) * 1e-3 #mm 
#    pressure_dist_micro = spl(diameter_micro_spl)
#    plt.figure()
#    plt.plot(pressure_dist_micro)
#    plt.xlabel('Level')
#    x_label_position = np.arange(len(diameter_micro_spl))
#    x_label = np.array([3, 4, 5, 6, 8, 9, 10])
#    plt.xticks(x_label_position, x_label)
#    plt.ylabel('Pressure (mmHg)')

## Create input pressure for second capillary level (level 7)
#p = (pressure_dist_micro[6-3]+pressure_dist_micro[8-3])/2
#pressure_dist_micro = np.insert(pressure_dist_micro, 4, p)
#
#if (0):
#    plt.figure()
#    plt.plot(pressure_dist_micro)
#    plt.xlabel('Level')
#    x_label_position = np.arange(len(pressure_dist_micro))
#    x_label = np.array([3, 4, 5, 6, 7, 8, 9, 10])
#    plt.xticks(x_label_position, x_label)
#    plt.ylabel('Pressure (mmHg)')
#    plt.title('Estimated microvessel input pressures')
#

## Use UnivariateSpline to obtain input pressures for microvessels - alternative method
#diameter_micro_spl = np.array([-30., -20., -10., -5.6, 5.6, 15., 30., 36.6]) * 1e-3 #mm 
#pressure_dist_micro = spl(diameter_micro_spl)
#plt.figure()
#plt.plot(pressure_dist_micro)
#plt.xlabel('Level')
#x_label_position = np.arange(len(diameter_micro_spl))
#x_label = np.array([3, 4, 5, 6, 7, 8, 9, 10])
#plt.xticks(x_label_position, x_label)
#plt.ylabel('Pressure (mmHg)')
#plt.title('Estimated microvessel input pressures')
#
#
## Split UnivariateSpline function into arterial and venous sides 
#diameter_micro_spl_art = np.array([30., 20., 10., 5.6]) * 1e-3
#pressure_dist_micro_art = spl_art(diameter_micro_spl_art)
#
#diameter_micro_spl_ven = np.array([5.6, 15., 30., 36.6]) * 1e-3
#pressure_dist_micro_ven = spl_ven(diameter_micro_spl_ven)
#
#pressure_dist_micro = np.concatenate((pressure_dist_micro_art, pressure_dist_micro_ven))
#
#plt.figure()
#plt.plot(pressure_dist_micro)
#plt.xlabel('Level')
#x_label_position = np.arange(len(diameter_micro_spl))
#x_label = np.array([3, 4, 5, 6, 7, 8, 9, 10])
#plt.xticks(x_label_position, x_label)
#plt.ylabel('Pressure (mmHg)')
#plt.title('Estimated microvessel input pressures')


# Kevin's method
pressure_dist_micro = np.array([59.8771, 56.6421, 48.2400, 43.5170, 39.4430, 30.8700, 26.7500, 24.8675])

pressure_dist_micro = np.append(pressure_dist_micro, 24.)

# Diverging microvessel flow
for i in range(3, 7):
    previous_flow = realistic_flow_matrix[0:(num_vessels_micro[i-3]//2), i-1]
    flow_values = np.array([])
    for j in range(0, len(previous_flow)):
        current_flow = np.repeat(previous_flow[j]/2, 2)
        flow_values = np.append(flow_values, current_flow)
    realistic_flow_matrix[0:num_vessels_micro[i-3], i] = flow_values
    
# Converging microvessel flow
n = 1
for i in range(7, 11): 
    symmetrical_flow = realistic_flow_matrix[:, i-n]
    realistic_flow_matrix[:, i] = symmetrical_flow
    n+=2 

# Store diameter values in matrix
for i in range(3, len(diameter_micro)+3):
    d = np.repeat(diameter_micro[i-3], num_vessels_micro[i-3])
    realistic_diameter_matrix[0:num_vessels_micro[i-3], i] = d
    
# Store viscosity values in matrix 
for i in range(3, len(viscosity_micro)+3):
    visc = np.repeat(np.round(viscosity_micro[i-3],2), num_vessels_micro[i-3])
    realistic_viscosity_matrix[0:num_vessels_micro[i-3], i] = visc

# Calculate resistance 
for i in range(3, len(diameter_micro)+3):
    flow_micro = realistic_flow_matrix[0:num_vessels_micro[i-3], i]
    resistance_micro = np.array([])
    for j in range(0, len(flow_micro)):
        r = calc_resistance(pressure_dist_micro[i-3], pressure_dist_micro[(i-3)+1], flow_micro[j])
        resistance_micro = np.append(resistance_micro, r)
    realistic_resistance_matrix[0:num_vessels_micro[i-3], i] = resistance_micro
    
# Calculate length from resistance
for i in range(3, len(diameter_micro)+3):
    r = realistic_resistance_matrix[0:num_vessels_micro[i-3], i]
    visc = np.repeat(viscosity_micro[i-3], num_vessels_micro[i-3])
    d = np.repeat(diameter_micro[i-3], num_vessels_micro[i-3])
    
    length_micro = np.array([])
    for j in range(0, len(r)):
        l = calc_length(r[j], visc[j], d[j])
        length_micro = np.append(length_micro, l)
    realistic_length_matrix[0:num_vessels_micro[i-3], i] = length_micro
    
for i in range(3, len(diameter_micro)+3):
    realistic_length_matrix[0:num_vessels_micro[i-3]//2, i] = realistic_length_matrix[num_vessels_micro[i-3]//2:num_vessels_micro[i-3], i]


############################################################################
# Calculate diameters and viscosities for top half of network using length (keeps length constant for all vessels in a level)

# Function for calculating resistance using diameter 
def microvessel_resistance(length, diameter):
    """ Calculate resistance for a given length and diameter. Length and diameter are in micrometres. """
    R = ((128*(estimated_viscosity(diameter)*1e-3)*length)/(np.pi*((diameter)**4))) * (0.00750062) *1e9
    return R

# For each microvessel level, find suitable diameter 
# Plot resistance against diameter for a range of diameters
for i in range(0+3, len(num_vessels_micro)+3):
    level = i
    resistance_value = realistic_resistance_matrix[0, level]
            
    level_length = realistic_length_matrix[0, level] * 1e3 #converted to micrometres 
    
    current_diameter = realistic_diameter_matrix[0, level] *1e3
    microvessel_diameters = np.arange(current_diameter, 100, 0.1) #micrometres 
    calculated_resistance_micro = np.array([])
    for n in range(0, len(microvessel_diameters)):
        calculated_resistance_micro = np.append(calculated_resistance_micro, microvessel_resistance(level_length, microvessel_diameters[n]))
    
    if (0):
        plt.figure()
        plt.plot(microvessel_diameters, calculated_resistance_micro)
        plt.xlabel('Diameter (\u03BCm)')
        plt.ylabel('Resistance (mmHgs/\u03BCl)')
        plt.title('Resistances for Microvessel Diameters')
        
    # Find diameter that gives closest required resistance value
    r_difference = np.absolute(calculated_resistance_micro - resistance_value) 
    min_diff_index = np.argmin(r_difference) #position of smallest difference
    d_micro = np.round(microvessel_diameters[min_diff_index], 1) * 1e-3 #convert back to mm 
    
    updated_d_micro = np.repeat(d_micro, num_vessels_micro[i-3]//2)
    realistic_diameter_matrix[0:num_vessels_micro[i-3]//2, i] = updated_d_micro

# Calculate corresponding viscosity 
for i in range(3, len(num_vessels_micro)+3):
    updated_d = realistic_diameter_matrix[0:num_vessels_micro[i-3]//2, i] * 1e3
    updated_visc_micro = np.round(estimated_viscosity(updated_d), 2)
    realistic_viscosity_matrix[0:num_vessels_micro[i-3]//2, i] = updated_visc_micro 

    
############################################################################


# %%
# Compartment 3: Veins 

num_vessels_ven = np.array([4, 2, 1])

start_pressure_ven = 24.
end_pressure_ven = 11. 

pressure_dist_ven = np.round(np.linspace(start_pressure_ven, end_pressure_ven, 4), 1) # Assume linear drop in pressure 
diameter_ven = np.array([0.36, 2.88, 4.5]) #mm
viscosity_ven = np.array([2.99, 2.99, 2.99])

# Store diameters in matrix
for i in range(11, len(diameter_ven)+11):
    d = np.repeat(diameter_ven[i-11], num_vessels_ven[i-11])
    realistic_diameter_matrix[0:num_vessels_ven[i-11], i] = d

# Store flow in matrix
n = 9
for i in range(11, 14):
    symmetrical_flow = realistic_flow_matrix[:, i-n]
    realistic_flow_matrix[:, i] = symmetrical_flow 
    n+=2 
    
# Calculate resistance using pressure difference and flow 
for i in range(11, len(diameter_ven)+11):
    flow_ven = realistic_flow_matrix[0:num_vessels_ven[i-11], i]
    resistance_ven = np.array([])
    for j in range(0, len(flow_ven)):
        r = calc_resistance(pressure_dist_ven[i-11], pressure_dist_ven[(i-11)+1], flow_ven[j])
        resistance_ven = np.append(resistance_ven, r)
    realistic_resistance_matrix[0:num_vessels_ven[i-11], i] = resistance_ven
 
# Calculate length from resistance
for i in range(11, len(diameter_ven)+11):
    r = realistic_resistance_matrix[0:num_vessels_ven[i-11], i]
    visc = np.repeat(viscosity_ven[i-11], num_vessels_ven[i-11])
    d = np.repeat(diameter_ven[i-11], num_vessels_ven[i-11])
    
    length_ven = np.array([])
    for j in range(0, len(r)):
        l = calc_length(r[j], visc[j], d[j])
        length_ven = np.append(length_ven, l)
    realistic_length_matrix[0:num_vessels_ven[i-11], i] = length_ven
    
# Set all vessels in a level to have the same length (use value for bottom half of network)
for i in range(11, len(diameter_ven)+11):
    realistic_length_matrix[0:num_vessels_ven[i-11]//2, i] = realistic_length_matrix[num_vessels_ven[i-11]//2:num_vessels_ven[i-11], i]

# Calculate diameters for top half of network using length (keeps length constistent for all vessels in a level)
for i in range(11, len(diameter_ven)+11):
    r = realistic_resistance_matrix[0:num_vessels_ven[i-11]//2, i]
    visc = np.repeat(viscosity_ven[i-11], num_vessels_ven[i-11]//2)
    l = realistic_length_matrix[0:num_vessels_ven[i-11]//2, i]
    
    diameter_a = np.array([])
    for j in range(0, len(r)):
        d = calc_diameter(r[j], visc[j], l[j])
        diameter_a = np.append(diameter_a, d)
    realistic_diameter_matrix[0:num_vessels_ven[i-11]//2, i] = diameter_a


realistic_viscosity_matrix[0, 11] = viscosity_ven[0]
realistic_viscosity_matrix[1, 11] = viscosity_ven[0]
realistic_viscosity_matrix[2, 11] = viscosity_ven[0]
realistic_viscosity_matrix[3, 11] = viscosity_ven[0]
realistic_viscosity_matrix[0, 12] = viscosity_ven[1]
realistic_viscosity_matrix[1, 12] = viscosity_ven[1]
realistic_viscosity_matrix[0, 13] = viscosity_ven[2]

# %% 
# Save top (path a) and bottom (path b) path values for diameter and length 

# Path a
diameter_path_a = realistic_diameter_matrix[0, :]
length_path_a = realistic_length_matrix[0, :]
viscosity_path_a = realistic_viscosity_matrix[0, :]

# Path b
path_b_index = [i-1 for i in branch_pattern]
diameter_path_b = np.array([])
length_path_b = np.array([])
viscosity_path_b = np.array([])
for n in range(0, len(path_b_index)):
    diameter_path_b = np.append(diameter_path_b, realistic_diameter_matrix[path_b_index[n], n])
    length_path_b = np.append(length_path_b, realistic_length_matrix[path_b_index[n], n])
    viscosity_path_b = np.append(viscosity_path_b, realistic_viscosity_matrix[path_b_index[n], n])
    
# Convert to micrometres - consistent units with Boas replication model 
diameter_path_a = diameter_path_a * 1e3
diameter_path_b = diameter_path_b * 1e3
length_path_a = length_path_a * 1e3
length_path_b = length_path_b * 1e3 

diameters = np.zeros([len(diameter_path_a), 2])
diameters[:, 0] = diameter_path_a
diameters[:, 1] = diameter_path_b

viscosities = np.zeros([len(viscosity_path_a), 2])
viscosities[:, 0] = viscosity_path_a
viscosities[:, 1] = viscosity_path_b

realistic_expected_pressure = np.concatenate((pressure_dist_art[:-1], pressure_dist_micro[1:-1], pressure_dist_ven))

if (1):
    print('Diameters for path a are: ' + str(diameter_path_a) + ' \u03BCm')
    print('Lengths for path a are: ' + str(length_path_a) + ' \u03BCm')
    print('Diameters for path b are: ' + str(diameter_path_b) + ' \u03BCm')
    print('Lengths for path b are: ' + str(length_path_b) + ' \u03BCm')
    print('Viscosities for path a are: ' + str(viscosity_path_a) + ' cP')
    print('Viscosities for path b are: ' + str(viscosity_path_b) + ' cP')
    print('ICA compliance is: \u03B2=' + str(beta_ICA))
    print('MCA compliance is: \u03B2=' + str(beta_MCA))