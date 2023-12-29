# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 11:07:18 2021

@author: kajal
"""

""" Class used to create vessel segments """
import numpy as np 
import numpy.matlib
from scipy.interpolate import interp1d

# Define pressures (mmHg)
minimum_pressure = 11. 
PIC = 7.

# Functions     
def parallel_resistance(resistance_0, resistance_1):
    """ Return total resistance for two resistors in parallel """ 
    total_resistance = (resistance_0 * resistance_1)/(resistance_0 + resistance_1)
    return total_resistance 


# Viscosity diameter plot (interpolation)
initial_diameter = np.array([30.5, 24.4, 19.5, 15.6, 12.5, 10., 8, 12., 15., 18.7, 23.4, 29.3, 36.6])
initial_viscosity = np.array([2.49, 2.34, 2.25, 2.20, 2.16, 2.12, 2.10, 2.15, 2.18, 2.22, 2.32, 2.51, 2.70]) #Convert to different units (x10**-3)
estimated_viscosity = interp1d(initial_diameter, initial_viscosity, fill_value="extrapolate") #in cP


class SingleSegment():
    """ Class for creating segments within a vessel """ 
    def __init__(self, input_pressure, output_pressure, diameter, length, viscosity, compliance, num_time_steps):
        self.input_pressure = input_pressure
        self.output_pressure = output_pressure
        self.diameter = diameter 
        self.length = length 
        self.viscosity = viscosity*1e-3
        self.compliance = compliance
        self.num_time_steps = num_time_steps 
        self.volume = np.pi * (self.diameter/2)**2 * self.length
        self.resistance = (8*self.viscosity*(self.length**3)*np.pi)/(self.volume**2) * (0.00750062) * 1e9 #convert from Pa to mmHg
#        self.resistance = self.resistance
        self.a0 =  self.volume/(self.input_pressure - PIC)**(1/self.compliance)
        
    def calculate_a0(self):
        """ Calculate a0 constant used in pressure-volume equation. """ 
        A_0 = self.volume/(self.input_pressure - PIC)**(1/self.compliance)
        self.a0 = A_0
        
    def change_pressure(self, new_pressure):
        """ Input new pressure, updates input pressure, volume, diameter and resistance for segment. """
        # if there is a change of pressure in the segment, update input_pressure, volume, resistance and diameter of segment
        if new_pressure != self.input_pressure:
            self.input_pressure = new_pressure
            
            updated_volume = self.a0 * ((self.input_pressure - PIC)**(1/self.compliance))
            self.volume = updated_volume 
            
            updated_diameter = (4*self.volume/np.pi*self.length)**(1/2)
            self.diameter = updated_diameter  # check units 
            
            # Only update viscosity (estimate_viscosity fn) for microvessels 
            if self.diameter <  200:
                updated_viscosity = estimated_viscosity(updated_diameter) * 1e-3
                self.viscosity = updated_viscosity
            else: 
                self.viscosity = self.viscosity

            updated_resistance = (8*self.viscosity*(self.length**3)*np.pi)/(self.volume**2)
            updated_resistance = updated_resistance * (0.00750062) *  1e9
            self.resistance = updated_resistance 
            
            
    def scale_pressures(self, remaining_resistance):
        # Calculate pressure drop across segment using remaining resistance in the network
        total_R = self.resistance + remaining_resistance
        pressure_drop = ((self.input_pressure - minimum_pressure)/total_R) * self.resistance
        perc_pressure_drop = pressure_drop/self.input_pressure
        scaling_factor = 1 - perc_pressure_drop 
        
        self.output_pressure = self.input_pressure * scaling_factor


class LongSegment(SingleSegment):
    """ Create segments within long segment/vessel."""
    def __init__(self, input_pressure, output_pressure, diameter, length, viscosity, compliance, num_time_steps):
        """ Initialise attributes of parent class """ 
        super().__init__(input_pressure, output_pressure, diameter, length, viscosity, compliance, num_time_steps)
        
        # number of time steps refers to how many time steps it takes to travel across the whole vessel, 
        # i.e how many segments in the vessel 
        
        self.input_pressure = np.array([self.input_pressure])
        self.output_pressure = np.array([self.output_pressure])
        self.resistance = np.array([self.resistance])
        self.volume = np.array([self.volume])
        
        self.diameter = np.array([self.diameter])
        
        self.viscosity = np.array([self.viscosity])
        
        self.flow = np.array([0])
        
        self.number_of_vessels = 1 
        
        # Create segments within long segment 
        self.seg = [] # list of individual segments within long segment 
        # create arrays for pressures, volumes and resistances in each individual segment within LongSegment object
        self.input_pressures = np.array([])
        self.output_pressures = np.array([])
        self.volumes = np.array([])
        self.diameters = np.array([])
        self.viscosities = np.array([])
        self.resistances = np.array([])
        for i in range(0, num_time_steps): 
            self.seg.append(SingleSegment(input_pressure, output_pressure, diameter, length/num_time_steps, viscosity, compliance, num_time_steps))
            self.input_pressures = np.append(self.input_pressures, self.seg[i].input_pressure)
            self.output_pressures = np.append(self.output_pressures, self.seg[i].output_pressure)
            self.volumes = np.append(self.volumes, self.seg[i].volume)
            self.diameters = np.append(self.diameters, self.seg[i].diameter)
            self.viscosities = np.append(self.viscosities, self.seg[i].viscosity)
            self.resistances = np.append(self.resistances, self.seg[i].resistance)
        
        self.flow_values = np.array([0])
        self.flow_matrix = np.zeros([1, num_time_steps]) #flow values for every segment within object 
        self.pressure_matrix = np.copy(self.input_pressures)
        self.pressure_matrix = np.append(self.pressure_matrix, self.output_pressure)
        
        self.flows = np.zeros([1, self.num_time_steps])
        self.resistances = np.reshape(self.resistances, (1, self.num_time_steps))
        
        self.pressures = np.append(self.input_pressures, self.output_pressure)
        
    def calculate_a0(self):
        """ Calculate a0 value for pressure, volume, compliance equation for every segment.""" 
        # Method to calculate a0 constant for every segment
        self.a0_values = np.array([])
        for i in range(0, self.num_time_steps):
            A_0 = self.volumes[i]/(self.input_pressures[i] - PIC)**(1/self.compliance)
            self.a0_values = np.append(self.a0_values, A_0)
        
        A_0 = self.volume/(self.input_pressure - PIC)**(1/self.compliance)
        self.a0 = A_0
                
#    def change_pressure(self, new_pressure):
#        # If new pressure input is different to current input pressure, update attributes 
#        self.input_pressures = self.input_pressures[0:-1] #update array of input pressures 
#        self.input_pressures = np.insert(self.input_pressures, 0, new_pressure)         
#        
#        if new_pressure != self.input_pressure: 
#            self.input_pressure = np.array([new_pressure]) #Input pressure is the input pressure for the vessel 
#
#        # Update values for individual segments within LongSegment 
#        # Calculate volumes and resistances of individual segments and find sum  
#        for i in range(0, len(self.input_pressures)):
#            self.seg[i].input_pressure = self.input_pressures[i]
#            self.seg[i].a0 = self.a0_values[i]
#            self.seg[i].volume = self.seg[i].a0 * ((self.seg[i].input_pressure - PIC)**(1/self.compliance))
#            self.volumes[i] = self.seg[i].volume
#            self.seg[i].diameter =  np.sqrt(((4*self.seg[i].volume)/(np.pi*self.seg[i].length)))
#            self.diameters[i] = self.seg[i].diameter
##            self.seg[i].viscosity = estimated_viscosity(self.seg[i].diameter) * 1e-3
#            self.viscosities[i] = self.seg[i].viscosity 
#            self.seg[i].resistance = (8*self.seg[i].viscosity*(self.seg[i].length**3)*np.pi)/(self.seg[i].volume**2) * (0.00750062) * 1e9
#            self.resistances[0][i] = self.seg[i].resistance
#       
#        total_updated_volume = np.sum(self.volumes)
#        self.volume = np.array([total_updated_volume])
#        
#        total_updated_resistance = np.sum(self.resistances)
#        self.resistance = np.array([total_updated_resistance])
#            
#        updated_diameter = np.sqrt(((4*self.volume)/(np.pi*self.length)))
#        self.diameter = updated_diameter  
#        
#        self.pressure_distribution = np.copy(self.input_pressures)
#        self.pressure_distribution = np.append(self.pressure_distribution, self.output_pressure)
#        self.pressure_matrix = np.vstack((self.pressure_matrix, self.pressure_distribution))
        
    def change_pressure(self, new_pressure):
        """ Input new pressure value for vessel, updates pressure, volume, diameter, viscosity, resistance attribute arrays. """
        # If new pressure input is different to current input pressure, update attributes 
        self.pressures = self.pressures[0:-1] #update array of input pressures 
        self.pressures = np.insert(self.pressures, 0, new_pressure)
                             
        self.input_pressures = np.copy(self.pressures)
        self.input_pressures = self.input_pressures[0:-1]
        self.output_pressures = np.copy(self.pressures)
        self.output_pressures = self.output_pressures[1:]
        
        self.input_pressure = np.array([self.input_pressures[0]])
        self.output_pressure = np.array([self.output_pressures[-1]])
        
        # Update values for individual segments within LongSegment 
        # Calculate volumes and resistances of individual segments and find sum  
        for i in range(0, len(self.input_pressures)):
            self.seg[i].input_pressure = self.input_pressures[i]
            self.seg[i].a0 = self.a0_values[i]
            self.seg[i].volume = self.seg[i].a0 * ((self.seg[i].input_pressure - PIC)**(1/self.compliance))
            self.volumes[i] = self.seg[i].volume
            self.seg[i].diameter =  np.sqrt(((4*self.seg[i].volume)/(np.pi*self.seg[i].length)))
            self.diameters[i] = self.seg[i].diameter
            
            
            # Only change viscosity values for microvessels, use estimated_viscosity() 
            if self.diameters[i] <  200:
                self.seg[i].viscosity = estimated_viscosity(self.seg[i].diameter) * 1e-3
                self.viscosities[i] = self.seg[i].viscosity 
            else:
                self.viscosities[i] = self.seg[i].viscosity
            
            self.seg[i].resistance = (8*self.seg[i].viscosity*(self.seg[i].length**3)*np.pi)/(self.seg[i].volume**2) * (0.00750062) * 1e9
            self.resistances[0][i] = self.seg[i].resistance
        
        total_updated_volume = np.sum(self.volumes)
        self.volume = np.array([total_updated_volume])
        
        total_updated_resistance = np.sum(self.resistances)
        self.resistance = np.array([total_updated_resistance])
            
        updated_diameter = np.sqrt(((4*self.volume)/(np.pi*self.length)))
        self.diameter = updated_diameter  
        
        self.pressure_distribution = np.copy(self.input_pressures)
        self.pressure_distribution = np.append(self.pressure_distribution, self.output_pressure)
        self.pressure_matrix = np.vstack((self.pressure_matrix, self.pressure_distribution))
                
        self.viscosity = self.viscosities[0]
        
    def calculate_flow_branch(self):
        self.flow_values = np.array([])
        self.pressure_distribution = np.copy(self.input_pressures)
        self.pressure_distribution = np.append(self.pressure_distribution, self.output_pressure)
        for i in range(0, self.num_time_steps):
            pressure_gradient = self.pressure_distribution[i] - self.pressure_distribution[i+1]
            flow_in_segment = pressure_gradient/self.resistances[i] #* 1e-9 #convert to m^3/s
            self.flow_values = np.append(self.flow_values, flow_in_segment) 
            
        self.flow_matrix = np.vstack((self.flow_matrix, self.flow_values))
        self.pressure_matrix = np.vstack((self.pressure_matrix, self.pressure_distribution))
        
        # output flow from middle segment
        #self.flow = np.array([self.flow_values[self.num_time_steps//2]])
        self.flow = np.array([self.flow_values[-1]]) #returns flow value for last segment in vessel 
         
    def scale_pressures(self, remaining_resistance):
        self.scaling_factor_ls = np.ones([len(self.input_pressures)]) #create array to store scaling factor values for each segment (acts on input pressure)
        for n in range(0, len(self.input_pressures)):
            total_remaining_resistance = np.sum(self.resistances[n:]) + remaining_resistance[0] # multiply sum of remaining resistance in current vessel by corresponding branch pattern value? 
            pressure_drop = ((self.input_pressures[n] - minimum_pressure)/total_remaining_resistance) * self.resistances[n]
            perc_pressure_drop = pressure_drop/self.input_pressures[n]
            self.scaling_factor_ls[n] = (1-perc_pressure_drop)
#            print(pressure_drop)
#            print(self.scaling_factor_ls[n])
#            print(self.resistances[n])
#            print(total_remaining_resistance)
        for n in range(0, len(self.input_pressures)):
            self.input_pressures[n] = self.input_pressures[n] * self.scaling_factor_ls[n]
        self.output_pressure = np.array([self.input_pressures[-1]])

    def update_flows(self):
        """ Updates self.flow attribute using vessel's middle segment value of flow. """
        #self.flow = np.array([self.flows[0][self.num_time_steps//2]]) #Set flow value for the object as middle segment value
        self.flow = np.array([self.flows[0][0]]) #first segment flow
        self.flow_values = np.append(self.flow_values, self.flow) #Add object flow value to array flow values
        self.flow_matrix = np.vstack((self.flow_matrix, self.flows)) #save flow value for each segment at every time point
       
    def calculate_pressures(self):
#        for n in range(0, len(self.pressures)-1):
#            if n == 0:
#                self.pressures = self.pressures[0:-1] #update array of input pressures 
#                self.pressures = np.insert(self.pressures, 0, new_pressure)
#                pressure_gradient = self.flows[0][0]*self.resistances[0][0]
#                P_out = self.pressures[0]-pressure_gradient
#                self.pressures[n+1] = P_out
#            else: 
#                pressure_gradient = self.flows[0][n]*self.resistances[0][n]
#                P_out = self.pressures[n] - pressure_gradient
#                self.pressures[n+1] = P_out
                
        for n in range(0, len(self.pressures)-1):
            pressure_gradient = self.flows[0][n]*self.resistances[0][n]
            P_out = self.pressures[n] - pressure_gradient
            self.pressures[n] = P_out # replace input pressures in self.pressures with output pressure for corresponding segment
        
        self.output_pressure = np.array([self.pressures[-2]]) #set object output pressure as value before the last value in self.pressures
        
#        self.input_pressures = np.copy(self.pressures)
#        self.input_pressures = self.input_pressures[0:-1]
#        self.output_pressures = np.copy(self.pressures)
#        self.output_pressures = self.output_pressures[1:]


#        # Update values for individual segments within LongSegment 
#        # Calculate volumes and resistances of individual segments and find sum  
#        for i in range(0, len(self.input_pressures)):
#            self.seg[i].input_pressure = self.input_pressures[i]
#            self.seg[i].a0 = self.a0_values[i]
#            self.seg[i].volume = self.seg[i].a0 * ((self.seg[i].input_pressure - PIC)**(1/self.compliance))
#            self.volumes[i] = self.seg[i].volume
#            self.seg[i].diameter =  np.sqrt(((4*self.seg[i].volume)/(np.pi*self.seg[i].length)))
#            self.diameters[i] = self.seg[i].diameter
##            self.seg[i].viscosity = estimated_viscosity(self.seg[i].diameter) * 1e-3
#            self.viscosities[i] = self.seg[i].viscosity 
#            self.seg[i].resistance = (8*self.seg[i].viscosity*(self.seg[i].length**3)*np.pi)/(self.seg[i].volume**2) * (0.00750062) * 1e9
#            self.resistances[0][i] = self.seg[i].resistance            
#
#        total_updated_volume = np.sum(self.volumes)
#        self.volume = np.array([total_updated_volume])
#        
#        total_updated_resistance = np.sum(self.resistances)
#        self.resistance = np.array([total_updated_resistance])
#            
#        updated_diameter = np.sqrt(((4*self.volume)/(np.pi*self.length)))
#        self.diameter = updated_diameter  
#        
#        self.pressure_matrix = np.vstack((self.pressure_matrix, self.pressures))
           
        


class DivergingBranchSegment(SingleSegment): 
    """ Class for creating closed branched (1-2) segment """ 
    def __init__(self, input_pressure, output_pressure, diameter, length, viscosity, compliance, num_time_steps):
        """ Initialise attributes of parent class """ 
        super().__init__(input_pressure, output_pressure, diameter, length, viscosity, compliance, num_time_steps)
        
        self.diameter_a = diameter[0]
        self.diameter_b = diameter[1]
        
        self.compliance_a = compliance[0]
        self.compliance_b = compliance[1]
        
        self.volume_a = np.pi * (self.diameter_a/2)**2 * self.length #total volume in vessel a of diverging segment
        self.volume_b = np.pi * (self.diameter_b/2)**2 * self.length #total volume in vessel b of diverging segment
        
#        self.diameter = np.array([self.diameter_a, self.diameter_b])

        self.flow = np.array([0, 0]) #hold current value of flow in object for vessels a and b 
        
        self.number_of_vessels = 2
                
        # Create segments within vessels a and b in branch 
        self.seg_a = []
        self.seg_b = [] 
        self.input_pressures_a = np.array([])
        self.input_pressures_b = np.array([])
        self.output_pressures_a = np.array([])
        self.output_pressures_b = np.array([])
        self.volumes_a = np.array([])
        self.volumes_b = np.array([])
        self.diameters_a = np.array([])
        self.diameters_b = np.array([])
        self.viscosities_a = np.array([])
        self.viscosities_b = np.array([])
        self.resistances_a = np.array([])
        self.resistances_b = np.array([])
        for i in range(0, num_time_steps):
            self.seg_a.append(SingleSegment(input_pressure, output_pressure, self.diameter_a, length/num_time_steps, viscosity, self.compliance_a, num_time_steps/num_time_steps))
            self.seg_b.append(SingleSegment(input_pressure, output_pressure, self.diameter_b, length/num_time_steps, viscosity, self.compliance_b, num_time_steps/num_time_steps))
            self.input_pressures_a = np.append(self.input_pressures_a, self.seg_a[i].input_pressure)
            self.input_pressures_b = np.append(self.input_pressures_b, self.seg_b[i].input_pressure)
            self.output_pressures_a = np.append(self.output_pressures_a, self.seg_a[i].output_pressure)
            self.output_pressures_b = np.append(self.output_pressures_b, self.seg_b[i].output_pressure)
            self.volumes_a = np.append(self.volumes_a, self.seg_a[i].volume)
            self.volumes_b = np.append(self.volumes_b, self.seg_b[i].volume)
            self.diameters_a = np.append(self.diameters_a, self.seg_a[i].diameter)
            self.diameters_b = np.append(self.diameters_b, self.seg_b[i].diameter)
            self.viscosities_a = np.append(self.viscosities_a, self.seg_a[i].viscosity)
            self.viscosities_b = np.append(self.viscosities_b, self.seg_b[i].viscosity)
            self.resistances_a = np.append(self.resistances_a, self.seg_a[i].resistance)
            self.resistances_b = np.append(self.resistances_b, self.seg_b[i].resistance)
        
        self.resistances = np.array([self.resistances_a, self.resistances_b])
        
        self.flows_a = np.zeros([self.num_time_steps])
        self.flows_b = np.zeros([self.num_time_steps])
        self.flows = np.array([self.flows_a, self.flows_b])
        
        self.flow_values_a = np.array([0])
        self.flow_values_b = np.array([0])
        
        # 2 output pressures, one for each vessel in branch
        self.output_pressure = np.array([self.output_pressures_a[-1], self.output_pressures_b[-1]])
        self.input_pressure = np.array([self.input_pressures_a[0], self.input_pressures_b[0]])
        
        self.pressures_a = np.append(self.input_pressures_a, self.output_pressures_a[-1])
        self.pressures_b = np.append(self.input_pressures_b, self.output_pressures_b[-1])

        # Store pressure and flow values across each segment for each branch within Diverging Branch Segment object 
        self.pressure_matrix_a = np.zeros([1, num_time_steps])
        self.pressure_matrix_b = np.zeros([1, num_time_steps])
        self.pressure_matrix_a[0, :] = self.input_pressures_a #initial (t=0) input pressures 
        self.pressure_matrix_b[0, :] = self.input_pressures_b
        
        self.flow_matrix_a = np.zeros([1, num_time_steps])
        self.flow_matrix_b = np.zeros([1, num_time_steps])

        # Sum of resistances, one value for each branch 
        self.resistance = np.array([np.sum(self.resistances_a), np.sum(self.resistances_b)])
        self.volume = np.array([np.sum(self.volumes_a), np.sum(self.volumes_b)])
        self.viscosity = np.array([self.viscosities_a[0], self.viscosities_b[0]])

    def calculate_a0(self):
        """ Calculate a0 value for pressure, volume, compliance equation for every segment for vessels a and b."""
        self.a0_values_a = np.array([])
        self.a0_values_b = np.array([])
        
        for i in range(0, self.num_time_steps):
            A_0_a = self.volumes_a[i]/(self.input_pressures_a[i] - PIC)**(1/self.compliance_a)
            self.a0_values_a = np.append(self.a0_values_a, A_0_a)
            
            A_0_b = self.volumes_b[i]/(self.input_pressures_b[i] - PIC)**(1/self.compliance_b)
            self.a0_values_b = np.append(self.a0_values_b, A_0_b)
        
        A_0 = self.volume/(self.input_pressure - PIC)**(1/self.compliance)
        self.a0 = A_0

    def change_pressure(self, new_pressure): # 2 values for new_pressure, one for each vessel     
        """ Input new pressure, updates input pressure, volume, diameter, viscosity and resistance for vessels a and b. """                                        
        self.pressures_a = self.pressures_a[0:-1]
        self.pressures_a = np.insert(self.pressures_a, 0, new_pressure[0])  
        
        self.pressures_b = self.pressures_b[0:-1]
        self.pressures_b = np.insert(self.pressures_b, 0, new_pressure[1])     
        
        self.pressure_matrix_a = np.vstack((self.pressure_matrix_a, self.input_pressures_a))
        self.pressure_matrix_b = np.vstack((self.pressure_matrix_b, self.input_pressures_b))
        
        self.input_pressure_a = np.array(self.pressures_a[0])
        self.output_pressure_a = np.array(self.pressures_a[-1])
        self.input_pressure_b = np.array(self.pressures_b[0])
        self.output_pressure_b = np.array(self.pressures_b[-1])
        
        self.input_pressures_a = np.copy(self.pressures_a)
        self.input_pressures_a = self.input_pressures_a[0:-1]
        
        self.input_pressures_b = np.copy(self.pressures_b)
        self.input_pressures_b = self.input_pressures_b[0:-1]
        
        self.input_pressure = np.array([self.input_pressure_a, self.input_pressure_b])
        self.output_pressure = np.array([self.output_pressure_a, self.output_pressure_b])
        
        self.pressures = np.array([self.pressures_a, self.pressures_b])


        # Calculate volumes and resistances of individual segments and find sum  
        for i in range(0, len(self.input_pressures_a)):
            self.seg_a[i].input_pressure = self.input_pressures_a[i]
            self.seg_a[i].a0 = self.a0_values_a[i]
            self.seg_a[i].volume = self.seg_a[i].a0 * ((self.seg_a[i].input_pressure - PIC)**(1/self.compliance_a))
            self.volumes_a[i] = self.seg_a[i].volume
            self.seg_a[i].diameter =  np.sqrt(((4*self.seg_a[i].volume)/(np.pi*self.seg_a[i].length)))
            self.diameters_a[i] = self.seg_a[i].diameter
                        
            # Only update viscosity (estimate_viscosity function) for microvessels
            if self.diameters_a[i] <  200:
                self.seg_a[i].viscosity = estimated_viscosity(self.seg_a[i].diameter) * 1e-3
                self.viscosities_a[i] = self.seg_a[i].viscosity 
            else:
                self.viscosities_a[i] = self.seg_a[i].viscosity 
                
            self.seg_a[i].resistance = (8*self.seg_a[i].viscosity*(self.seg_a[i].length**3)*np.pi)/(self.seg_a[i].volume**2) * (0.00750062) * 1e9
            self.resistances_a[i] = self.seg_a[i].resistance
            
            self.seg_b[i].input_pressure = self.input_pressures_b[i]
            self.seg_b[i].a0 = self.a0_values_b[i]
            self.seg_b[i].volume = self.seg_b[i].a0 * ((self.seg_b[i].input_pressure - PIC)**(1/self.compliance_b))
            self.volumes_b[i] = self.seg_b[i].volume
            self.seg_b[i].diameter =  np.sqrt(((4*self.seg_b[i].volume)/(np.pi*self.seg_b[i].length)))
            self.diameters_b[i] = self.seg_b[i].diameter
            
            # Only update viscosity (estimate_viscosity function) for microvessels
            if self.diameters_b[i] <  200:
                self.seg_b[i].viscosity = estimated_viscosity(self.seg_b[i].diameter) * 1e-3
                self.viscosities_b[i] = self.seg_b[i].viscosity 
            else:
                self.viscosities_b[i] = self.seg_b[i].viscosity 
                
            self.seg_b[i].resistance = (8*self.seg_b[i].viscosity*(self.seg_b[i].length**3)*np.pi)/(self.seg_b[i].volume**2) * (0.00750062) * 1e9
            self.resistances_b[i] = self.seg_b[i].resistance

       
        total_updated_volume_a = np.sum(self.volumes_a)
        total_updated_volume_b = np.sum(self.volumes_b)
        self.volume = np.array([total_updated_volume_a, total_updated_volume_b])
        self.resistances = np.array([self.resistances_a, self.resistances_b])
        
#        total_updated_resistance = np.sum(self.resistances_a) + np.sum(self.resistances_b)
#        self.resistance = total_updated_resistance
        
        total_updated_resistance_a = np.sum(self.resistances_a)
        total_updated_resistance_b = np.sum(self.resistances_b)
        self.resistance = np.array([total_updated_resistance_a, total_updated_resistance_b])
            
        updated_diameter_a = np.sqrt(((4*total_updated_volume_a)/(np.pi*self.length)))
        updated_diameter_b = np.sqrt(((4*total_updated_volume_b)/(np.pi*self.length)))
        self.diameter = np.array([updated_diameter_a, updated_diameter_b])  # check units 
        
    
    def scale_pressures(self, remaining_resistance):
        self.scaling_factor_a = np.ones([len(self.input_pressures_a)])
        self.scaling_factor_b = np.ones([len(self.input_pressures_b)])

        for n in range(0, len(self.input_pressures_a)):
#            current_resistance_a = np.sum(self.resistances_a[n:])
#            current_resistance_b = np.sum(self.resistances_b[n:])
#            current_branch_resistance = parallel_resistance(current_resistance_a, current_resistance_b)
            
#            total_remaining_resistance_a = current_branch_resistance + remaining_resistance[0]
            total_remaining_resistance_a = np.sum(self.resistances_a[n:]) + remaining_resistance[0]
            pressure_drop_a = ((self.input_pressures_a[n] - minimum_pressure)/total_remaining_resistance_a) * self.resistances_a[n]
            perc_pressure_drop_a = pressure_drop_a/self.input_pressures_a[n]
            self.scaling_factor_a[n] = (1-perc_pressure_drop_a)
#            print(self.scaling_factor_a[n])
#            print(self.resistances_a[n])
#            print(total_remaining_resistance_a)
            
#            total_remaining_resistance_b = current_branch_resistance + remaining_resistance[1]
            total_remaining_resistance_b = np.sum(self.resistances_b[n:]) + remaining_resistance[1]
            pressure_drop_b = ((self.input_pressures_b[n] - minimum_pressure)/total_remaining_resistance_b) * self.resistances_b[n]
            perc_pressure_drop_b = pressure_drop_b/self.input_pressures_b[n]
            self.scaling_factor_b[n] = (1-perc_pressure_drop_b)            
            
        for n in range(0, len(self.input_pressures_a)):
            self.input_pressures_a[n] = self.input_pressures_a[n] * self.scaling_factor_a[n]
            self.input_pressures_b[n] = self.input_pressures_b[n] * self.scaling_factor_b[n]
        
        self.output_pressure = np.array([self.input_pressures_a[-1], self.input_pressures_b[-1]])
        
    def calculate_flow_branch(self):
        self.flow_values_a = np.array([])
        self.flow_values_b = np.array([])
        self.pressure_dist_a = np.copy(self.input_pressures_a)
        self.pressure_dist_a = np.append(self.pressure_dist_a, self.output_pressure[0])
        self.pressure_dist_b = np.copy(self.input_pressures_b)
        self.pressure_dist_b = np.append(self.pressure_dist_b, self.output_pressure[1])
        for t in range(0, self.num_time_steps):
            pressure_gradient_a = self.pressure_dist_a[t] - self.pressure_dist_a[t+1]
            flow_in_segment_a = pressure_gradient_a/self.resistances_a[t] #* 1e-9
            self.flow_values_a = np.append(self.flow_values_a, flow_in_segment_a) 
            
            pressure_gradient_b = self.pressure_dist_b[t] - self.pressure_dist_b[t+1]
            flow_in_segment_b = pressure_gradient_b/self.resistances_b[t] #* 1e-9
            self.flow_values_b = np.append(self.flow_values_b, flow_in_segment_b) 

        self.flow_matrix_a = np.vstack((self.flow_matrix_a, self.flow_values_a))
        self.flow_matrix_b = np.vstack((self.flow_matrix_b, self.flow_values_b))
        
        #self.flow = np.array([self.flow_values_a[self.num_time_steps//2], self.flow_values_b[self.num_time_steps//2]])
        self.flow = np.array([self.flow_values_a[-1], self.flow_values_b[-1]])
        
    def update_flows(self):
        """ Update self.flow attribute, use middle segment flow values for vessels a and b. """ 
        self.flows_a = np.copy(self.flows[0])
        self.flow_values_a = np.append(self.flow_values_a, self.flows_a[self.num_time_steps//2]) #pick middle segment flow value

        self.flows_b = np.copy(self.flows[1])
        self.flow_values_b = np.append(self.flow_values_b, self.flows_b[self.num_time_steps//2])        
        
        self.flow_matrix_a = np.vstack((self.flow_matrix_a, self.flows_a)) #add segment flow values to matrix 
        self.flow_matrix_b = np.vstack((self.flow_matrix_b, self.flows_b)) 
        
        #self.flow = np.array([self.flows_a[self.num_time_steps//2], self.flows_b[self.num_time_steps//2]])
        self.flow = np.array([self.flows_a[0], self.flows_b[0]]) #first segment flow

    def calculate_pressures(self):
        self.flows_a = self.flows[0]
        self.flows_b = self.flows[1]
#        for n in range(0, len(self.input_pressures_a)-1):
#            if n == 0:
#                self.input_pressures_a = self.input_pressures_a[0:-1]
#                self.input_pressures_a = np.insert(self.input_pressures_a, 0, new_pressure[0])
#                pressure_gradient_a  = self.flows_a[0] * self.resistances_a[0]
#                P_out_a = self.input_pressures_a[0] - pressure_gradient_a
#                self.input_pressures_a[n+1] = P_out_a 
#            
#                self.input_pressures_b = self.input_pressures_b[0:-1]
#                self.input_pressures_b = np.insert(self.input_pressures_b, 0, new_pressure[1])
#                pressure_gradient_b = self.flows_b[0]*self.flows_b[0]
#                P_out_b = self.input_pressures_b[0] - pressure_gradient_b
#                self.input_pressures_b[n+1] = P_out_b
#                
#            else: 
#                pressure_gradient_a = self.flows_a[n] * self.resistances_a[n]
#                P_out_a = self.input_pressures_a[n] - pressure_gradient_a
#                self.input_pressures_a[n+1] = P_out_a
#                
#                pressure_gradient_b = self.flows_b[n] * self.resistances_b[n]
#                P_out_b = self.input_pressures_b[n] - pressure_gradient_b
#                self.input_pressures_b[n+1] = P_out_b 
#                
#            self.output_pressure = np.array([self.input_pressures_a[-1], self.input_pressures_b[-1]])
#            self.input_pressure = np.array([self.input_pressures_a[0], self.input_pressures_b[0]])
        
#        for n in range(0, len(self.pressures_a)-1):
#            if n == 0: 
#                self.pressures_a = self.pressures_a[0:-1]
#                self.pressures_a = np.insert(self.pressures_a, 0, new_pressure[0])
#                pressure_gradient_a = self.flows[0][0]*self.resistances[0][0]
#                P_out_a = self.pressures_a[0]-pressure_gradient_a
#                self.pressures_a[n+1] = P_out_a 
#                
#                self.pressures_b = self.pressures_b[0:-1]
#                self.pressures_b = np.insert(self.pressures_b, 0, new_pressure[1])
#                pressure_gradient_b = self.flows[1][0] * self.resistances[1][0]
#                P_out_b = self.pressures_b[0]-pressure_gradient_b
#                self.pressures_b[n+1] = P_out_b 
#                
#            else: 
#                pressure_gradient_a = self.flows[0][n]*self.resistances[0][n]
#                P_out_a = self.pressures_a[n] - pressure_gradient_a
#                self.pressures_a[n+1] = P_out_a 
#                
#                pressure_gradient_b = self.flows[1][n]*self.resistances[1][n]
#                P_out_b = self.pressures_b[n] - pressure_gradient_b 
#                self.pressures_b[n+1] = P_out_b 
#                
#        self.output_pressure_a = np.array(self.pressures_a[-1])
#        self.input_pressure_a = np.array(self.pressures_a[0])
#        self.output_pressure_b = np.array(self.pressures_b[-1])
#        self.input_pressure_b = np.array(self.pressures_b[0])
#        
#        self.input_pressures_a = np.copy(self.pressures_a)
#        self.input_pressures_a = self.input_pressures_a[0:-1]
#        
#        self.input_pressure = np.array([self.input_pressure_a, self.input_pressure_b])
#        self.output_pressure = np.array([self.output_pressure_a, self.output_pressure_b])
#        
#        self.pressures = np.array([self.pressures_a, self.pressures_b])
#
#        # Calculate volumes and resistances of individual segments and find sum  
#        for i in range(0, len(self.input_pressures_a)):
#            self.seg_a[i].input_pressure = self.input_pressures_a[i]
#            self.seg_a[i].a0 = self.a0_values_a[i]
#            self.seg_a[i].volume = self.seg_a[i].a0 * ((self.seg_a[i].input_pressure - PIC)**(1/self.compliance_a))
#            self.volumes_a[i] = self.seg_a[i].volume
#            self.seg_a[i].diameter =  np.sqrt(((4*self.seg_a[i].volume)/(np.pi*self.seg_a[i].length)))
#            self.diameters_a[i] = self.seg_a[i].diameter
##            self.seg_a[i].viscosity = estimated_viscosity(self.seg_a[i].diameter) * 1e-3
#            self.viscosities_a[i] = self.seg_a[i].viscosity 
#            self.seg_a[i].resistance = (8*self.seg_a[i].viscosity*(self.seg_a[i].length**3)*np.pi)/(self.seg_a[i].volume**2) * (0.00750062) * 1e9
#            self.resistances_a[i] = self.seg_a[i].resistance
#            
#            self.seg_b[i].input_pressure = self.input_pressures_b[i]
#            self.seg_b[i].a0 = self.a0_values_b[i]
#            self.seg_b[i].volume = self.seg_b[i].a0 * ((self.seg_b[i].input_pressure - PIC)**(1/self.compliance_b))
#            self.volumes_b[i] = self.seg_b[i].volume
#            self.seg_b[i].diameter =  np.sqrt(((4*self.seg_b[i].volume)/(np.pi*self.seg_b[i].length)))
#            self.diameters_b[i] = self.seg_b[i].diameter
##            self.seg_b[i].viscosity = estimated_viscosity(self.seg_b[i].diameter) * 1e-3
#            self.viscosities_b[i] = self.seg_b[i].viscosity 
#            self.seg_b[i].resistance = (8*self.seg_b[i].viscosity*(self.seg_b[i].length**3)*np.pi)/(self.seg_b[i].volume**2) * (0.00750062) * 1e9
#            self.resistances_b[i] = self.seg_b[i].resistance
#
#        total_updated_volume_a = np.sum(self.volumes_a)
#        total_updated_volume_b = np.sum(self.volumes_b)
#        self.volume = np.array([total_updated_volume_a, total_updated_volume_b])
#        self.resistances = np.array([self.resistances_a, self.resistances_b])
#        
##        total_updated_resistance = np.sum(self.resistances_a) + np.sum(self.resistances_b)
##        self.resistance = total_updated_resistance
#        
#        total_updated_resistance_a = np.sum(self.resistances_a)
#        total_updated_resistance_b = np.sum(self.resistances_b)
#        self.resistance = np.array([total_updated_resistance_a, total_updated_resistance_b])
#            
#        updated_diameter_a = np.sqrt(((4*total_updated_volume_a)/(np.pi*self.length)))
#        updated_diameter_b = np.sqrt(((4*total_updated_volume_b)/(np.pi*self.length)))
#        self.diameter = np.array([updated_diameter_a, updated_diameter_b])  # check units 
        
        for n in range(0, len(self.pressures_a)-1):
            pressure_gradient_a = self.flows_a[n] * self.resistances_a[n]
            P_out_a = self.pressures_a[n] - pressure_gradient_a
            self.pressures_a[n] = P_out_a
            
            pressure_gradient_b = self.flows_b[n] * self.resistances_b[n]
            P_out_b = self.pressures_b[n] - pressure_gradient_b
            self.pressures_b[n] = P_out_b             
        
        self.output_pressure = np.array([self.pressures_a[-2], self.pressures_b[-2]])
        
            
class ConvergingBranchSegment(SingleSegment): 
    """ Class for creating closed branched (1-2) segment """ 
    def __init__(self, input_pressure, output_pressure, diameter, length, viscosity, compliance, num_time_steps):
        """ Initialise attributes of parent class """ 
        super().__init__(input_pressure, output_pressure, diameter, length, viscosity, compliance, num_time_steps)
        
        self.output_pressure = np.array([self.output_pressure, self.output_pressure]) #output pressure per vessel 
 
        self.diameter_a = diameter[0]
        self.diameter_b = diameter[1]
                                                 #to work with x, y, z indexing 
        self.compliance_a = compliance[0]
        self.compliance_b = compliance[1]
        
        self.volume_a = np.pi * (self.diameter_a/2)**2 * self.length #total volume in vessel a of diverging segment
        self.volume_b = np.pi * (self.diameter_b/2)**2 * self.length #total volume in vessel b of diverging segment

#        self.diameter = np.array([self.diameter_a, self.diameter_b])

        self.flow = np.array([0, 0])

        self.number_of_vessels = 2
        
        # Create segments within branches a and b
        self.seg_a = []
        self.seg_b = [] 
        self.input_pressures_a = np.array([])
        self.input_pressures_b = np.array([])
        self.output_pressures_a = np.array([])
        self.output_pressures_b = np.array([])
        self.volumes_a = np.array([])
        self.volumes_b = np.array([])
        self.diameters_a = np.array([])
        self.diameters_b = np.array([])
        self.viscosities_a = np.array([])
        self.viscosities_b = np.array([])
        self.resistances_a = np.array([])
        self.resistances_b = np.array([])
        for i in range(0, num_time_steps):
            self.seg_a.append(SingleSegment(input_pressure, output_pressure, self.diameter_a, length/num_time_steps, viscosity, self.compliance_a, num_time_steps/num_time_steps))
            self.seg_b.append(SingleSegment(input_pressure, output_pressure, self.diameter_b, length/num_time_steps, viscosity, self.compliance_b, num_time_steps/num_time_steps))
            self.input_pressures_a = np.append(self.input_pressures_a, self.seg_a[i].input_pressure)
            self.input_pressures_b = np.append(self.input_pressures_b, self.seg_b[i].input_pressure)
            self.output_pressures_a = np.append(self.output_pressures_a, self.seg_a[i].output_pressure)
            self.output_pressures_b = np.append(self.output_pressures_b, self.seg_b[i].output_pressure)
            self.volumes_a = np.append(self.volumes_a, self.seg_a[i].volume)
            self.volumes_b = np.append(self.volumes_b, self.seg_b[i].volume)
            self.diameters_a = np.append(self.diameters_a, self.seg_a[i].diameter)
            self.diameters_b = np.append(self.diameters_b, self.seg_b[i].diameter)
            self.viscosities_a = np.append(self.viscosities_a, self.seg_a[i].viscosity)
            self.viscosities_b = np.append(self.viscosities_b, self.seg_b[i].viscosity)
            self.resistances_a = np.append(self.resistances_a, self.seg_a[i].resistance)
            self.resistances_b = np.append(self.resistances_b, self.seg_b[i].resistance)
            
    
        self.resistances = np.array([self.resistances_a, self.resistances_b])
        
        self.flows_a = np.zeros([self.num_time_steps])
        self.flows_b = np.zeros([self.num_time_steps])
        self.flows = np.array([self.flows_a, self.flows_b])
        self.flow_values_a = np.array([0])
        self.flow_values_b = np.array([0])

        self.pressure_matrix_a = np.zeros([1, num_time_steps])
        self.pressure_matrix_b = np.zeros([1, num_time_steps])
        self.pressure_matrix_a[0, :] = self.input_pressures_a
        self.pressure_matrix_b[0, :] = self.input_pressures_b
        
        self.flow_matrix_a = np.zeros([1, num_time_steps])
        self.flow_matrix_b = np.zeros([1, num_time_steps])
        
        self.resistance = np.array([np.sum(self.resistances_a), np.sum(self.resistances_b)])
        self.volume = np.array([np.sum(self.volumes_a), np.sum(self.volumes_b)])

        self.input_pressure = np.array([self.input_pressures_a[0], self.input_pressures_b[0]])
        self.output_pressure = np.array([self.output_pressures_a[-1], self.output_pressures_b[-1]])
        
        self.pressures_a = np.append(self.input_pressures_a, self.output_pressures_a[-1])
        self.pressures_b = np.append(self.input_pressures_b, self.output_pressures_b[-1])
        
    def calculate_a0(self):
        """ Calculate a0 value for pressure, volume, compliance equation for every segment for vessels a and b."""
        self.a0_values_a = np.array([])
        self.a0_values_b = np.array([])
        for i in range(0, self.num_time_steps):
            A_0_a = self.volumes_a[i]/(self.input_pressures_a[i] - PIC)**(1/self.compliance_a)
            self.a0_values_a = np.append(self.a0_values_a, A_0_a)
            
            A_0_b = self.volumes_b[i]/(self.input_pressures_b[i] - PIC)**(1/self.compliance_b)
            self.a0_values_b = np.append(self.a0_values_b, A_0_b)
        
        A_0 = self.volume/(self.input_pressure - PIC)**(1/self.compliance)
        self.a0 = A_0

    def change_pressure(self, new_pressure):   # two values for new_pressure
        """ Input new pressure, updates input pressure, volume, diameter, viscosity and resistance for vessels a and b. """                                        

        self.pressures_a = self.pressures_a[0:-1]
        self.pressures_a = np.insert(self.pressures_a, 0, new_pressure[0])  
        
        self.pressures_b = self.pressures_b[0:-1]
        self.pressures_b = np.insert(self.pressures_b, 0, new_pressure[1])  
        
        self.pressure_matrix_a = np.vstack((self.pressure_matrix_a, self.input_pressures_a))
        self.pressure_matrix_b = np.vstack((self.pressure_matrix_b, self.input_pressures_b))

        self.input_pressures_a = np.copy(self.pressures_a)
        self.input_pressures_a = self.input_pressures_a[0:-1]
        self.output_pressures_a = np.copy(self.pressures_a)
        self.output_pressures_a = self.output_pressures_a[1:]

        self.input_pressures_b = np.copy(self.pressures_b)
        self.input_pressures_b = self.input_pressures_b[0:-1]
        self.output_pressures_b = np.copy(self.pressures_b)
        self.output_pressures_b = self.output_pressures_b[1:]
        
        self.input_pressure = np.array([self.pressures_a[0], self.pressures_b[0]]) 
        output_pressure = (((self.pressures_a[-1]*self.volumes_a[-1]) + (self.pressures_b[-1]*self.volumes_b[-1]))/(self.volumes_a[-1] + self.volumes_b[-1]))
        self.output_pressure = np.array([output_pressure, output_pressure])
    
        self.pressures = np.array([self.pressures_a, self.pressures_b])


        # Calculate volumes and resistances of individual segments and find sum  
        for i in range(0, len(self.input_pressures_a)):
            self.seg_a[i].input_pressure = self.input_pressures_a[i]
            self.seg_a[i].a0 = self.a0_values_a[i]
            self.seg_a[i].volume = self.seg_a[i].a0 * ((self.seg_a[i].input_pressure - PIC)**(1/self.compliance_a))
            self.volumes_a[i] = self.seg_a[i].volume
            self.seg_a[i].diameter =  np.sqrt(((4*self.seg_a[i].volume)/(np.pi*self.seg_a[i].length)))
            self.diameters_a[i] = self.seg_a[i].diameter
            if self.diameters_a[i] <  200:
                self.seg_a[i].viscosity = estimated_viscosity(self.seg_a[i].diameter) * 1e-3
                self.viscosities_a[i] = self.seg_a[i].viscosity 
            else:
                self.viscosities_a[i] = self.seg_a[i].viscosity 
 
            self.seg_a[i].resistance = (8*self.seg_a[i].viscosity*(self.seg_a[i].length**3)*np.pi)/(self.seg_a[i].volume**2) * (0.00750062) * 1e9
            self.resistances_a[i] = self.seg_a[i].resistance
            
            self.seg_b[i].input_pressure = self.input_pressures_b[i]
            self.seg_b[i].a0 = self.a0_values_b[i]
            self.seg_b[i].volume = self.seg_b[i].a0 * ((self.seg_b[i].input_pressure - PIC)**(1/self.compliance_b))
            self.volumes_b[i] = self.seg_b[i].volume
            self.seg_b[i].diameter =  np.sqrt(((4*self.seg_b[i].volume)/(np.pi*self.seg_b[i].length)))
            self.diameters_b[i] = self.seg_b[i].diameter
            if self.diameters_b[i] <  200:
                self.seg_b[i].viscosity = estimated_viscosity(self.seg_b[i].diameter) * 1e-3
                self.viscosities_b[i] = self.seg_b[i].viscosity 
            else:
                self.viscosities_b[i] = self.seg_b[i].viscosity 
 
            self.seg_b[i].resistance = (8*self.seg_b[i].viscosity*(self.seg_b[i].length**3)*np.pi)/(self.seg_b[i].volume**2) * (0.00750062) * 1e9
            self.resistances_b[i] = self.seg_b[i].resistance


        total_updated_volume_a = np.sum(self.volumes_a) 
        total_updated_volume_b = np.sum(self.volumes_b)
        self.volume = np.array([total_updated_volume_a, total_updated_volume_b])
        self.resistances = np.array([self.resistances_a, self.resistances_b])

        
        total_updated_resistance_a = np.sum(self.resistances_a)
        total_updated_resistance_b = np.sum(self.resistances_b)
        self.resistance = np.array([total_updated_resistance_a, total_updated_resistance_b])
            
        updated_diameter_a = np.sqrt(((4*total_updated_volume_a)/(np.pi*self.length)))
        updated_diameter_b = np.sqrt(((4*total_updated_volume_b)/(np.pi*self.length)))
        self.diameter = np.array([updated_diameter_a, updated_diameter_b])
    
    def scale_pressures(self, remaining_resistance):
        self.scaling_factor_a = np.ones([len(self.input_pressures_a)])
        self.scaling_factor_b = np.ones([len(self.input_pressures_b)])

        for n in range(0, len(self.input_pressures_a)):
#            current_resistance_a = np.sum(self.resistances_a[n:])
#            current_resistance_b = np.sum(self.resistances_b[n:])
#            current_branch_resistance = parallel_resistance(current_resistance_a, current_resistance_b)

#            total_remaining_resistance_a = current_branch_resistance + remaining_resistance[0]
            total_remaining_resistance_a = np.sum(self.resistances_a[n:]) + remaining_resistance[0]
            pressure_drop_a = ((self.input_pressures_a[n] - minimum_pressure)/total_remaining_resistance_a) * self.resistances_a[n]
            perc_pressure_drop_a = pressure_drop_a/self.input_pressures_a[n]
            self.scaling_factor_a[n] = (1-perc_pressure_drop_a)
#            print(self.scaling_factor_a[n])
#            print(self.resistances_a[n])
#            print(total_remaining_resistance_a)

            
#            total_remaining_resistance_b = current_branch_resistance + remaining_resistance[1]
            total_remaining_resistance_b = np.sum(self.resistances_b[n:]) + remaining_resistance[1]
            pressure_drop_b = ((self.input_pressures_b[n] - minimum_pressure)/total_remaining_resistance_b) * self.resistances_b[n]
            perc_pressure_drop_b = pressure_drop_b/self.input_pressures_b[n]
            self.scaling_factor_b[n] = (1-perc_pressure_drop_b)            
            
        for n in range(0, len(self.input_pressures_a)):
            self.input_pressures_a[n] = self.input_pressures_a[n] * self.scaling_factor_a[n]
            self.input_pressures_b[n] = self.input_pressures_b[n] * self.scaling_factor_b[n]
            
#            self.output_pressures_a[n] = self.input_pressures_a[n] * self.scaling_factor_a[n]
#            self.output_pressures_b[n] = self.input_pressures_b[n] * self.scaling_factor_b[n]
            
        output_pressure = (((self.input_pressures_a[-1]*self.volumes_a[-1]) + (self.input_pressures_b[-1]*self.volumes_b[-1]))/(self.volumes_a[-1] + self.volumes_b[-1]))
        self.output_pressure = np.array([output_pressure, output_pressure])
        
    def calculate_flow_branch(self):
        self.flow_values_a = np.array([])
        self.flow_values_b = np.array([])
        self.pressure_dist_a = np.copy(self.input_pressures_a)
        self.pressure_dist_a = np.append(self.pressure_dist_a, self.output_pressure)
        self.pressure_dist_b = np.copy(self.input_pressures_b)
        self.pressure_dist_b = np.append(self.pressure_dist_b, self.output_pressure)
        for t in range(0, self.num_time_steps):
            pressure_gradient_a = self.pressure_dist_a[t] - self.pressure_dist_a[t+1]
            flow_in_segment_a = pressure_gradient_a/self.resistances_a[t] #* 1e-9
            self.flow_values_a = np.append(self.flow_values_a, flow_in_segment_a) 
            
            pressure_gradient_b = self.pressure_dist_b[t] - self.pressure_dist_b[t+1]
            flow_in_segment_b = pressure_gradient_b/self.resistances_b[t] #* 1e-9
            self.flow_values_b = np.append(self.flow_values_b, flow_in_segment_b) 

        self.flow_matrix_a = np.vstack((self.flow_matrix_a, self.flow_values_a))
        self.flow_matrix_b = np.vstack((self.flow_matrix_b, self.flow_values_b))

        #self.flow = np.array([self.flow_values_a[self.num_time_steps//2], self.flow_values_b[self.num_time_steps//2]])
        self.flow = np.array([self.flow_values_a[-1], self.flow_values_b[-1]])


    def update_flows(self):
        """ Update self.flow attribute, use middle segment flow values for vessels a and b. """ 

        self.flows_a = np.copy(self.flows[0])
        self.flow_values_a = np.append(self.flow_values_a, self.flows_a[self.num_time_steps//2])

        self.flows_b = np.copy(self.flows[1])
        self.flow_values_b = np.append(self.flow_values_b, self.flows_b[self.num_time_steps//2])        
        
        self.flow_matrix_a = np.vstack((self.flow_matrix_a, self.flows_a))
        self.flow_matrix_b = np.vstack((self.flow_matrix_b, self.flows_b)) 
        
        #self.flow = np.array([self.flows_a[self.num_time_steps//2], self.flows_b[self.num_time_steps//2]])
        self.flow = np.array([self.flows_a[0], self.flows_b[0]]) #first segment flow


    def calculate_pressures(self):
        self.flows_a = self.flows[0]
        self.flows_b = self.flows[1]
#        for n in range(0, len(self.input_pressures_a)-1):
#            if n == 0:
#                self.input_pressures_a = self.input_pressures_a[0:-1]
#                self.input_pressures_a = np.insert(self.input_pressures_a, 0, new_pressure[0])
#                pressure_gradient_a  = self.flows_a[0] * self.resistances_a[0]
#                P_out_a = self.input_pressures_a[0] - pressure_gradient_a
#                self.input_pressures_a[n+1] = P_out_a 
#            
#                self.input_pressures_b = self.input_pressures_b[0:-1]
#                self.input_pressures_b = np.insert(self.input_pressures_b, 0, new_pressure[1])
#                pressure_gradient_b = self.flows_b[0]*self.flows_b[0]
#                P_out_b = self.input_pressures_b[0] - pressure_gradient_b
#                self.input_pressures_b[n+1] = P_out_b
#                
#            else: 
#                pressure_gradient_a = self.flows_a[n] * self.resistances_a[n]
#                P_out_a = self.input_pressures_a[n] - pressure_gradient_a
#                self.input_pressures_a[n+1] = P_out_a
#                
#                pressure_gradient_b = self.flows_b[n] * self.resistances_b[n]
#                P_out_b = self.input_pressures_b[n] - pressure_gradient_b
#                self.input_pressures_b[n+1] = P_out_b 
#                
#        output_pressure = (((self.input_pressures_a[-1]*self.volumes_a[-1]) + (self.input_pressures_b[-1]*self.volumes_b[-1]))/(self.volumes_a[-1] + self.volumes_b[-1]))
#        self.output_pressure = np.array([output_pressure, output_pressure])
#        self.input_pressure = np.array([self.input_pressures_a[0], self.input_pressures_b[0]])
        
#        for n in range(0, len(self.pressures_a)-1):
#            if n == 0: 
#                self.pressures_a = self.pressures_a[0:-1]
#                self.pressures_a = np.insert(self.pressures_a, 0, new_pressure[0])
#                pressure_gradient_a = self.flows[0][0]*self.resistances[0][0]
#                P_out_a = self.pressures_a[0]-pressure_gradient_a
#                self.pressures_a[n+1] = P_out_a 
#                
#                self.pressures_b = self.pressures_b[0:-1]
#                self.pressures_b = np.insert(self.pressures_b, 0, new_pressure[1])
#                pressure_gradient_b = self.flows[1][0] * self.resistances[1][0]
#                P_out_b = self.pressures_b[0]-pressure_gradient_b
#                self.pressures_b[n+1] = P_out_b 
#                
#            else: 
#                pressure_gradient_a = self.flows[0][n]*self.resistances[0][n]
#                P_out_a = self.pressures_a[n] - pressure_gradient_a
#                self.pressures_a[n+1] = P_out_a 
#                
#                pressure_gradient_b = self.flows[1][n]*self.resistances[1][n]
#                P_out_b = self.pressures_b[n] - pressure_gradient_b 
#                self.pressures_b[n+1] = P_out_b 
#                
#        self.input_pressures_a = np.copy(self.pressures_a)
#        self.input_pressures_a = self.input_pressures_a[0:-1]
#        self.output_pressures_a = np.copy(self.pressures_a)
#        self.output_pressures_a = self.output_pressures_b[1:]
#
#        self.input_pressures_b = np.copy(self.pressures_b)
#        self.input_pressures_b = self.input_pressures_b[0:-1]
#        self.output_pressures_b = np.copy(self.pressures_b)
#        self.output_pressures_b = self.output_pressures_b[1:]
#        
#        self.input_pressure = np.array([self.pressures_a[0], self.pressures_b[0]]) 
#        output_pressure = (((self.pressures_a[-1]*self.volumes_a[-1]) + (self.pressures_b[-1]*self.volumes_b[-1]))/(self.volumes_a[-1] + self.volumes_b[-1]))
#        self.output_pressure = np.array([output_pressure, output_pressure])
#    
#        self.pressures = np.array([self.pressures_a, self.pressures_b])
#
#        
#        # Calculate volumes and resistances of individual segments and find sum  
#        for i in range(0, len(self.input_pressures_a)):
#            self.seg_a[i].input_pressure = self.input_pressures_a[i]
#            self.seg_a[i].a0 = self.a0_values_a[i]
#            self.seg_a[i].volume = self.seg_a[i].a0 * ((self.seg_a[i].input_pressure - PIC)**(1/self.compliance_a))
#            self.volumes_a[i] = self.seg_a[i].volume
#            self.seg_a[i].diameter =  np.sqrt(((4*self.seg_a[i].volume)/(np.pi*self.seg_a[i].length)))
#            self.diameters_a[i] = self.seg_a[i].diameter
##            self.seg_a[i].viscosity = estimated_viscosity(self.seg_a[i].diameter) * 1e-3
#            self.viscosities_a[i] = self.seg_a[i].viscosity 
#            self.seg_a[i].resistance = (8*self.seg_a[i].viscosity*(self.seg_a[i].length**3)*np.pi)/(self.seg_a[i].volume**2) * (0.00750062) * 1e9
#            self.resistances_a[i] = self.seg_a[i].resistance
#            
#            self.seg_b[i].input_pressure = self.input_pressures_b[i]
#            self.seg_b[i].a0 = self.a0_values_b[i]
#            self.seg_b[i].volume = self.seg_b[i].a0 * ((self.seg_b[i].input_pressure - PIC)**(1/self.compliance_b))
#            self.volumes_b[i] = self.seg_b[i].volume
#            self.seg_b[i].diameter =  np.sqrt(((4*self.seg_b[i].volume)/(np.pi*self.seg_b[i].length)))
#            self.diameters_b[i] = self.seg_b[i].diameter
##            self.seg_b[i].viscosity = estimated_viscosity(self.seg_b[i].diameter) * 1e-3
#            self.viscosities_b[i] = self.seg_b[i].viscosity 
#            self.seg_b[i].resistance = (8*self.seg_b[i].viscosity*(self.seg_b[i].length**3)*np.pi)/(self.seg_b[i].volume**2) * (0.00750062) * 1e9
#            self.resistances_b[i] = self.seg_b[i].resistance
#   
#        total_updated_volume_a = np.sum(self.volumes_a) 
#        total_updated_volume_b = np.sum(self.volumes_b)
#        self.volume = np.array([total_updated_volume_a, total_updated_volume_b])
#        self.resistances = np.array([self.resistances_a, self.resistances_b])
#        
#        total_updated_resistance_a = np.sum(self.resistances_a)
#        total_updated_resistance_b = np.sum(self.resistances_b)
#        self.resistance = np.array([total_updated_resistance_a, total_updated_resistance_b])
#            
#        updated_diameter_a = np.sqrt(((4*total_updated_volume_a)/(np.pi*self.length)))
#        updated_diameter_b = np.sqrt(((4*total_updated_volume_b)/(np.pi*self.length)))
#        self.diameter = np.array([updated_diameter_a, updated_diameter_b])

        for n in range(0, len(self.pressures_a)-1):
            pressure_gradient_a = self.flows_a[n] * self.resistances_a[n]
            P_out_a = self.pressures_a[n] - pressure_gradient_a
            self.pressures_a[n] = P_out_a
            
            pressure_gradient_b = self.flows_b[n] * self.resistances_b[n]
            P_out_b = self.pressures_b[n] - pressure_gradient_b
            self.pressures_b[n] = P_out_b      
            
        output_pressure = (((self.pressures_a[-2]*self.volumes_a[-1]) + (self.pressures_b[-2]*self.volumes_b[-1]))/(self.volumes_a[-1] + self.volumes_b[-1]))
        self.output_pressure = np.array([output_pressure, output_pressure])
   
#class SinkSegment(SingleSegment): 
#    """ Class for creating sink segment """
#    def __init__(self, input_pressure, output_pressure, diameter, length, viscosity, compliance = 1, num_time_steps = 1):
#        """ Initialise attributes of parent class """ 
#        super().__init__(input_pressure, output_pressure, diameter, length, viscosity, compliance, num_time_steps)
#    
#        self.volume = np.array([self.volume])
#        self.resistance = np.array([self.resistance])
#        self.input_pressure = np.array([self.input_pressure])
#        self.output_pressure = np.array([self.output_pressure])
#        self.diameter = np.array([self.diameter])
#        
#        self.flow_matrix = np.zeros([1, num_time_steps]) 
#        self.flows = np.zeros([1, self.num_time_steps])
#        
#        self.flow = np.array([0])
#        
#        self.number_of_vessels = 1
#
#    def calculate_a0(self):
#        A_0 = self.volume/(self.input_pressure - PIC)**(1/self.compliance)
#        self.a0 = A_0
#        
#    def change_pressure(self, new_pressure):
#        # if there is a change of pressure in the segment
#        if new_pressure != self.input_pressure:
#            self.input_pressure = np.array([new_pressure])
#            
#            updated_volume = self.a0 * ((self.input_pressure - PIC)**(1/self.compliance))
#            self.volume = np.array([updated_volume]) 
#
#            updated_diameter = np.sqrt(((4*self.volume)/(np.pi*self.length)))
#            self.diameter = updated_diameter  
#            
##            updated_viscosity = estimated_viscosity(self.diameter) * 1e-3
##            self.viscosity = updated_viscosity
#
#            updated_resistance = (8*self.viscosity*(self.length**3)*np.pi)/(self.volume**2)
#            updated_resistance = updated_resistance * (0.00750062) *  1e9
#            self.resistance = np.array([updated_resistance]) 
#            
#    def scale_pressures(self, remaining_resistance):
#        self.output_pressure = np.array([minimum_pressure])
#
#    def calculate_flow_branch(self):
#        pressure_gradient = self.input_pressure - self.output_pressure
#        sink_flow = pressure_gradient/self.resistance #* 1e-9
##        self.sink_flow_matrix = np.vstack((self.sink_flow_matrix, sink_flow))
#        
#        self.flow = np.array([sink_flow])

# Treating last vessel in network as a sink
#class SinkSegment(SingleSegment): 
#    """ Class for creating sink segment - use as last object for Boas replication """
#    def __init__(self, input_pressure, output_pressure, diameter, length, viscosity, compliance, num_time_steps):
#        """ Initialise attributes of parent class """ 
#        super().__init__(input_pressure, output_pressure, diameter, length, viscosity, compliance, num_time_steps)
#    
#        self.input_pressure = np.array([self.input_pressure])
#        self.output_pressure = np.array([self.output_pressure])
#        self.resistance = np.array([self.resistance])
#        self.volume = np.array([self.volume])
#        
#        self.diameter = np.array([self.diameter])
#        
#        self.flow = np.array([0])
#        
#        self.number_of_vessels = 1 
#        
#        # Create segments within long segment 
#        self.seg = [] # list of individual segments within long segment 
#        # create arrays for pressures, volumes and resistances in each individual segment within LongSegment object
#        self.input_pressures = np.array([])
#        self.output_pressures = np.array([])
#        self.volumes = np.array([])
#        self.diameters = np.array([])
#        self.viscosities = np.array([])
#        self.resistances = np.array([])
#        for i in range(0, num_time_steps): 
#            self.seg.append(SingleSegment(input_pressure, output_pressure, diameter, length/num_time_steps, viscosity, compliance, num_time_steps))
#            self.input_pressures = np.append(self.input_pressures, self.seg[i].input_pressure)
#            self.output_pressures = np.append(self.output_pressures, self.seg[i].output_pressure)
#            self.volumes = np.append(self.volumes, self.seg[i].volume)
#            self.diameters = np.append(self.diameters, self.seg[i].diameter)
#            self.viscosities = np.append(self.viscosities, self.seg[i].viscosity)
#            self.resistances = np.append(self.resistances, self.seg[i].resistance)
#     
#        self.sink_flow_matrix = np.zeros([1, num_time_steps]) #flow values for every segment within object 
#        self.sink_pressure_matrix = np.copy(self.input_pressures)
#        self.sink_pressure_matrix = np.append(self.sink_pressure_matrix, self.output_pressure)
# 
#
#    def calculate_a0(self):
#        # Method to calculate a0 constant for every segment
#        self.a0_values = np.array([])
#        for i in range(0, self.num_time_steps):
#            A_0 = self.volumes[i]/(self.input_pressures[i] - PIC)**(1/self.compliance)
#            self.a0_values = np.append(self.a0_values, A_0)
#        
#        A_0 = self.volume/(self.input_pressure - PIC)**(1/self.compliance)
#        self.a0 = A_0
#
#    def change_pressure(self, new_pressure):
#        # If new pressure input is different to current input pressure, update attributes 
#        self.input_pressures = self.input_pressures[0:-1] #update array of input pressures 
#        self.input_pressures = np.insert(self.input_pressures, 0, new_pressure)         
#        
#        if new_pressure != self.input_pressure: 
#            self.input_pressure = np.array([new_pressure]) #update input_pressure for Sink Segment object
#
#        # Update values for individual segments within Sink Segment 
#        # Calculate volumes and resistances of individual segments and find sum  
#        for i in range(0, len(self.input_pressures)):
#            self.seg[i].input_pressure = self.input_pressures[i]
#            self.seg[i].a0 = self.a0_values[i]
#            self.seg[i].volume = self.seg[i].a0 * ((self.seg[i].input_pressure - PIC)**(1/self.compliance))
#            self.volumes[i] = self.seg[i].volume
#            self.seg[i].diameter =  np.sqrt(((4*self.seg[i].volume)/(np.pi*self.seg[i].length)))
#            self.diameters[i] = self.seg[i].diameter
##            self.seg[i].viscosity = estimated_viscosity(self.seg[i].diameter) * 1e-3
#            self.viscosities[i] = self.seg[i].viscosity 
#            self.seg[i].resistance = (8*self.seg[i].viscosity*(self.seg[i].length**3)*np.pi)/(self.seg[i].volume**2) * (0.00750062) * 1e9
#            self.resistances[i] = self.seg[i].resistance
#       
#        total_updated_volume = np.sum(self.volumes)
#        self.volume = np.array([total_updated_volume])
#        
#        total_updated_resistance = np.sum(self.resistances)
#        self.resistance = np.array([total_updated_resistance])
#            
#        updated_diameter = np.sqrt(((4*self.volume)/(np.pi*self.length)))
#        self.diameter = updated_diameter  
#        
#    def calculate_flow_branch(self):
#        self.flow_values = np.array([])
#        self.pressure_distribution = np.copy(self.input_pressures)
#        self.pressure_distribution = np.append(self.pressure_distribution, self.output_pressure)
#        for i in range(0, self.num_time_steps):
#            pressure_gradient = self.pressure_distribution[i] - self.pressure_distribution[i+1]
#            flow_in_segment = pressure_gradient/self.resistances[i] #* 1e-9 #convert to m^3/s
#            self.flow_values = np.append(self.flow_values, flow_in_segment) 
#            
#        self.sink_flow_matrix = np.vstack((self.sink_flow_matrix, self.flow_values))
#        self.sink_pressure_matrix = np.vstack((self.sink_pressure_matrix, self.pressure_distribution))
#        
#        # output flow from middle segment
#        #self.flow = np.array([self.flow_values[self.num_time_steps//2]])
#        self.flow = np.array([self.flow_values[-1]]) #returns flow value for last segment in vessel 
#       
#    def scale_pressures(self, remaining_resistance):
#        self.scaling_factor_ls = np.ones([len(self.input_pressures)])
#        for n in range(0, len(self.input_pressures)):
#            total_remaining_resistance = np.sum(self.resistances[n:]) + remaining_resistance[0]
#            pressure_drop = ((self.input_pressures[n] - minimum_pressure)/total_remaining_resistance) * self.resistances[n]
#            perc_pressure_drop = pressure_drop/self.input_pressures[n]
#            self.scaling_factor_ls[n] = (1-perc_pressure_drop)
##            print(self.scaling_factor_ls[n])
##            print(self.resistances[n])
##            print(total_remaining_resistance)
#        for n in range(0, len(self.input_pressures)):
#            self.input_pressures[n] = self.input_pressures[n] * self.scaling_factor_ls[n]
#        self.output_pressure = np.array([minimum_pressure])
        