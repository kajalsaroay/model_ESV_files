# Thesis title
A Dynamic Computational Model of Pulsatile Brain Blood Flow 

# Purpose 
Develop a model to look at effects of compliance and pulsatility on blood flow in a network of cerebral blood vessels. 

# How to Use 
## List of required files: 
* vessels.py - contains all classes and methods used to create input network of vessels. 
* model_ESV.py - creates storage matrices, functions to map between level and matrix positions, initial (t=0) values, runs main loop, saves output values as text files.
* pressure_wave.txt - example of an input pressure wave.  
* ESV_exec.py - Use with model_ESV_function.py for required input files. Change file paths for input file directory and compliances for chosen files. 
* network text file e.g. input_network.txt - Holds information about number of vessels per level, branch pattern for levels, type of vessels, all baseline parameters. 

*vessels.LongSegment(Input pressure, Output pressure, Diameter, Length, Blood Viscosity, Compliance (beta value), Number of Segments)* 

*vessels.DivergingBranchSegment(Input pressure, Output pressure, np.array([Diameter vessel a, Diameter vessel b]), Length, Blood Viscosity, np.array([Compliance vessel a, Compliance vessel b]), Number of Segments)* 

## Other files:
* network_input_replace.py - changes compliance values in generic text files.  
* plausible_vessel_network_dev.py - used for calculating parameters for plausible vessel network. Input chosen diameter and blood velocity values for every level, output correct lengths for every level. Calculates parameters for three compartments (arteries, microvessels and veins). 


### model_ESV.py (iterating through time steps)
1. Change pressure - New pressure in every vessel (determined by the pressure wave and how far it has travelled across the network), update volume and resistance values at every time point. If in level 0, use pressure wave as new pressure, if in every other level, new pressure comes from previous vessel at previous time point (output_pressure_matrix). 
2. ESV calculation - Reduce network to Equivalent Single Vessel (ESV). Calculate total value of resistance for every segment, calculate pressure drop across every segment, calculate flow through every segment. Distribute flow in corresponding segments in main network. 
3. Update object's flow attribute.
4. Calculate change in pressure across segments using flow and resistance values - Update output_pressure_matrix. 
