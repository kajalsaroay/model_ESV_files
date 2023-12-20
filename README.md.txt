#Project/Thesis title
A multi-scale computational model of brain blood flow: enhancing MRI for assessing brain health.

#Purpose  
Develop a model to look at effects of compliance and pulsatility on blood flow in a network of cerebral blood vessels.

# Project Description 
Input network describing chosen vessel parameters and chosen pressure wave.  
Output simulated pressure and flow distribution values for chosen input network. 

# How to Use 
List of required files: 
* vessels.py - contains all classes and methods used in Boas_replication_ESV.py 
* Boas_replication_ESV.py - creates storage matrices, functions to map between level and matrix positions, initial (t=0) values, runs main loop, saves output values as text files 
* ESV_exec.py - runs Boas_replication_ESV.py for required input files. Change file paths for input file directory and compliances for chosen files. 
* network text file e.g. Boas_network_100000.txt - Holds information about number of vessels per level, branch pattern for levels, type of vessels, all baseline parameters. 
_vessels.LongSegment(Input pressure, Output pressure, Diameter, Length, Blood Viscosity, Compliance (beta value), Number of Segments)_
_vessels.DivergingBranchSegment(Input pressure, Output pressure, np.array([Diameter vessel a, Diameter vessel b]), Length, Blood Viscosity, np.array([Compliance vessel a, Compliance vessel b]), Number of Segments)_

Other files:
* network_input_replace.py - changes compliance values in generic text files 
* calculate_vessel_lengths.py - used for calculating parameters for realistic network. Input chosen diameter and blood velocity values for every level, output correct lengths for every level. 
* python files to create figures 

### Boas_replication_ESV.py (main loop) 
1. Change pressure - New pressure in every vessel, update volume and resistance values. If in level 0, use pressure wave as new pressure, if in every other level, new pressure comes from previous vessel at previous time point (output_pressure_matrix). 
2. ESV calculation - Reduce network to Equivalent Single Vessel (ESV). Calculate total value of resistance for every segment, calculate pressure drop across every segment, calculate flow through every segment. Distribute flow in corresponding segments in main network. 
3. Update object's flow attribute
4. Calculate change in pressure across segments using flow and resistance values - Update output_pressure_matrix 

