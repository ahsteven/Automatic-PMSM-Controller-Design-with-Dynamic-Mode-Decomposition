
### Identifying and compensationg PMSM harmonics with Dynamic Mode Decomposition

This code is supplementary material for the paper:
#### "Automated Predictive Current Harmonic Compensation on the PMSM using Dynamic Mode Decomposition"
#### authors Adam Stevens and Iqbal Husain

There are two main files to run the examples:
#### Harmonic_Identification_DMD_vs_FFT.m
Reads in actual data from a PMSM for harmonic identification.
The algorithm compares DMD to the FFT.
The file calls the dmd function to run the DMD analysis.
The DMD script is from Scott Dawson.

For more information, see Dawson, Hemati, Williams & Rowley,
"Characterizing and correcting for the effect of sensor noise in the dynamic mode decomposition", 
Experiments in Fluids, 2016.

The other main file is :
#### run_harmonic_compensation.m
This runs a simulink simulation which contains the DMD harmonic compensation algorithm.
It calls the function:
##### get_comp_vecs_algo_multi.m
This function generates a table of DMD vectors for each operating speed of the motor

##### ASHE_DEMO.m 
shows an implementation of the ASHE algorithm that was used as a harmonic compensation benchmark.
For more information see:

Blasko, Vladimir. "A novel method for selective harmonic elimination 
in power electronic equipment." IEEE transactions on Power 
Electronics 22, no. 1 (2007): 223-228.

