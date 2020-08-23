# EEGSimulationStudy
Source Code for Paper:
Imperatori, L. S., Betta, M., Cecchetti, L., Canales-Johnson, A., Ricciardi, E., Siclari, F., ... &amp; Bernardi, G. (2019). EEG functional connectivity metrics wPLI and wSMI account for distinct types of brain functional interactions. Scientific reports, 9(1), 1-15. https://www.nature.com/articles/s41598-019-45289-7

The code is fundamentally based on the Berlin Brain Connectivity Framework:
Haufe, S., & Ewald, A. (2019). A simulation framework for benchmarking EEG-based brain connectivity estimation methodologies. Brain topography, 32(4), 625-642. https://link.springer.com/article/10.1007/s10548-016-0498-y

It also refers to the Chaotic Systems Toolbox (MATLAB) and the PALM (Permutation Analysis of Linear Models) software: 
https://it.mathworks.com/matlabcentral/fileexchange/1597-chaotic-systems-toolbox
https://www.sciencedirect.com/science/article/pii/S1053811914000913.

To replicate the plots from the paper, first generate simulated EEG datasets with different underlying source dynamics and compute wPLI and wSMI connectivity metrics using the main_script.m in the "GenerateEEGDatasets" folder, followed by computing topographic and whole-brain accuracy in the "TestPerformance" folder.

Note that it is a very computationally expensive process that requires access to computational clusters.
