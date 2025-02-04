This folder contains all the data and scripts used for Fig.3

MieStripes.fsp - Ansys Lumerical FDTD file with defined geometry and sweeps
MieStripes_RefX/RefY/X/Y - Folders with Ansys Lumerical FDTD files containing data from finished simulations, and videos
MieStripes_RefX/RefY/X/Y.mat - .mat files containing simulation data (for loading to Matlab)
BLS_signal2D.m - Matlab script for BLS signal calculation using the reciprocity theorem
DataLoader.m - Matlab script for easier data loading from .mat files (used in BLS_signal2D.m)
SpinWaveGreen.m - Matlab script that calculates the density of states for BLS signal calculations (used in BLS_signal2D.m)
SpinWaveToolkit.py - python module for calculating spin wave characteristics (used in BLS_signal2D.m)
inferno.m - colormap for Matlab (used in BLS_signal2D.m)
Results.mat - .mat file containing results from the script BLS_signal2D.m
Figure_1.png - calculated BLS signal for A = 150 nm 
(we used to describe the periodicity with half period A/2 = a = 75 nm)
Figure_2.png - calculated dynamic susceptibility tensor
Figure_3.png - calculated transfer function
