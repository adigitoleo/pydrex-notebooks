Auxiliary material for Paper 2008GC002032R

Prediction of anisotropy from flow models: A comparison of three methods

Einat Lev and Bradford H. Hager 
Department of Earth, Atmospheric and Planetary Science Massachusetts Institute of
Technology Cambridge, Massachusetts, USA


Introduction 

Auxiliary material for this manuscript includes Figure S1, source code for the program FedRex 
(version 1.1) and source code for the program DirectorsEvolution used in experiments described 
in the text.

Lev, E., and B. H. Hager (2008), Prediction of anisotropy from flow models: 
A comparison of three methods, Geochem. Geophys. Geosyst., 9, Q07014, doi:10.1029/2008GC002032.


1. 2008gc002032-fs01.eps Figure S1. Percentage of \textit{Underworld} execution
time spent on time-integration of the orientation of a set of
directors (diamonds). Circles -- how much of the total time
spent on time-integration (mostly advection of the directors in
space and updating the orientation) was spent on updating the
orientations. We performed this test using a 2D Rayleigh-Taylor
instability model with a grid resolution of 64x64 elements, and
80 (red, green) or 120 (yellow, blue) directors per element.

2. 2008gc002032-FedRex_V1.1.zip Software S1. Source code for the program FedRex used
in the experiments described in the text. FedRex is used
to calculate LPO orientation and magnitude for olivine (or
olivine+enstitite) aggregates advected through a given velocity
field. FedRex is a modification of D-Rex (Kaminski and Ribe
2004), and includes several new and useful features. This ZIP
file containing Fortran90 (*.f90) files, a README file
(FedRex.readme), and input examples (input.dat* and *.00001.dat
velocity field files). This code was compiled and tested only
on a Linux machine using Intel Fortran90 (ifort) compiler.

3) 2008gc002032-DirectorsEvolution.zip Software S2. Source code for the program
DirectorsEvolution used in experiments described in the
text. This version of the code is essentially 2D and
steady-state only. A time-dependent three-dimensional version
can be obtained from the authors. This program mimics the
orientation evolution of directors as done in Underworld
(following M\"{u}hlhaus et al., 2004). This is a ZIP archive
containing three Matlab files (*.m), and three input examples
(*.dat).
