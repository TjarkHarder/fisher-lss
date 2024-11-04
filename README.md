# README

I have developed this code ('fisher-lss') as part of my Master's Thesis at Heidelberg University under supervision by Luca Amendola at the Institute of Theoretical Physics (ITP).

The code calculates the Fisher information matrix of model independent variables of the one-loop Power Spectrum and tree-level Bispectrum for a Euclid like survey (or any other survey, given the correct fiducials).
To extend the usual Gaussian approximation of the PP, BB and PB covariance matrices, the code includes the next higher order terms, allowing for correlation between different power spectrum and bispectrum modes in the respective covariance matrix.

The code is mostly finished and running, but needs testing for the B only and P+B Fisher matrices and their forecasts. It also requires some optimisation in some areas, which I am currently working on.


## Example

The 'example.c' file in the '/example' directory has been used to produce the preliminary results in the '/output/spec' and '/output/fish' directories. To get the contraints on the variables (i.e. the '...-err.dat' files) + the triangle plots in the '/plots' directory, one can use the 'fisher-lss-plots.py' script.


## Requirements
The code requires the GNU Scientific Library (GSL) and SPLINTER by Bjarne Grimstad (https://github.com/bgrimstad/splinter) to run. To be able to run 'make fisher-lss', one should make sure the libraries for GSL and SPLINTER are properly linked in 'Makefile'. The usage of the created static library is shown in '/example'.


## Documentation
In progress...
