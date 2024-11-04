# README

This code accompanies my master's thesis (see the '/doc' directory) in which I investigate non-Gaussian corrections to the estimators of the (one-loop) power spectrum and (tree-level) bispectrum in a Fisher forecast for a Euclid-like survey.


## Requirements

For anyone interested enough in using this code, one should first make sure that the following requirements are met

- OpenMP
- GNU Scientific Library (GSL) (see https://www.gnu.org/software/gsl/)
- Basic Linear Algebra Subsystems (BLAS)
- Linear Algebra PACKage (LAPACK)
- cmake, g++, gcc, ... (build_essential)

Additionally, SPLINTER by Bjarne Grimstad (see https://github.com/bgrimstad/splinter) is used for several interpolation functions, as well as Cuba by Thomas Hahn (see https://feynarts.de/cuba/) for several bin-average integrals. These libraries are already included in the '/external' directory and are also compiled by following the instructions outlined below.

One important note regarding Cuba, I have slightly modified the source code (line 72 in Cuba-4.2.2/src/common/Fork.c) to forcefully disable parallisation in Cuba, as it seems to mess with OpenMp.


## Compile

In order to compile the code into a static library, one may use the command

- make fisher-lss

To showcase the basic useage of the library, the 'example.c' file in the '/example' directory shows the majority of the functions needed to produce the (numerical) results in my master's thesis. TO compile the example, one may simply call

- make example

or also

- make all

if one wishes to compile the example automatically after the library. In order to undo previous make calls, one can simply do so by calling

- make clean
