# README

This code accompanies my master's thesis (see the '/doc' directory) in which I investigate non-Gaussian corrections to the estimators of the (one-loop) power spectrum and (tree-level) bispectrum in a Fisher forecast for a Euclid-like survey.


## Requirements

- OpenMP
- GNU Scientific Library (GSL) (see https://www.gnu.org/software/gsl/)
- Basic Linear Algebra Subsystems (BLAS)
- Linear Algebra PACKage (LAPACK)
- cmake, g++, gcc, ... (build_essential)

Additionally, SPLINTER by Bjarne Grimstad (see https://github.com/bgrimstad/splinter) is used for several interpolation functions, as well as Cuba by Thomas Hahn (see https://feynarts.de/cuba/) is used for several bin-average integrals. These libraries are already included in the '/external' directory and are also compiled by following the instructions outlined above.

One important note regarding Cuba, I have slightly modified the source code (line 72 in Cuba-4.2.2/src/common/Fork.c) to forcefully disable parallisation in Cuba, as it seems to mess with OpenMp.


## Compile

For anyone interested enough in using this code, the 'example.c' file in the '/example' directory shows the basic usage in order to produce the results in my master's thesis.

For this, one may use the provided Makefile to create a static library and compile the example by running

- make all

If one only wishes to create the library, or only the example, run

- make fisher-lss
- make example

In order to undo the previous step, simply run

- make clean
