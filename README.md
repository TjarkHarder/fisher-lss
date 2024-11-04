# README

This code accompanies my master's thesis (see the '/doc' directory) in which I investigate non-Gaussian corrections to the estimators of the (one-loop) power spectrum and (tree-level) bispectrum in a Fisher forecast for a Euclid-like survey.


## Example

For anyone interested enough in using this code, the 'example.c' file in the '/example' directory shows the basic usage in order to produce the results in my master's thesis.


## Compile

In order to use the code, one should run


## Requirements

- OpenMP
- GNU Scientific Library (GSL) (see https://www.gnu.org/software/gsl/)
- Basic Linear Algebra Subsystems (BLAS)
- Linear Algebra PACKage (LAPACK)

Additionally, SPLINTER by Bjarne Grimstad (see https://github.com/bgrimstad/splinter) is used for several interpolation functions, as well as Cuba by Thomas Hahn (see https://feynarts.de/cuba/) is used for several bin-average integrals. These libraries are already included in the '/external' directory and are also compiled by following the instructions outlined above.

One important note regarding Cuba, I have slightly modified the source code (line 72 in Cuba-4.2.2/src/common/Fork.c) to forcefully disable parallisation in Cuba, as it seems to mess with OpenMp.
