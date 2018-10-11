# README #

Alphahouse contains most of the multipurpose functions required for alphasuite software.

It contains three folders, 


BioComputational - for procedures and functionality that is tied to biology
Utilities - for procedures that aren't tied to biology in any way
tests - unit tests written in pfunit

There is a standard `main.f90` in the main src/ folder. This can effectively be used to prototype code from alphahouse.

If you build the cmake file in this directory, an alphahouse executable will be built with the `main.f90` file, otherwise it will be left out. 

