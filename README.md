# ShockWaves

This is a collection of functions that implements various equations for shock wave properties. 
This package is written primarily as a tool for calculating shock wave properties in the Julia REPL.
Currently, only properties for normal shock waves and moving shock waves are implemented, though support for oblique shock waves would be nice to have.

## Usage

Functions in this package generally have their first argument be a struct of a certain type.
This is done to allow dispatch to the correct function based on the type of shock wave you are wanting to calculate.
So, for example, to calculate the density ratio across a normal shock wave with a Mach number of 1.2 and a ratio of specific heats of 1.4, you would use
`DensityRatio( NormalShockWave(), 1.4, 1.2)`
or, for a moving shock wave with a ratio of specific heats of 1.4 and a pressure ratio of 1.1
`DensityRatio( MovingShockWave(), 1.4, 1.1)`
Alternatively, if you know the Mach number of your moving shock wave but not the pressure ratio, you can use
`DensityRatioFromMachNumber( MovingShockWave(), 1.4, 1.2 )`
where 1.4 is again the ratio of specific heats and 1.2 is the Mach number.

As always, the built-in Julia help and tab completion can be used to discover what functions are available:
```
using ShockWaves

?ShockWaves
```
and information on the functions, including arguments, source equations, and other implementation information can be found similarly:
```
?DensityRatio
```

## Notes


All equations come from Anderson, John D. "Modern Compressible Flow with Historical Perspective" 3rd ed. (2003) McGraw-Hill. ISBN: O-01-242443-5

Specifically, relationships come from
* Normal Shocks: Chapter 3
* Moving Shocks: Chapter 7

While I have done all I can to set up proper tests and ensure these functions are implemented properly, mistakes can (and do) happen. If you find a function that is not behaving as expected, feel free to submit a pull request.
