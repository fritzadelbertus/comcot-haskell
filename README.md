# Porting Project: COMCOT in Haskell

## Introduction
COMCOT (Cornel Multi-grid Coupled Tsunami Model) is a numerical modeling tool designed for simulating tsunamis, including their generation, propagation, and inundation. The model was initially developed at Cornell University and has been widely used for research and practical applications in tsunami science.

COMCOT is available in [Fortran](https://github.com/AndybnACT/comcot-gfortran), but in this project the COMCOT model will be implemented in Haskell, a Pure-Functional Programming Language. The goal of this project is to help the writer understand and implement concepts related to Functional Programming.

The writer also believes that Pure-Functional Programming has similar syntax writing with how mathematicians write formulas, and so the writer believes writing a numerical model in Haskell must make the code resembles mathematical formula better than Fortran.

## Table of Contents
1. [Introduction](#installation-guide)
2. [How To Use](#how-to-use)
    1. [Usage Recommendation](#usage-recommendation)
    2. [Environment Requirements](#environment-requirements)
    3. [Installation Guide](#installation-guide)
    4. [Running The Program](#running-the-program)
3. [Key Functional Programming Concepts](#key-functional-programming-concepts)
4. [Limitation](#limitation)


## How To Use

Before explaining more on the project, here are a quick steps on how to use the program. The writer assumes the reader have read COMCOT 1.7 Manual and understand how to run haskell program (.hs file extension) in the terminal.

### Usage Recommendation

Based on the [results](#results) the writer recommend the reader to use this project to understand how the COMCOT model is implemented in code, and how Functional Programming implemented the COMCOT model. Thus, it is more useful for academic and theoretic research. For practical research and simulation the writer highly recommend using the [Fortran](https://github.com/AndybnACT/comcot-gfortran) version.

### Environment Requirements
1. GHC version 9.4.8

### Dependencies
Some haskell package dependencies:
1. `Numeric` for `showFFloat`
2. `Data.List` for `transpose`
3. `Data.Maybe` for `fromJust` and `isNothing`

### Installation Guide
1. Setup a directory
2. In terminal, ```git clone https://github.com/fritzadelbertus/comcot-haskell.git```

### Running the Program
1. Please read program [limitation](#limitation) to understand which features are ported and supported in this haskell version.
2. Make sure the required inputs for tsunami modeling is in the directory.
3. Modify the `comcot.ctl`  specificly only modify the `General Parameters for Simulation` part, the `Parameters for Fault Model (Segment 01)` part, and the `Parameters for 1st-level grid -- layer 01` part.
4. Open terminal and write `ghci comcot.hs`
5. The output is available in the `output` directory

## Key Functional Programming Concepts
Below listed concepts in Functional Programming implemented in this project.

1. Pure Functions, most of the program functionality is written in pure function which has no side-effects.
2. Higher-Order Function, `map` and `fold` functions was implemented in list manipulations. Custom higher-order function are also created to achieve point-free style writing.
3. Partial Application, some higher-order fucntions are partialy aplicated to produce a more specific function.
4. Function Composition, other function are alse composed which also result in point-free style writing.
5. Lazy Evaluation, this results in interesting performance result ([See Result]())

## Limitation
This is a 60 hours targeted project across 6 weeks, and the writer admits his incapability to port the full COMCOT features from Fortran to Haskell. Below listed the features which are supported:

### General Parameters
1. Output format is only Timeseries, Z-Max not available
2. Hotstart feature not available
3. Initial Condition only supports Fault Model (Okada) and File (in XYZ format)
4. Boundary Condition only supports Open (Radiation)

### Fault Model Parameters
1. Multi-Fault not available
2. Data format option only support XYZ format

### Wavemaker Parameters
Not available

### Landslide Parameters
Not available

### Layer Parameters
1. Multi-Layer not available
2. Coordinate system only support Cartesian
3. Governing equation only support Linear SWE
4. Bottom Friction and Manning Roughness not available
5. Layer output only support Z values
6. Data Format only support XYZ format

### Other limitation
`fault_multi.ctl` and `landslide.ctl` not available

## Results

## References
Wang, X. (2009). User manual for COMCOT version 1.7.
Andy. COMCOT in Fortran [github](https://github.com/AndybnACT/comcot-gfortran)


