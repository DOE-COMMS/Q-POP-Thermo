![GitHub release version](https://img.shields.io/github/v/release/DOE-COMMS/Q-POP-Thermo?color=%2350C878&include_prereleases)
![License](https://img.shields.io/github/license/DOE-COMMS/Q-POP-Thermo)
![GitHub Size](https://img.shields.io/github/repo-size/DOE-COMMS/Q-POP-Thermo)
![HitCount](https://hits.dwyl.com/DOE-COMMS/Q-POP-Thermo.svg?style=flat-square&show=unique)
![HitCount](https://img.shields.io/endpoint?url=https%3A%2F%2Fhits.dwyl.com%2FDOE-COMMS%2FQ-POP-Thermo.json&label=total%20hits&color=pink)

# Q-POP-Thermo: A Thermodynamic Solver for Ferroic and Quantum Materials

## Overview

Q-POP-Thermo was developed as part of U.S. Department of Energy Award No. DE-SC0020145 for solving the thermodynamic approximations of ferroelectric systems, both in bulk and thin film boundary conditions. The software allows researchers not fluent in programming to generate accurate approximations of ferroelectric properties and domain structures. The long-term goal of Q-POP-Thermo is to provide speedy and accurate approximations of material systems and their microstructures to reduce the number of phase-field simulations to eludicate material phase diagrams. 

At current, the software calculate monodomain and polydomain approximations of ferroelectric thin films and bulk single crystals under various boundary conditions and external fields (i.e. thermal, mechanically clamped, mechanically stressed, electric fields, etc.). Additional systems, such as metal-insulator transitions (MIT, e.g. Vanadium Oxide (VO<sub>2</sub>)), structural transitions, electronic transitions, and ferroelastics are in development. 

## Acknowlegments
# <a name="Thanks"></a>

This development of the software was supported as part of the Computational Materials Sciences Program funded by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, under Award No. DE-SC0020145.

1. [Motativation](#Motivate)
2. [Using Q-POP-Thermo](#Using)
2. [Examples](#Examples)
3. [Future Works](#Future)
3. [Citing Q-POP-Thermo](#Citing)
4. [Licensing](#License)
5. [Library and Module Information](#LibraryInfo)

## Motivation 
# <a name="Motivate"></a>

Currently, there exists a great deal of various thermodynamic solvers for calculating the properties of materials. For example, there exist multiple programs for implementing the CALPHAD model for deriving phase diagrams and properties of common material systems. However, no such program exists for Ferroic materials. Q-POP-Thermo seeks to change that. Providing researchers with an adaptable and expendable framework, in Python, to calculate the properties and stable states of ferroelectric thin films and single crystal systems. For more information on Q-POP-Thermo and the implementation of the necessary equations and boundary conditions, please see the provided manuscript. For citation information please see [Citing Q-POP-Thermo](#Citing).

## Using Q-POP-Thermo
# <a name="Using"></a>

The code is available for researchers and individuals interested in extending the code to other systems or boundary conditions. In addition, the program can be executed as is to study different ferroelectric material systems by calling Q-POP-Thermo from the command line like any other python program, e.g.

`$ python3 qpopthermo.py`

## Examples
# <a name="Examples"></a>

In the provided white paper, we present 5 separate case studies to prove the validity of the Q-POP-Thermo program. Each of these case studies can be viewed in the *Examples* directory located in this repo.

## Future of Q-POP-Thermo
# <a name="Future"></a>

There are a variety of additions that we intend to make to Q-POP-Thermo in the future. The following Development Plan outlines the proposed additions.

| Version Number | Proposed Date | Software Addition(s) |
| ------ | ----- | ------ |
| 1.1  | Winter 2022 | Metal-Insulator Transitions, Multi-layer Architectures, Multi-Domain Calculations |
| 1.2  | Spring 2023   | Energy Surfaces, Minimum Switching Pathways |
| 1.3  | Summer 2023  | Chemical Thermodynamics |
| 1.4  | Fall 2023  | High-Throughput Preprocessing, High-Throughput Post-processing, Parallelization of Solving Routines |

## Citing Q-POP-Thermo
# <a name="Citing"></a>

If you utilize Q-POP-Thermo in your research, please cite the following manuscript: 

Zorn, J.A., et al. "Q-POP-Thermo: A general-purpose thermodynamics solver for ferroelectric materials" *Computer Physics Communications* **[275]** (2022) 108302

## Licensing
# <a name="License"></a>

This software is licensed under the MIT License.

## Library and Module Information
# <a name="LibraryInfo"></a>

The following versions are necessary to be installed in order for Q-POP-Thermo to work properly.

| Library | Version |
| ---- | ----- |
| Numpy | 1.20.3 |
| Pandas | 1.1.3 |
| Sympy | 1.8 |
| Scipy | 1.6.3 |

