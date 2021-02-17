# Introduction
 
This repository includes a set of MATLAB codes that add different functionalities to the open-source software SimVascular (www.simvascular.org). 

# Adding catheter to the inlet boundary condition
~~~
Language: MATLAB

Main file: inletBCT.m
~~~
### Table of Contents
* [Description](#description)
* [Input Files](#input-files)
* [Example](#example)
* [Reference](#reference)
## Description 
comming soon after the manuscript is accepted for publication.

## Input Files
The MATLAB function inletBCT.m needs two input .csv file.
1) flowrate.csv: A .csv file that contains the inlet flow rate (pulsatile).
2) inlet_coordinates.csv: A .csv file with seven columns that contain the coordinates of the mesh nodes in the inlet plane. The columns are in the following order:
GlobalNodeID | VelocityX | VelocityY | VelocityZ | XCoordinate | YCoordinate | ZCoordinate
where VelocityX, VelocityY, and VelocityZ are the velocity components in x, y, z direction (obtained from a steady-state simulation, e.g. mesh independency steady), respectively. XCoordinate, YCoordinate, and ZCoordinate are the (x, y, z) coordinates of the mesh nodes.

### How to prepare inlet_coordinate.csv file

## Example
Description coming soon ...

## Reference
Please cite the following manuscript:

Taebi, A., Berk, S., Roncali, E. Realistic boundary conditions in SimVascular through inlet catheter modeling, Under review.
