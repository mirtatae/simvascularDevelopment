# Introduction
 
This repository includes a set of MATLAB codes that add different functionalities to the open-source software SimVascular (www.simvascular.org). 

# Adding catheter to the inlet boundary condition

**Language: _MATLAB_**

**Main file: _inletBCT.m_**

### Table of Contents
* [Description](#description)
* [Input Files](#input-files)
  * [How to prepare inlet_coordinate.csv file](#How-to-prepare-inlet_coordinate.csv-file)
* [Example](#example)
* [Reference](#reference)
## Description 
Coming soon after the manuscript (see [Reference](#reference)) is accepted for publication.

## Input Files
The MATLAB function _inletBCT.m_ needs two input .csv file.

- **flowrate.csv**: A .csv file that contains the inlet flow rate (pulsatile). See an [example](https://github.com/mirtatae/simvascularDevelopment/blob/master/example/flowrate.csv).
- **inlet_coordinates.csv**: A .csv file with seven columns that contain the coordinates of the mesh nodes in the inlet plane. See an [example](https://github.com/mirtatae/simvascularDevelopment/blob/master/example/inlet_coordinates.csv). The columns are in the following order:

| GlobalNodeID | VelocityX | VelocityY | VelocityZ | XCoordinate | YCoordinate | ZCoordinate |
| ------------ | --------- | --------- | --------- | ----------- | ----------- | ----------- |
|     ...      |           |           |           |             |             |             |

where VelocityX, VelocityY, and VelocityZ are the velocity components in x, y, z direction (obtained from a steady-state simulation, e.g. mesh independency steady), respectively. XCoordinate, YCoordinate, and ZCoordinate are the (x, y, z) coordinates of the mesh nodes.

### How to prepare inlet_coordinate.csv file

The _inlet_coordinate.csv_ file can be created in different ways. As long as the spreadsheet consists of seven columns as described above, it can be used as the input of the _inletBCT.m_. Here, I describe a way that can be used to create this file. In this method, we run a short simulation and then extract the mesh node information in the open-source software [ParaView](www.paraview.org).

1) Create a simulation job in SimVascular and carry out a short steady-state simulation. Maeke sure that you use the same mesh file that you want to use in the unsteady simulation with catheter. If you have carried out a mesh independency study (to select an appropriate mesh size) that simulation should be enough.
2) Convert the last time step (last restart file).
3) Open the .vtu file in ParaView. Under the "Representation" dropdown menue, select "Surface With Edges". Now, you can see the mesh elements and nodes
4) Use "Select Points On", "Add Selection", and "Subtract Selection" to select all the mesh nodes in the inlet plane. Make sure that **only** the inlet plane mesh nodes are selected.
5) Split the view and add a "SpreadSheetView". Under the "Attribute" dropdown menue, select "Point Data". Click "Show only selected elements". Under the "Toggle column visibility" dropdown menue, select only "GlobalNodeID", "Points", and "Velocity". Now, click "Export Spreadsheet" and save it as _inlet_coordinate.csv_.

## Example
Clone this repository. Set the Matlab directory to this clone. Run the follwoing in the Matlab command line:
```
inletBCT
```

## Reference
Please cite the following manuscript:

Taebi, A., Berk, S., Roncali, E. Realistic boundary conditions in SimVascular through inlet catheter modeling, Under review.
