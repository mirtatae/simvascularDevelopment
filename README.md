# Introduction
 
This repository includes a set of MATLAB codes that add different functionalities to the open-source software SimVascular (www.simvascular.org). 

# Adding catheter to the inlet boundary condition

### Table of Contents
* [Description](#description)
* [Input Files](#input-files)
  * [How to prepare inlet_coordinate.csv file](#how-to-prepare-inlet_coordinatecsv-file)
* [Example](#example)
* [Reference](#reference)
## Description 
This pipeline maps a volumetric flow rate of interest to the inlet plane, which contains a catheter with an adjustable radius, wall thickness, and eccentricity. The pipeline's output is a boundary condition file (_bct.dat_) compatible with the software SimVascular.

More description coming soon after the manuscript (see [Reference](#reference)) is accepted for publication.

## Input Files
The MATLAB function _inletBCT.m_ needs two input .csv file.

- **flowrate.csv**: A .csv file that contains the inlet flow rate (pulsatile). See an [example](https://github.com/mirtatae/simvascularDevelopment/blob/master/example/flowrate.csv).
- **inlet_coordinates.csv**: A .csv file with seven columns that contain the coordinates of the mesh nodes in the inlet plane. See an [example](https://github.com/mirtatae/simvascularDevelopment/blob/master/example/inlet_coordinates.csv). The columns are in the following order:

| GlobalNodeID | VelocityX | VelocityY | VelocityZ | XCoordinate | YCoordinate | ZCoordinate |
| ------------ | --------- | --------- | --------- | ----------- | ----------- | ----------- |
|     ...      |           |           |           |             |             |             |

where VelocityX, VelocityY, and VelocityZ are the velocity components in x, y, z direction (obtained from a steady-state simulation, e.g. mesh independency steady), respectively. XCoordinate, YCoordinate, and ZCoordinate are the (x, y, z) coordinates of the mesh nodes. The columns corresponding to _velocity_ are used to identify the mesh nodes corresponding to the arterial wall (zero velocity).

### How to prepare inlet_coordinate.csv file

The _inlet_coordinate.csv_ file can be created in different ways. As long as the spreadsheet consists of seven columns as described above, it can be used as the input of the _inletBCT.m_. Here, I describe a way that can be used to create this file. In this method, we run a short simulation and then extract the mesh node information in the open-source software [ParaView](www.paraview.org).

1) Create a simulation job in SimVascular and carry out a short steady-state simulation. Maeke sure that you use the same mesh file that you want to use in the unsteady simulation with catheter. If you have carried out a mesh independency study (to select an appropriate mesh size) that simulation should be enough.
2) Convert the last time step (last restart file).
3) Open the .vtu file in ParaView. Under the "Representation" dropdown menue, select "Surface With Edges". Now, you can see the mesh elements and nodes
4) Use "Select Points On", "Add Selection", and "Subtract Selection" to select all the mesh nodes in the inlet plane. Make sure that **only** the inlet plane mesh nodes are selected.
5) Split the view and add a "SpreadSheetView". Under the "Attribute" dropdown menue, select "Point Data". Click "Show only selected elements". Under the "Toggle column visibility" dropdown menue, select only "GlobalNodeID", "Points", and "Velocity". Now, click "Export Spreadsheet" and save it as _inlet_coordinate.csv_.

## Example
Clone this repository. Set the MATLAB directory to this clone. Adjust the values of the parameters in the "switches" and "parameter definition" sections of the _inletBCT.m_. Set appropriate values for _directory_, _filename_ (corresponding to the _flowrate.csv_), and _filename2_ (corresponding to the _inlet_coordinate.csv_).
```ruby
directory = 'example\';
filename = 'flowrate.csv';
filename2 = 'inlet_coordinates.csv';
```

Run the follwoing in the MATLAB command line:
```
inletBCT
```
This will generate a _bct.dat_ file. After you "Create Data Files for Simulation" in SimVascular, you need to replace the _bct.dat_ file that is created in the simulation directory with the _bct.dat_ file that is generated using this MATLAB code. Now, you are ready to run your CFD simulation with a catheter in the inlet plane.

## Reference
Please cite the following manuscript:

Taebi, A., Berk, S., Roncali, E. Realistic boundary conditions in SimVascular through inlet catheter modeling, Under review.

## Developed In

[Roncali Lab](https://roncalilab.engineering.ucdavis.edu/) at @ucdavis

| <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/f/f3/The_University_of_California_Davis.svg/500px-The_University_of_California_Davis.svg.png" width="100"> | <img src="https://uploads-ssl.webflow.com/5f71f6ba15ef4216be8dd209/5f7619583a504af1f2b64115_logo-p-500.png" width="100"> |
|------------|-------------|
