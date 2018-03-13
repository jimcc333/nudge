# NUDGE: NUclear Database GEneration software
This software package takes a black-box model and algorithmically selects inputs to run and populate in a database. The created database is intended to be used to train an arbitrary machine learning model. NUDGE tries to select inputs to maximize the information obtained about the underlying blackbox model.

## Running NUDGE
To run nudge, run nudge.py in the nudge/ directory from terminal:

>$ python /src/nudge.py -h

##

 Copyright 2016 Cem Bagdatlioglu

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

## Naming and standards:
The database folder should contain the following. The default names can be changed from the objects.py file.
* basecase.py				xsgen input file containing base-case values
* inputs.txt				file containing database inputs
* /SR_Inputs 				folder containing all screening inputs
* /FR_Inputs 				folder containing all full run inputs
* /SR_Outputs				folder for all screening output libraries
* /FR_Outputs				folder for all full library outputs

## Terms:
- Library number: indicated as [number]. Unique number for input-output pair. Starts at zero.
- Library progress: screening:[0:1), full=1
- Screening library: A library that's run in a short time and that has curtailed outputs
- Varied inputs: Names of inputs that libraries get interpolated on
- Coordinates: the normalized ([0,1]) input array with only the varied inputs, ordered alphabetically
- Voronoi cell: The hyper-dimensional "volume" made by points closest to target point


## Workflow:
1 Start UI and read command line arguments

2 Initialize database

  - If there are inputs, read them; or create the input folders
  
  - Attempt to read the output of an input if available
  
3 Screening
  - Run basecase as screening run
  - Estimate total time, ask if ok to proceed (improve estimate in the background)
  - Monte-Carlo inv-norml dist point sampling
  - Multi-d domain cropping
  - Scout topography map
  
4 Exploration
  - Run basecase, space-filling points
  
5 Exploitation
  - Find highest scored points and inputs near them
  - Run new points
  - Estimate max error
  - Repeat until stop criteria met

## Notes:
- Folder structure and naming chosen to be simple and intuitive. This enables users to copy-paste their existing libraries and easily allow NUDGE to use it in a given database
- xsgen inputs include void and cladding radius, NUDGE also uses thickness in inputs and some workflow
- Creation of new library: 1) Generate input file, 2) Initiate library, 3) Add library object to database
- Constant inputs will be assigned the value in basecase, which will be the 0th library in database
- Dicts in the xsgen input file (initial heavy metal) should be written so that each item is in a new line
- During the Voronoi cell volume calculation, best points to use as inputs during the next-batch are saved too
