#!/usr/bin/python3
#
#  main.py
#  
#  Copyright 2016 cem <cem@cem-VirtualB>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  Naming and standards:
#	The database folder should contain:
#		- /SR_Inputs 				folder containing all scouting inputs
#			- [number].py			input file for scout library [number]
#		- /FR_Inputs 				folder containing all full run inputs
#			- [number].py			input file for full library [number]
#		- /SR_Outputs				folder for all scouting output libraries
#			- /build-[number] 		this number must match the one in SR_Inputs
#				- /brightlite0		created by xsgen
#					- [nucid].txt	output for [nucid] nuclide of input [number]
#		- /FR_Outputs				folder for all full library outputs
#			- /build-[number] 		this number must match the one in FR_Inputs
#				- /brightlite0		created by xsgen
#					- [nucid].txt	output for [nucid] nuclide of input [number]
#		- basecase.py				xsgen input file containing base-case values
#		- inputs.txt				file containing database inputs
#  
#
#	Terms:
#		- Library number: indicated as [number]. Unique number for input-output pair. Starts at zero.
#		- Library progress: scout:[0:1), full=1
#		- Scout library: A library thats run in a short time and that has curtailed outputs
#		- Metric: Names of inputs that libraries get interpolated on
#		- Coordinates: the normalized ([0,1]) metrics with only the varied ones so that dbase dimensions match coordinate dimensions
#		- Neighborhood: Determined by inputs, the "closest" libs to a given lib (for gradient estimation)
#
#	Flags:
#		-g (guided): start NUDGE in guided mode 
#		-d (database): used for the database path
#		-h (help):  help screen
#
#
#	Notes:
#		- Folder structure and naming chosen to be simple and intuitive. This enables users to copy-paste
#		  their existing libraries and easily allow NUDGE to use it in a given database
#		- xsgen inputs include void and cladding radius, NUDGE also uses thickness in inputs and some workflow
#		- Dicts in the xsgen input file (initial heavy metal) should be written so that each item is in a new line
#
#
#
#TODO: updatemetrics fails when a variable isnt varied by a divide by zero error

from objects import *
from interface import *
import os

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main(args):
	os.system('cls' if os.name == 'nt' else 'clear')
	
	screen = Screen()
	
	# Check if help is requested
	if '-h' in args:
		screen.HelpScreen()
		return
	
	# Initialize screen
	screen.InitScreen()
	
	# Take user inputs
	
	# Initialize paths
	try:
		paths = PathNaming(database_path = args[args.index('-d')+1])
	except ValueError:
		usr_path = input('No database path found, please enter full path to database: \n')
		paths = PathNaming(database_path = usr_path)
	# Check xsgen run command
	try:
		paths.xsgen_command = args[args.index('-x')+1]
	except ValueError:
		pass
	
	# Initiate database
	database = DBase(paths, 2)	
	database.UpdateData()
	database.UpdateMetrics()
	screen.UpdateInfo(database)
	
	#os.system(str(paths.xsgen_command))
	
	#database.UpdateNeigbors()
	#database.UpdateNeigbors(slib=0)
	
	#database.PCA()
	
	#database.EstLib(database.slibs[8:10], 	t_enrichment=0.01, t_cool_density=0.99)
	
		
	

	"""
	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111, projection='3d')
	ax1.scatter(database.data_mat[:,0], database.data_mat[:,1], database.data_mat[:,2])
	ax1.set_xlabel("Neutron Production")
	ax1.set_ylabel("Neutron Destruction")
	ax1.set_zlabel("Burnup")
	
	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111, projection='3d')
	ax2.scatter(database.pca_mat.Y[:,0], database.pca_mat.Y[:,1], database.pca_mat.Y[:,2])
	ax2.set_xlabel("PC1")
	ax2.set_ylabel("PC2")
	ax2.set_zlabel("PC3")
	
	fig3 = plt.figure()
	ax3 = fig3.add_subplot(111, projection='3d')
	# Size by PC1 value
	print(database.pca_mat.Y[:,0])
	s = (database.pca_mat.Y[:,0]+2) **2 * 10	
	print(s)
	ax3.scatter(database.fuel_cell_radius, database.enrichment, database.clad_cell_thickness, s=s)
	ax3.set_xlabel("Fuel Radius")
	ax3.set_ylabel("Enrichment")
	ax3.set_zlabel("Clad Thickness")
	plt.show()
	"""
	
	input('')
	screen.PrintAt(colors.reset,y=screen.lines)
		
	
	return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
