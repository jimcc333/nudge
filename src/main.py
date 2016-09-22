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
#		- Library progress: screening:[0:1), full=1
#		- Scout library: A library thats run in a short time and that has curtailed outputs
#		- Metric: Names of inputs that libraries get interpolated on
#		- Coordinates: the normalized ([0,1]) metrics with only the varied ones so that dbase dimensions match coordinate dimensions
#		- Neighborhood: Determined by inputs, the "closest" libs to a given lib (for gradient estimation)
#		- Voronoi cell:
#
#
#	Workflow:
#		1 Start UI and read command line arguments
#		2 Initialize database
#			- If there are inputs, read them; or create the input folders
#				- Attempt to read the output of an input if available
#		3 Screening
#			- Run basecase as screening run
#			- Estimate total time, ask if ok to proceed (improve estimate in the background)
#			- Monte-Carlo inv-norml dist point sampling
#				- Multi-d domain cropping
#				- Scout topography map
#		4 Exploration
#			- Run basecase, space-filling points
#		5 Exploitation
#			- Find highest scored points and inputs near them
#			- Run new points
#			- Estimate max error
#			- Repeat until stop criteria met
#
#	Flags:
#		-m (manual): start NUDGE in manual mode 
#		-d (database): used for the database path
#		-h (help):  help screen
#		-x (xsgen): command to run xsgen
#
#
#	Notes:
#		- Folder structure and naming chosen to be simple and intuitive. This enables users to copy-paste
#		  their existing libraries and easily allow NUDGE to use it in a given database
#		- xsgen inputs include void and cladding radius, NUDGE also uses thickness in inputs and some workflow
#		- Dicts in the xsgen input file (initial heavy metal) should be written so that each item is in a new line
#		- During the Voronoi cell volume calculation, best points to use as inputs during the next-batch are saved too
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
	
	# Manual mode check
	if '-m' not in args:
		# Take user inputs
		# 	Initialize paths
		try:
			paths = PathNaming(database_path = args[args.index('-d')+1])
		except ValueError:
			usr_path = input('No database path found, please enter full path to database: \n')
			paths = PathNaming(database_path = usr_path)
		# 	Check xsgen run command
		try:
			paths.xsgen_command = args[args.index('-x')+1]
		except ValueError:
			pass
		
		# Initiate database
		database = DBase(paths)		# Read all available inputs and outputs in the folder
		database.UpdateMetrics()	# Update database data about the library inputs, outputs, and states
		screen.UpdateInfo(database)	# Print new info on screen
		
		database.Print()
	else:
		# Manual mode
		usr_path = 'db2/'
		paths = PathNaming(database_path = usr_path)
		database = DBase(paths)	
		database.UpdateMetrics()
		database.Exploration(screening = True)
	
	
	'''
	calccount = 0
	for i in range(calccount):
		database.EstVoronoi()
		
	x = []
	y = []
	for lib in database.slibs:
		x.append(lib.normalized['fuel_cell_radius'])
		y.append(lib.normalized['enrichment'])
	
	fig, ax = plt.subplots()
	ax.set_xlim([0,1])
	ax.set_ylim([0,1])
	sizes = database.EstVoronoi(s_mult=500)
	sizes = [i*100000 for i in sizes]
	#sizes = 1
	ax.scatter(x, y, s=sizes, c=sizes)
	ax.grid(True)
	#fig.tight_layout()
	plt.show()
	'''
	
	#os.system(str(paths.xsgen_command))
	
	#database.UpdateNeigbors()
	#database.UpdateNeigbors(slib=0)
	
	#database.PCA()
	
	#database.EstLib(database.slibs[8:10], 	t_enrichment=0.01, t_cool_density=0.99)
	
		
	

	'''
	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111, projection='3d')	
	x = []
	y = []
	z = []
	for lib in database.slibs:
		x.append(lib.normalized['fuel_cell_radius'])
		y.append(lib.normalized['enrichment'])
		z.append(lib.normalized['fuel_density'])
	print('Calculating Voronoi volumes')
	sizes = database.EstVoronoi(s_mult=100)
	sizes = [i*100000 for i in sizes]
	max_size = sizes.index(max(sizes))
	print('Calculating neighbors')
	database.UpdateNeigbors(slib=max_size)
	colors = [1 for i in range(len(sizes))]
	colors[max_size] = 10
	for lib_i in database.slib_neighbors[max_size].lib_numbers:
		colors[lib_i] = 5
	max_p = [x[max_size], y[max_size], z[max_size]]
	max_size = sizes[max_size]
	ax1.scatter(x,y,z, s=sizes, c=colors)
	ax1.scatter(max_p[0], max_p[1], max_p[2], s=max_size, c=10, marker='*')
	ax1.set_xlim([0,1])
	ax1.set_ylim([0,1])
	ax1.set_zlim([0,1])
	plt.show()
	'''
	"""
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
	
	print('\n-TheEnd-')
	#input('')
	screen.PrintAt(colors.reset,y=screen.lines)
		
	
	return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
