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
#		- /SR_Inputs 				(folder containing all scouting inputs)
#			- [number].txt			(input file for scout library [number])
#		- /FR_Inputs 				(folder containing all full run inputs)
#			- [number].txt			(input file for full library [number])
#		- SR_Outputs				(folder for all scouting output libraries)
#			- /build-[number] 		(this number must match the one in SR_Inputs)
#				- /brightlite0		(created by xsgen)
#					- [nucid].txt	(output for [nucid] nuclide of input [number])
#		- FR_Outputs				(folder for all full library outputs)
#			- /build-[number] 		(this number must match the one in FR_Inputs)
#				- /brightlite0		(created by xsgen)
#					- [nucid].txt	(output for [nucid] nuclide of input [number])
#		- database_inputs.py		(file containing database inputs)
#  
#
#	Terms:
#		- Scout library: A library thats run in a short time and that has curtailed outputs
#		- Metric: Values that libraries get interpolated on
#  
#TODO: let users input changes in xsgen-dependent constants (such as brightlite0 folder)

from objects import *
import os

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main(args):
	print("hello world")
	
	# inputs to be taken from user
	databasepath = "/home/cem/nudge/db_dbtest1/"
	
	database = DBase("database1")
	
	#TODO: first read inputs and create a lib for each ID
	#	at lib creation attempt to read lib output
	#	at the end of reading all inputs read outputs and do consistency check
	for filename in os.listdir(databasepath):
		if filename[:6] == "Inputs":
			print("Reading libraries")
			for inputfile in os.listdir(databasepath + "Inputs/"):
				# Initialize library
				print(inputfile, '  ', inputfile.split())
				lib_number = [int(s) for s in inputfile.split() if s.isdigit()]
				inputlib = Library(databasepath, databasepath + "Inputs/" + inputfile, lib_number[-1])
				
				# Add lib to database
				if not database.Exists(inputlib.number):
					database.AddSLib(inputlib)
					print(" Added lib #" + str(inputlib.number) + " to database")
			break
	
	
	print("Total libraries: " + str(len(database.slibs)))
	database.UpdateData()
	#database.Print()
	database.PCA()
	
	database.UpdateMetrics();
	database.EstLib(database.slibs[3:7], t_fuel_radius=1, t_fuel_density=0.3, t_cool_density=0, t_enrichment=0.8)
	
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
	
	return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
