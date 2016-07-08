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
#		- /Inputs 				(folder containing all inputs)
#			- [tag][number].txt		(input file for library [number])
#		- /build-[tag][number] 
#			- /brightlite0
#				- [nucid].txt	(output for [nucid] nuclide of input [number])
#  
#  

from objects import *
import os

import numpy as np
from matplotlib.mlab import PCA as mlabPCA

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
				inputlib = Library(databasepath, databasepath + "Inputs/" + inputfile, int(inputfile[:-3]))
				
				# Add lib to database
				if not database.Exists(inputlib.number):
					database.AddSLib(inputlib)
					print(" Added lib #" + str(inputlib.number) + " to database")
				
				break
			break
	
	
	print("Total libraries: " + str(len(database.slibs)))
	#database.UpdateData()	
	#database.PCA()
	
	"""
	pca_mat = mlabPCA(database.data_mat)
	
	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111, projection='3d')
	ax1.scatter(database.data_mat[:,0], database.data_mat[:,1], database.data_mat[:,2])
	ax1.set_xlabel("Neutron Production")
	ax1.set_xlabel("Neutron Destruction")
	ax1.set_zlabel("Burnup")
	
	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111, projection='3d')
	ax2.scatter(pca_mat.Y[:,0], pca_mat.Y[:,1], pca_mat.Y[:,2])
	ax2.set_xlabel("PC1")
	ax2.set_xlabel("PC2")
	ax2.set_zlabel("PC3")
	plt.show()"""

	return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
