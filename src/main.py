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
#  
#  
#  

import os
from objects import *
from parser import *

import numpy as np
from matplotlib.mlab import PCA as mlabPCA

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main(args):
	print "hello world"
	
	databasepath = "/home/cem/nudge/db_dbtest1/"
	database = objects.DBase("database1")
	
	libraries = 0
	for filename in os.listdir(databasepath):
		if filename[0:5] == "build": 
			libraries += 1
			lib = objects.Library(databasepath + filename, libraries)
			database.Add(lib)
	
	print "Total libraries: ", libraries
	database.UpdateData()	
	database.PCA()
	
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
	plt.show()
	
	return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
