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

def main(args):
	print "hello world"
	
	databasepath = "/home/cem/nudge/db_dbtest1/"
	libraries = 0
	
	for filename in os.listdir(databasepath):
		
		if filename[0:5] == "build": 
			libraries += 1
			inputfile = databasepath + filename + "/brightlite0/922350.txt"
			
			try:
				doc = open(inputfile, "r")
				
			except IOError:
				print "Could not open ", inputfile
				libraries -= 1
				continue
			
			for line in doc.readlines():
				items = line.split()
				# print items, " ", items[0]
				if items[0] == "NEUT_PROD":
					print "detected!"
		break
	
	print "Total libraries: ", libraries
	
	
	lib1 = objects.Library(path = "path", number = 1)
	
	print "blabla"
	print lib1.fuel_cell_radius
	
	return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
