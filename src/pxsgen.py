#!/usr/bin/python3
#
#  
#  Copyright 2016 Cem Bagdatli
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
#  Inputs:
#  	1. enrichment
#
#
#  
#  Notes:
#	- Does not use fuel composition: just enrichment
#	- All inputs are [0,1]
#
#
#
#

from objects import xsgenParams, Library
import os
import random
import numpy as np


def main(args):
	os.system('cls' if os.name == 'nt' else 'clear')
	print('Placeholder xsgen program')
	
	
	# Inputs and constants
	lib = Library('blank', 'blank', args[1], 1, False)
	
	outputs = []
	
	print('\n-TheEnd-')	
	return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
