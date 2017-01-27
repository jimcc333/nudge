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
#   The database folder should contain:
#       - /SR_Inputs 				folder containing all screening inputs
#           - [number].py			input file for scout library [number]
#       - /FR_Inputs 				folder containing all full run inputs
#           - [number].py			input file for full library [number]
#       - /SR_Outputs				folder for all screening output libraries
#           - /build-[number] 		this number must match the one in SR_Inputs
#               - /brightlite0		created by xsgen
#                   - [nucid].txt	output for [nucid] nuclide of input [number]
#       - /FR_Outputs				folder for all full library outputs
#           - /build-[number] 		this number must match the one in FR_Inputs
#               - /brightlite0		created by xsgen
#                   - [nucid].txt	output for [nucid] nuclide of input [number]
#       - basecase.py				xsgen input file containing base-case values
#       - inputs.txt				file containing database inputs
#
#
#   Terms:
#       - Library number: indicated as [number]. Unique number for input-output pair. Starts at zero.
#       - Library progress: screening:[0:1), full=1
#       - Screening library: A library that's run in a short time and that has curtailed outputs
#       - Metric: Names of inputs that libraries get interpolated on
#       - Coordinates: the normalized ([0,1]) input array with only the varied inputs, order based on sorting
#       - Neighborhood: Determined by inputs, the "closest" libs to a given lib (for gradient estimation)
#       - Voronoi cell: The hyper-dimensional "volume" made by points closest to target point
#
#
#   Workflow:
#       1 Start UI and read command line arguments
#       2 Initialize database
#           - If there are inputs, read them; or create the input folders
#               - Attempt to read the output of an input if available
#       3 Screening
#           - Run basecase as screening run
#           - Estimate total time, ask if ok to proceed (improve estimate in the background)
#           - Monte-Carlo inv-norml dist point sampling
#               - Multi-d domain cropping
#               - Scout topography map
#       4 Exploration
#           - Run basecase, space-filling points
#       5 Exploitation
#           - Find highest scored points and inputs near them
#           - Run new points
#           - Estimate max error
#           - Repeat until stop criteria met
#
#   Flags:
#       -m (manual): start NUDGE in manual mode
#       -d (database): used for the database path
#       -h (help):  help screen
#       -x (xsgen): command to run xsgen
#
#
#   Notes:
#       - Folder structure and naming chosen to be simple and intuitive. This enables users to copy-paste
#           their existing libraries and easily allow NUDGE to use it in a given database
#       - xsgen inputs include void and cladding radius, NUDGE also uses thickness in inputs and some workflow
#       - Creation of new library: 1) Generate input file, 2) Initiate library, 3) Add library object to database
#       - Constant inputs will be assigned the value in basecase, which will be the 0th library in database
#       - Dicts in the xsgen input file (initial heavy metal) should be written so that each item is in a new line
#       - During the Voronoi cell volume calculation, best points to use as inputs during the next-batch are saved too
#
#
#
"""" easy copy paste:
import os
from objects import PathNaming
from dbase import DBase
paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\1\\')
database = DBase(paths)
database.update_metrics()
database.plot()
"""

import os

from objects import PathNaming
from dbase import DBase

from pxsgen import *


def main(args):
    os.system('cls' if os.name == 'nt' else 'clear')
    print('----------------- NUDGE: NUclear Database GEneration software -----------------')

    # Check if help is requested
    if '-h' in args:
        print('--- Help ---')
        print('NUDGE is a global surrogate modeling software. It is used to create a database')
        print('of input-output pairs to be used for interpolation for quick estimation of simulations')
        print('Flags:')
        print(' -h: help screen')
        print(' -d [PATH]: generate database at [PATH] based on the input file')
        print(' -explore [PATH]: add one exploration point to database at [PATH]')
        print(' -exploit [PATH]: add one exploration point to database at [PATH]')
        print(' -errors [PATH]: display the recorded errors for database at [PATH], if no [PATH] is given will attempt '
              'to use current directory')
        print()
        print('The database folder should have two files:')
        print(' - basecase.py: xsgen input file containing base-case values')
        print(' - inputs.txt: file containing database inputs')
        return

    if args[1] == '-explore':
        try:
            print('Adding 1 new sample in the database at ' + args[2] + ' using exploration method')
        except IndexError:
            print('Please provide the path of the database')
            return

        paths = PathNaming(os.name, args[2])
        database = DBase(paths)
        database.update_metrics()
        database.exploration()
        return

    if args[1] == '-exploit':
        try:
            print('Adding 1 new sample in the database at ' + args[2] + ' using exploitation method')
        except IndexError:
            print('Please provide the path of the database')
            return

        paths = PathNaming(os.name, args[2])
        database = DBase(paths)

        database.update_metrics()
        database.exploitation()
        return

    if args[1] == '-errors':
        try:
            database_path = args[2]
        except IndexError:
            database_path = os.getcwd()
        print('Reading errors in', database_path)
        read_error_outputs(database_path)
        return

    # Database path
    if '-d' in args:
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\factor00\\', 12, 20, 20, processes=4, record_errors=False)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\factor03\\', 12, 20, 20, processes=4, record_errors=False, exploit_method='guided')
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\factor05\\', 12, 20, 20, processes=4, record_errors=False, exploit_method='guided')
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\factor07\\', 12, 20, 20, processes=4, record_errors=False, exploit_method='guided')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\factor00\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\factor03\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\factor05\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\factor07\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\factor00\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\factor03\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\factor05\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\factor07\\')
        paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\factor00\\0\\')
        database = DBase(paths)
        database.update_metrics()
        database.plot()
        del database
        paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\factor05\\0\\')
        database = DBase(paths)
        database.update_metrics()
        database.plot()
        del database
        return
    # Manual mode check
    if '-m' in args:
        print('Begin database analysis')
        paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\factor00\\0\\')
        database = DBase(paths)
        database.update_metrics()
        database.plot()
        del database
        paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\factor05\\0\\')
        database = DBase(paths)
        database.update_metrics()
        database.plot()
        del database
        paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\factor07\\0\\')
        database = DBase(paths)
        database.update_metrics()
        database.plot()
        return


    print('\n-TheEnd-')

    return 0



