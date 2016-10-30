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
import shutil
from multiprocessing import Pool

from objects import PathNaming
from dbase import DBase

from pxsgen import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm
import numpy as np

def main(args):
    os.system('cls' if os.name == 'nt' else 'clear')

    # Check if help is requested
    if '-h' in args:
        return

    # Database path
    if '-d' in args:



        paths = PathNaming(os.name, database_path='C:\\Users\\cb39852\\Documents\\nudge\\2D_random\\')
        database = DBase(paths)
        database.update_metrics()
        database.plot()
        return

        paths = PathNaming(os.name, database_path='C:\\Users\\cb39852\\Documents\\nudge\\2D_exploit\\')
        database = DBase(paths)
        database.update_metrics()
        database.build(20, 40, print_progress=True)

        paths = PathNaming(os.name, database_path='C:\\Users\\cb39852\\Documents\\nudge\\2D_random\\')
        database = DBase(paths)
        database.update_metrics()
        database.random_selection(count=63, print_progress=True)

    # Manual mode check
    if '-m' in args:
        print('Begin database repeat mode')
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\test\\', 3, 10, 3)
        return
    else:

        # Manual mode

        pool = Pool(processes=7)
        paths = ['C:\\Users\\Cem\\Documents\\nudge\\1\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\2\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\3\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\4\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\5\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\6\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\7\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\8\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\9\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\10\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\11\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\12\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\13\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\14\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\15\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\16\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\17\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\18\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\19\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\20\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\21\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\22\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\23\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\24\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\25\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\26\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\27\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\28\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\29\\',
                 'C:\\Users\\Cem\\Documents\\nudge\\30\\']
        explorations = [97 for i in range(10)] + [72 for i in range(10)] + [47 for i in range(10)]
        exploitations = [0 for i in range(10)] + [25 for i in range(10)] + [50 for i in range(10)]

        if len(paths) == len(explorations) and len(paths) == len(exploitations):
            pass
        else:
            print('Input error dude')
            return

        pool.starmap(database_thread, zip(paths, explorations, exploitations))

    print('\n-TheEnd-')

    return 0


def repeat_databases(source_path, database_count, exploration_count, exploitation_count, random_count=0, processes=7):
    # Generate threading lists
    paths = PathNaming(os.name, database_path=source_path)
    database_paths = [source_path + str(i) + paths.slash for i in range(database_count)]
    explorations = [exploration_count for i in range(database_count)]
    exploitations = [exploitation_count for i in range(database_count)]
    randoms = [random_count for i in range(database_count)]

    # Make a new folder for each database and place the input files in it
    for i in range(database_count):
        os.mkdir(database_paths[i])
        shutil.copy(source_path + paths.base_input, database_paths[i])
        shutil.copy(source_path + paths.dbase_input, database_paths[i])

    # Run databases
    pool = Pool(processes=7)
    pool.starmap(database_thread, zip(database_paths, explorations, exploitations, randoms))

    return


def database_thread(database_path, exploration_count, exploitation_count, random_count):
    paths = PathNaming(os.name, database_path=database_path)
    database = DBase(paths)
    database.update_metrics()
    if random_count > 0:
        database.random_selection(random_count)
    else:
        database.build(exploration_count, exploitation_count)

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
