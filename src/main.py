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
paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\11\\')
database = DBase(paths)
database.update_metrics()
database.plot()
"""

import os
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
        paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\1\\')
        database = DBase(paths)
        database.update_metrics()
        #database.build(20, 20)
        database.find_error(print_result=True)
        database.plot()
        return

    # Manual mode check
    if '-m' not in args:
        #TODO: have non-manual mode
        pass
    else:

        # Manual mode

        pool = Pool(processes=6)
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



        """
        # Plot data
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        X = np.arange(0, 1, 0.01)
        Y = np.arange(0, 1, 0.01)
        X, Y = np.meshgrid(X, Y)

        R = np.sqrt(X ** 2 + Y ** 2)

        Z = []
        for i in range(len(X)):
            vector = []
            for j in range(len(Y)):
                vector.append(f6(X[i][j], Y[i][j], 0.5))
            Z.append(vector)

        #Z = np.sin(R)

        print('done')
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.set_zlim(-1.01, 1.01)

        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()
        """
        '''
        # Find correlation matrix
        a = np.array([x1, x2, ..., xd, y1, y2, y3])
        corr = np.corrcoef(a)
        print(corr)
        '''

    print('\n-TheEnd-')

    return 0


def database_thread(database_path, exploration_count, exploitation_count):
    paths = PathNaming(os.name, database_path=database_path)
    database = DBase(paths)
    database.update_metrics()
    database.build(exploration_count, exploitation_count)

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
