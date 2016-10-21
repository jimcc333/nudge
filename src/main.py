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

import os

from objects import PathNaming
from dbase import DBase


def main(args):
    os.system('cls' if os.name == 'nt' else 'clear')

    # Check if help is requested
    if '-h' in args:
        return

    # Manual mode check
    if '-m' not in args:
        #TODO: have non-manual mode
        pass
    else:
        # Manual mode

        paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\1\\')
        database = DBase(paths)
        database.update_metrics()

        database.build(50, 0)

        paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\2\\')
        database = DBase(paths)
        database.update_metrics()

        database.build(40, 10)

        paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\3\\')
        database = DBase(paths)
        database.update_metrics()

        database.build(30, 20)

        paths = PathNaming(os.name, database_path='C:\\Users\\Cem\\Documents\\nudge\\4\\')
        database = DBase(paths)
        database.update_metrics()

        database.build(25, 25)


        '''
        # Plot data
        x = []
        y = []
        z = []
        for i, lib in enumerate(database.flibs):
            x.append(lib.normalized['fuel_density'])
            y.append(lib.normalized['clad_density'])
            z.append(lib.normalized['cool_density'])
            print(i, lib.rank)
        x.append(next_lib.furthest_point[0])
        y.append(next_lib.furthest_point[1])

        fig, ax = plt.subplots()
        ax.set_xlim([-0.1,1.1])
        ax.set_ylim([-0.1,1.1])
        labels = [i for i in range(len(x))]
        #ax.scatter(rx, ry, s=20, c='g')
        #ax.scatter(p_cand[0], p_cand[1], s=200, c='r')
        ax.grid(True)
        fig.tight_layout()
        color_l = []

        for i, txt in enumerate(labels):
            #label_i = 'p'+str(txt)+', v:'+str(sizes[i])+'\n nonlin:'+str(round(database.flibs[i].neighborhood.nonlinearity,1))
            label_i = 'p'+str(txt)
            if i < len(database.flibs):
                color_l.append(database.flibs[i].neighborhood.nonlinearity)
            else:
                label_i = 'NEXT'
            ax.annotate(label_i, (x[i],y[i]), xytext = (x[i]-0.07,y[i]+0.05))
        ax.scatter(x, y, s=200, c='b')
        plt.show()
        '''

        '''
        # Find correlation matrix
        a = np.array([x1, x2, ..., xd, y1, y2, y3])
        corr = np.corrcoef(a)
        print(corr)
        '''

    print('\n-TheEnd-')

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
