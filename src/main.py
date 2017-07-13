#!/usr/bin/python3
#
#  main.py
#
#  Copyright 2016 Cem Bagdatlioglu <cem@cem-VirtualB>
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
#       - basecase.py				xsgen input file containing base-case values
#       - inputs.txt				file containing database inputs
#   The folder will also have:
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
#
#
#   Terms:
#       - Library number: indicated as [number]. Unique number for input-output pair. Starts at zero.
#       - Library progress: screening:[0:1), full=1
#       - Screening library: A library that's run in a short time and that has curtailed outputs
#       - Varied inputs: Names of inputs that libraries get interpolated on
#       - Coordinates: the normalized ([0,1]) input array with only the varied inputs, ordered alphabetically
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

from dbase import DBase
from repeat import *


def main(args):
    os.system('cls' if os.name == 'nt' else 'clear')
    print('----------------- NUDGE: NUclear Database GEneration software -----------------')

    if '-d' in args:
        # will try functions (f3, f8, f9) and voronoi adjusters (0.0, 0.3, 0.5, 0.8)
        print('begin')
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f33\\', 28, processes=7, exploit_method='guided', record_errors=False, add_new=True)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f35\\', 28, processes=7, exploit_method='guided', record_errors=False, add_new=True)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f38\\', 28, processes=7, exploit_method='guided', record_errors=False, add_new=True)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f3f\\', 28, processes=7, record_errors=False, add_new=True)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f83\\', 28, processes=7, exploit_method='guided', record_errors=False, add_new=True)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f85\\', 28, processes=7, exploit_method='guided', record_errors=False, add_new=True)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f88\\', 28, processes=7, exploit_method='guided', record_errors=False, add_new=True)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f8f\\', 28, processes=7, record_errors=False, add_new=True)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f93\\', 28, processes=7, exploit_method='guided', record_errors=False, add_new=True)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f95\\', 28, processes=7, exploit_method='guided', record_errors=False, add_new=True)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f98\\', 28, processes=7, exploit_method='guided', record_errors=False, add_new=True)
        repeat_databases('C:\\Users\\Cem\\Documents\\nudge\\f9f\\', 28, processes=7, record_errors=False, add_new=True)

        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f33\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f35\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f38\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f3f\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f83\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f85\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f88\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f8f\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f93\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f95\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f98\\')
        find_errors('C:\\Users\\Cem\\Documents\\nudge\\f9f\\')

        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f33\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f35\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f38\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f3f\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f83\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f85\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f88\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f8f\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f93\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f95\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f98\\')
        read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\f9f\\')

        return

    if '-md' in args:

        path = 'C:\\Users\\Cem\\Documents\\nudge\\3D\\'
        database = DBase(path)
        database.exploration(False)
        print()
        del database
        path = 'C:\\Users\\Cem\\Documents\\nudge\\3D1\\'
        database = DBase(path)
        database.timer(57, 0)
        print()
        del database
        path = 'C:\\Users\\Cem\\Documents\\nudge\\3D2\\'
        database = DBase(path)
        database.timer(57, 0)
        print()
        del database
        path = 'C:\\Users\\Cem\\Documents\\nudge\\3D3\\'
        database = DBase(path)
        database.timer(57, 0)
        print()
        del database
        path = 'C:\\Users\\Cem\\Documents\\nudge\\3D4\\'
        database = DBase(path)
        database.timer(57, 0)
        print()
        del database
        return

    if '-2' in args:
        # 2d: 90; 5d: 65; 8d: 35
        count = 12

        path = 'C:\\Users\\Cem\\Documents\\nudge\\2d\\'
        clear_databases(path)
        repeat_databases(path, count, record_errors=False)

        times = []
        for i in range(count):
            database = DBase(path + str(i) + '\\')
            times.append(database.timer(0, 250))
            del database
            print(times[-1])
            print('************************************')

        print(np.mean(times, axis=0))
        print('---------------------------------------------------------------------------------------------------')
        print(np.mean(times, axis=1))

        # path = 'C:\\Users\\Cem\\Documents\\nudge\\3d\\'
        # clear_databases(path)
        # repeat_databases(path, count, record_errors=False)
        #
        # for i in range(count):
        #     database = DBase(path + str(i) + '\\')
        #     database.timer(0, 70)
        #     print()
        #     del database
        #
        # path = 'C:\\Users\\Cem\\Documents\\nudge\\5d\\'
        # clear_databases(path)
        # repeat_databases(path, count, record_errors=False)
        #
        # for i in range(count):
        #     database = DBase(path + str(i) + '\\')
        #     database.timer(0, 60)
        #     print()
        #     del database
        #
        # path = 'C:\\Users\\Cem\\Documents\\nudge\\8d\\'
        # clear_databases(path)
        # repeat_databases(path, count, record_errors=False)
        #
        # for i in range(count):
        #     database = DBase(path + str(i) + '\\')
        #     database.timer(0, 40)
        #     print()
        #     del database


        return


    # Check if help is requested
    if '-h' in args:
        print('--- Help ---')
        print('NUDGE is a global surrogate modeling software. It is used to create a database')
        print('of input-output pairs to be used for interpolation for quick estimation of simulations')
        print('Flags:')
        print(' -h                                        : help screen')
        print(' -d [explore_count] [exploit_count] [PATH] : generate database at [PATH] based on the input file')
        print(' -plot [PATH]                              : plot 2D database at [PATH]')
        print(' -est [PATH]                               : plot 2D database estimate at [PATH]')
        print(' -diff [PATH]                              : plot 2D database difference at [PATH]')
        print(' -explore [PATH]                           : add one exploration point to database at [PATH]')
        print(' -exploit [PATH]                           : add one exploration point to database at [PATH]')
        print(' -errors [PATH]                            : display the recorded errors at [PATH]')
        print(' -study [PATH] [x] [p]                     : perform a database study from [PATH] by repeating')
        print('                                             the database [x] times using [p] processors')
        print()
        print('The database folder should have two files:')
        print(' - basecase.py: xsgen input file containing base-case values')
        print(' - inputs.txt: file containing database inputs')
        return

    if args[1] == '-build':

        database = DBase('C:\\Users\\Cem\\Documents\\nudge\\moxexplore\\')
        database.build(print_progress=True, record_errors=False)

        print('Complete!')
        return

    if args[1] == '-explore':
        try:
            print('Adding 1 new sample in the database at ' + args[2] + ' using exploration method')
        except IndexError:
            print('Please provide the path of the database')
            return

        database = DBase(args[2])
        database.exploration()
        database.run_pxsgen(False)
        database.plot(mark_last=True)
        return

    if args[1] == '-exploit':
        try:
            print('Adding 1 new sample in the database at ' + args[2] + ' using exploitation method')
        except IndexError:
            print('Please provide the path of the database')
            return

        database = DBase(args[2])
        database.exploitation()
        database.run_pxsgen(False)
        database.plot(mark_last=True)
        return

    if args[1] == '-m':
        counter = 0
        database = DBase('C:\\Users\\Cem\\Documents\\nudge\\mox50\\')
        print('Read database')
        for lib in database.flibs:
            inputs = []
            inputs.append(0.00004343 + lib.inputs.xsgen['fuel_cell_radius'] * 2.17148E-05)
            inputs.append(0.00003729 + lib.inputs.xsgen['void_cell_radius'] * 0.000018644)
            inputs.append(0.00081248 + lib.inputs.xsgen['clad_cell_radius'] * 0.00040624)
            inputs.append(0.00038604 + lib.inputs.xsgen['unit_cell_height'] * 0.00019302)
            inputs.append(0.00013993 + lib.inputs.xsgen['fuel_density'] * 0.000069964)
            inputs.append(0.00010561 + lib.inputs.xsgen['clad_density'] * 0.000052804)
            total = sum(inputs)
            inputs.append(0.0232930 - total)
            ipfile = database.base_file
            ipfile = ipfile.replace('aa1', str(inputs[0]))
            ipfile = ipfile.replace('aa2', str(inputs[1]))
            ipfile = ipfile.replace('aa3', str(inputs[2]))
            ipfile = ipfile.replace('aa4', str(inputs[3]))
            ipfile = ipfile.replace('aa5', str(inputs[4]))
            ipfile = ipfile.replace('aa6', str(inputs[5]))
            ipfile = ipfile.replace('aa7', str(inputs[6]))
            with open('C:\\Users\\Cem\\Documents\\nudge\\mox50\\inputs\\' + str(counter) + '.py', 'w') as openfile:
                openfile.write(ipfile)
            openfile.close()
            counter += 1
        return

    if args[1] == '-errors':
        try:
            database_path = args[2]
        except IndexError:
            database_path = os.getcwd()
        print('Reading errors in', database_path)
        read_error_outputs(database_path)
        return

    if args[1] == '-d' and len(args) == 5:
        print('Building database', args[4], 'with', args[2], 'exploration and', args[3], 'exploitation samples.')
        database = DBase(args[4])
        database.build(int(args[2]), int(args[3]), print_progress=True)
        database.plot()
        return

    if args[1] == '-plot':
        print('Plotting database', args[2])
        database = DBase(args[2])
        database.update_coordinates()
        database.plot(numbers=True)
        return

    if args[1] == '-est':
        print('Plotting estimate from database', args[2])
        database = DBase(args[2])
        database.update_metrics()
        database.plot_estimate()
        return

    if args[1] == '-diff':
        print('Plotting difference of database', args[2])
        database = DBase(args[2])
        database.update_metrics()
        database.plot_estimate(diff=True)
        return

    if args[1] == '-study':
        print('Beginning database study from directory:', args[2])
        print('Repeating', args[3], 'times, using', args[4], 'processors')
        repeat_databases(args[2], args[3], 0, 0, processes=args[4])
        return


