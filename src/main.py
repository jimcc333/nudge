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
        print('Help request')
        delete_after('C:\\software\\nudge\\f8_120_s15\\', 119)
        return

    if '-e' in args:
        read_error_outputs('C:\\software\\nudge\\f8_240_s0\\')
        read_error_outputs('C:\\software\\nudge\\f8_120_s10\\')
        read_error_outputs('C:\\software\\nudge\\f8_120_s15\\')

        return

        for path in ['C:\\software\\nudge\\f8_300\\0\\', 'C:\\software\\nudge\\f8_300_s10\\0\\', 'C:\\software\\nudge\\f8_300_s15\\0\\']:
            paths = PathNaming(os.name, database_path=path)
            database = DBase(paths)
            database.update_metrics()
            database.plot(est_errors=True)

        return
        return

    # Database path
    if '-d' in args:
        # find_errors('C:\\software\\nudge\\f8_240_s0\\')
        find_errors('C:\\software\\nudge\\f8_120_s10\\', exclude_after=120)
        # read_error_outputs('C:\\software\\nudge\\f8_120_s15\\')

        return

        # repeat_databases('C:\\software\\nudge\\f8_240_s0\\', 18, 120, 120, processes=6, record_errors=False)
        for i in range(6):
            repeat_databases('C:\\software\\nudge\\f8_120_s10\\', 12, 10, 10, processes=6, record_errors=False)
        for i in range(4):
            repeat_databases('C:\\software\\nudge\\f8_120_s15\\', 12, 15, 15, processes=6, record_errors=False)
        for i in range(2):
            repeat_databases('C:\\software\\nudge\\f8_120_s10\\', 12, 30, 30, processes=6, record_errors=False)

        return
    # Manual mode check
    if '-m' in args:
        print('Begin database analysis')
        paths = PathNaming(os.name, database_path='C:\\software\\nudge\\f8_120_s10\\0\\')
        database = DBase(paths)
        database.update_metrics()
        database.plot(points=True)
        return
        database.build(10,10, print_progress=True, record_errors=False)
        return
        database.estimate_error(plot=True)
        return
        database.plot_estimate(diff=True, exclude_after=50, abs_max=0.3)
        database.plot_estimate(diff=True, exclude_after=150, abs_max=0.3)
        database.plot_estimate(diff=True, abs_max=0.3)

        return

        it_range = list(range(100))
        it_range = it_range[10:]
        for i in it_range:
            database.estimate_error(exclude_after=i, print_result=True)
            continue

        return

    random_data = read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\random1\\')
    explore_only = read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\0_2001\\')
    exploit_40 = read_error_outputs('C:\\Users\\Cem\\Documents\\nudge\\40_160\\')
    """
    plt.plot(random_data[:-1])
    plt.plot(explore_only[:-1])
    plt.plot(exploit_20[:-1])
    plt.plot(exploit_40[:-1])
    plt.plot(exploit_60[:-1])
    plt.plot(exploit_80[:-1])

    plt.legend(['Random', 'Explore Only', 'Exploit 20', 'Exploit 40', 'Exploit 60', 'Exploit 80'])
    x_max = int(min([len(random_data), len(explore_only), len(exploit_20)])) - 1
    y_max = max([max(random_data), max(explore_only), max(explore_only)])
    y_max *= 1.05
    plt.axis([15, x_max, 0, 0.7])
    plt.xlabel('# of points in database')
    plt.ylabel('Error (%)')
    plt.show()
    """
    print('\n-TheEnd-')

    return 0


def repeat_databases(source_path, database_count, exploration_count, exploitation_count, random_count=0, processes=7,
                     exploit_method='furthest', record_errors=True):
    # Generate threading lists
    paths = PathNaming(os.name, database_path=source_path)
    database_paths = [source_path + str(i) + paths.slash for i in range(database_count)]
    explorations = [exploration_count for i in range(database_count)]
    exploitations = [exploitation_count for i in range(database_count)]
    randoms = [random_count for i in range(database_count)]
    exploit_method = [exploit_method for i in range(database_count)]
    record_errors = [record_errors for i in range(database_count)]

    # Make a new folder for each database and place the input files in it
    for i in range(database_count):
        if not os.path.isdir(database_paths[i]):
            os.mkdir(database_paths[i])
            shutil.copy(source_path + paths.base_input, database_paths[i])
            shutil.copy(source_path + paths.dbase_input, database_paths[i])

    # Run databases
    pool = Pool(processes=processes)
    pool.starmap(database_thread, zip(database_paths, explorations, exploitations, randoms, exploit_method,
                                      record_errors))

    return


# Goes through a database study and builds errors for each database
def find_errors(source_path, find_all=False, exclude_after=None):
    try:
        folders = os.listdir(source_path)
    except FileNotFoundError:
        print('The database study directory', source_path, 'does not exist')
        return
    slash = '\\' if os.name == 'nt' else '/'

    for folder_name in folders:
        try:
            paths = PathNaming(os.name, database_path=source_path+folder_name+slash)
            database = DBase(paths)
            print('Finding errors of', database.paths.database_path)
            database.update_metrics()
            database.estimate_error(exclude_after=exclude_after)
            database.find_error()
            database.write_errors()
            del database
        except (FileNotFoundError, NotADirectoryError):
            continue

    return


# Reads errors of databases inside a folder
def read_error_outputs(source_path):
    slash = '\\' if os.name == 'nt' else '/'
    folders = os.listdir(source_path)
    file_count = 0
    max_errors = []
    min_errors = []
    mean_errors = []
    real_max = []
    real_errors = []

    for folder_name in folders:
        try:
            doc = open(source_path + folder_name + slash + 'errors.txt', "r")
            file_count += 1
            lines = doc.readlines()
            max_errors.append([float(i) for i in lines[1][1:-2].split()])
            min_errors.append([float(i) for i in lines[3][1:-2].split()])
            mean_errors.append([float(i) for i in lines[5][1:-2].split()])
            real_max.append([float(i) for i in lines[7][1:-2].split()])
            real_errors.append([float(i) for i in lines[9][1:-2].split()])

        except FileNotFoundError:
            continue
    if file_count == 0:
        return []

    # plt.plot(np.mean(real_errors, axis=0))
    # plt.show()

    print(source_path)
    print(np.mean(max_errors, axis=0))
    print()
    print(np.mean(mean_errors, axis=0))
    print()
    print(np.mean(real_max, axis=0))
    print()
    print(np.mean(real_errors, axis=0))

    return np.mean(real_errors, axis=0)


# Deletes all libraries in a database study after the specified number
def delete_after(database_path, number):
    path = PathNaming(os.name, database_path=database_path)
    folders = os.listdir(database_path)

    for folder_name in folders:
        try:
            ip_path = database_path + folder_name + path.slash + path.FR_Input_folder
            op_path = database_path + folder_name + path.slash + path.FR_Output_folder
            libraries = os.listdir(ip_path)
            for lib in libraries:
                lib_i = int(lib.split('.')[0])
                if lib_i > number:
                    print('Deleting', lib)
                    os.remove(ip_path + path.slash + lib)
                    os.remove(op_path + path.slash + lib)
        except FileNotFoundError:
            continue

def database_thread(database_path, exploration_count, exploitation_count, random_count, exploit_method, record_errors):
    paths = PathNaming(os.name, database_path=database_path)
    database = DBase(paths)
    database.update_metrics()
    if random_count > 0:
        database.random_selection(random_count)
    else:
        database.build(exploration_count, exploitation_count, record_errors=record_errors,
                       exploit_method=exploit_method)

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
