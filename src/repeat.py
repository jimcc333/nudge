import os
import shutil
from distutils.dir_util import copy_tree
from multiprocessing import Pool

import numpy as np
from math import floor
import matplotlib.pyplot as plt

from objects import PathNaming
from dbase import DBase


# Repeats databases from the same inputs and basecase
def repeat_databases(source_path, database_count, exploration_count=0, exploitation_count=0, random_count=0, processes=7,
                     add_new=False, exploit_method='furthest', record_errors=True, start_number=0):
    # Check source path
    if not os.path.isdir(source_path):
        print('Unable to find source:' + source_path)
        return

    # Generate threading lists
    if add_new:
        items = os.listdir(source_path)
        databases = []
        for item in items:
            try:
                databases.append(int(item))
            except ValueError:
                continue
        if len(databases) != 0:
            start_number = max(databases) + 1
    end_number = start_number + database_count
    paths = PathNaming(os.name, database_path=source_path)
    database_paths = [source_path + str(i) + paths.slash for i in range(start_number, end_number)]
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

    # If there's existing libraries to begin with, copy them as well
    if os.path.isdir(source_path + paths.slash + paths.FR_Input_folder + paths.slash):
        source_dir = source_path + paths.slash + paths.FR_Input_folder + paths.slash
        for i in range(database_count):
            copy_tree(source_dir, database_paths[i] + paths.slash + paths.FR_Input_folder)
    if os.path.isdir(source_path + paths.slash + paths.FR_Output_folder + paths.slash):
        source_dir = source_path + paths.slash + paths.FR_Output_folder + paths.slash
        for i in range(database_count):
            copy_tree(source_dir, database_paths[i] + paths.slash + paths.FR_Output_folder)

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
    paths = PathNaming(os.name)

    for folder_name in folders:
        if os.path.exists(source_path+folder_name+slash+paths.error_file):
            print('Errors file exists, skipping', source_path+folder_name)
            continue
        try:
            database = DBase(source_path+folder_name+slash)
            print('Finding errors of', database.paths.database_path)
            database.update_metrics()
            database.estimate_error(exclude_after=exclude_after)
            database.find_error()
            database.write_errors()
            del database
        except (FileNotFoundError, NotADirectoryError, RuntimeError):
            continue

    return


# Reads errors of databases inside a folder
def read_error_outputs(source_path, final_error_only=False):
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
    if len(real_errors[0]) != len(mean_errors[0]):
        print('Warning. Errors inconsistent.')

    prev_count = len(real_errors[0])
    for errors in real_errors[1:]:
        if len(errors) != prev_count:  # numpy will catch this but this error message will make it clearer
            print('Warning. Number of reported errors vary between databases')

    if final_error_only:
        print(source_path, 'Databases:', len(real_errors))
        print('  Final real error:', round(np.mean(real_errors, axis=0)[-1], 3), '+/-',
              round(np.std(real_errors, axis=0)[-1], 4))
        return

    print(source_path, 'Databases:', len(real_errors))
    print(np.mean(max_errors, axis=0))
    print(source_path)
    print('estimated max:', np.mean(max_errors, axis=0))
    print()
    print('estimated errors:', np.mean(mean_errors, axis=0))
    print()
    print('real max:', np.mean(real_max, axis=0))
    print()
    print('real errors:', np.mean(real_errors, axis=0))
    print('real errors STD:', np.std(real_errors, axis=0))
    print('-----------------------------------------------------------------------------------------')


# Plots heat maps
def plot_heat_map(source_path, resolution=50):
    print('Plotting heat map of', source_path)
    slash = '\\' if os.name == 'nt' else '/'
    folders = os.listdir(source_path)
    counter = 0

    # Generate a grid
    colors = np.zeros((resolution, resolution))

    # Go through each database and its samples assuming all inputs [0,1]
    print('Reading databases:', end='  ')
    for folder_name in folders:
        try:
            database = DBase(PathNaming(os.name, database_path=source_path + folder_name + slash))
            database.update_metrics()
            samples = [coords for coords in database.lib_inputs]
            counter += 1
            print(folder_name, end=' ')
            for i in samples:
                x = floor(resolution*i[0])
                y = floor(resolution*i[1])
                if x == resolution:
                    x -= 1
                if y == resolution:
                    y -= 1
                colors[y, x] += 1
        except (NotADirectoryError, FileNotFoundError):
            continue
    print('(', counter, ' databases total)', sep='')

    colors = np.divide(colors, colors.max())
    fig, ax = plt.subplots()
    ax.set_xlim([-0.01, 1.01])
    ax.set_ylim([-0.01, 1.01])
    plt.imshow(colors, extent=(0, 1, 0, 1), origin='lower', interpolation='nearest')
    plt.show()


# Deletes all libraries in a database study (or one database) after the specified number
def delete_after(study_path, number, database_path=None):
    if database_path is not None:
        path = PathNaming(os.name, database_path=database_path)
        ip_path = database_path + path.slash + path.FR_Input_folder
        try:
            op_path = database_path + path.slash + path.FR_Output_folder
            libraries = os.listdir(ip_path)
            for lib in libraries:
                lib_i = int(lib.split('.')[0])
                if lib_i > number:
                    print('Deleting lib', lib, 'in', database_path)
                    os.remove(ip_path + path.slash + lib)
                    os.remove(op_path + path.slash + lib)
        except FileNotFoundError:
            return
        return

    path = PathNaming(os.name, database_path=study_path)
    folders = os.listdir(study_path)
    for folder_name in folders:
        try:
            ip_path = study_path + folder_name + path.slash + path.FR_Input_folder
            op_path = study_path + folder_name + path.slash + path.FR_Output_folder
            libraries = os.listdir(ip_path)
            for lib in libraries:
                lib_i = int(lib.split('.')[0])
                if lib_i > number:
                    print('Deleting', lib)
                    os.remove(ip_path + path.slash + lib)
                    os.remove(op_path + path.slash + lib)
        except FileNotFoundError:
            continue


# Runs a database for threading in database study
def database_thread(database_path, exploration_count=0, exploitation_count=0, random_count=0, exploit_method='furthest',
                    record_errors=True):
    database = DBase(database_path)
    if random_count > 0:
        database.random_selection(random_count)
    else:
        database.build(exploration_to_add=exploration_count, exploitation_to_add=exploitation_count,
                       record_errors=record_errors, exploit_method=exploit_method)
