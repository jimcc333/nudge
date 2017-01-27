import os
import shutil
from multiprocessing import Pool

import numpy as np

from objects import PathNaming
from dbase import DBase


# Repeats databases from the same inputs and basecase
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


# Runs a database for threading in database study
def database_thread(database_path, exploration_count, exploitation_count, random_count, exploit_method, record_errors):
    paths = PathNaming(os.name, database_path=database_path)
    database = DBase(paths)
    database.update_metrics()
    if random_count > 0:
        database.random_selection(random_count)
    else:
        database.build(exploration_count, exploitation_count, record_errors=record_errors,
                       exploit_method=exploit_method)