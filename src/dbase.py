import os
import copy
import itertools
import random
import subprocess

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import PCA as mlabPCA
from scipy.spatial import distance
from scipy.interpolate import griddata

from library import Library
from objects import xsgenParams, Neighborhood
from pxsgen import burnup_maker, prod_maker, dest_maker, main

"""
A class that handles all generated libraries

Variables:
    inputs							# Database inputs
    self.ip_ranges = xsgenParams()	# Ranges of varied inputs
    self.varied_ips
    slibs, flibs
    complete_slibs, complete_flibs
    flib_neighbors, neighbor_score


"""


class DBase:

    def __init__(self, paths):
        # Read database, assuming software may have been interrupted and
        # the folder may have some inputs and outputs

        # Database constants
        self.paths = paths
        self.varied_ips = []            # the varied inputs
        self.ip_ranges = xsgenParams()  # Ranges of varied inputs
        self.proj_threshold = 0.0001    # Projection check (exploration) threshold
        self.dimensions = None          # Number of varied inputs for database (assigned when inputs are read)
        self.basecase = None            # The basecase Library object (assigned when basecase is read)
        self.base_file = None           # The basecase as string (assigned when basecase is read)

        # Database inputs
        self.inputs = {
            'max_error': None,          # In [%]
            'max_time': 100,            # In [hour]
            'scout_frac': 10,           # Weight of screening time allocation
            'explore_frac': 40,         # Weight of exploration time allocation
            'exploit_frac': 50,         # Weight of exploitation time allocation
        }

        # Database libraries
        self.slibs = []                 # Screening libraries (with and without outputs)
        self.flibs = []                 # Full libraries (with and without outputs)

        # Database analysis
        self.lib_inputs = []            # (d,n) matrix of varied input normalized coordinates (n: number of flibs)
        self.lib_outputs = []           # (n) vector of library outputs
        self.voronoi_sizes = []         # Voronoi cell sizes of points in the database
        self.np_prods = None            # numpy array of neutron production values
        self.np_dests = None            # numpy array of neutron destruction values
        self.np_BUs = None              # numpy array of burnup values
        self.data_mat = []              # numpy data matrix of neutron prod/dest and BU
        self.pca_mat = None             # numpy PCA matrix
        self.est_error_max = []         # Estimated maximum database error vector
        self.est_error_min = []         # Estimated minimum database error vector
        self.est_error_mean = []        # Estimated mean database error vector
        self.database_error = []        # "True" database error vector calculated using generated test points

        # Check that database path exists
        if not os.path.isdir(paths.database_path):
            error_message = 'The database path does not exist. Looking for: ' + paths.database_path
            raise RuntimeError(error_message)

        # Check that database input file exists
        if not os.path.exists(paths.database_path + paths.dbase_input):
            error_message = 'The database input file does not exist. Looking for: ' \
                            + paths.database_path + paths.dbase_input
            raise RuntimeError(error_message)

        # Check that basecase input exists
        if not os.path.exists(paths.database_path + paths.base_input):
            error_message = 'The database base-case input file does not exist. Looking for: ' \
                            + paths.database_path + paths.base_input
            raise RuntimeError(error_message)

        self.name = paths.database_name

        # Read database input file
        self.read_input(paths.database_path + paths.dbase_input)

        # Read basecase input
        self.read_base(paths.database_path + paths.base_input)

        # Check to see if there's a screening library input folder
        if os.path.exists(paths.database_path + paths.SR_Input_folder):
            tot_sfiles = len(os.listdir(paths.database_path + paths.SR_Input_folder))
        # If the input screening library folder doesn't exist create it
        else:
            os.mkdir(paths.database_path + paths.SR_Input_folder)
            tot_sfiles = 0

        # Check to see if there's a full library input folder
        if os.path.exists(paths.database_path + paths.FR_Input_folder):
            tot_ffiles = len(os.listdir(paths.database_path + paths.FR_Input_folder))
        # If the input full library folder doesn't exist create it
        else:
            os.mkdir(paths.database_path + paths.FR_Input_folder)
            tot_ffiles = 0

        tot_sr_libs = 0
        tot_fr_libs = 0

        # If there are files SR_Inputs, read them
        if tot_sfiles > 0:
            for ip_number in range(tot_sfiles):
                    ip_path = paths.database_path + paths.SR_Input_folder + paths.slash + str(ip_number) + '.py'
                    #TODO: op_path will need to be updated to work with xsgen
                    op_path = paths.database_path + paths.SR_Output_folder + paths.slash + str(ip_number) + '.py'
                    #+ paths.xsgen_prefix + paths.sr_prefix + str(ip_number) + '/' + paths.xsgen_op_folder

                    if os.path.exists(ip_path):
                        inputlib = Library(database_path=paths.database_path, op_path=op_path, ip_path=ip_path,
                                           number=ip_number, scout=True)
                        self.slibs.append(copy.deepcopy(inputlib))
                        tot_sr_libs += 1
                    else:
                        # Stops reading inputs if they're not sequential. This will overwrite inputs
                        break

        # If there are files FR_Inputs, read them
        if tot_ffiles > 0:
            for ip_number in range(tot_ffiles):
                ip_path = paths.database_path + paths.FR_Input_folder + paths.slash + str(ip_number) + '.py'
                # TODO: op_path will need to be updated to work with xsgen
                op_path = paths.database_path + paths.FR_Output_folder + paths.slash + str(ip_number) + '.py'
                # + paths.xsgen_prefix + paths.sr_prefix + str(ip_number) + '/' + paths.xsgen_op_folder

                if os.path.exists(ip_path):
                    inputlib = Library(database_path=paths.database_path, op_path=op_path, ip_path=ip_path,
                                       number=ip_number, scout=True)
                    self.flibs.append(copy.deepcopy(inputlib))
                    tot_fr_libs += 1
                else:
                    # Could continue here instead of breaking, but at this point this is better
                    break

        # If there's no output folder, create it
        if not os.path.exists(paths.database_path + paths.SR_Output_folder):
            os.mkdir(paths.database_path + paths.SR_Output_folder)
        if not os.path.exists(paths.database_path + paths.FR_Output_folder):
            os.mkdir(paths.database_path + paths.FR_Output_folder)

    # Once normalized coords are generated, pass here to add
    #TODO: outputting in correct units
    def add_lib(self, norm_coords, screening):
        # Adding a library will:
        # - Convert normalized coords to correct units
        # - Generate the input file
        # - Add to flib or slib list in this object

        # Dimension consistency check (note that this is technically not sufficient)
        if len(norm_coords) != self.dimensions:
            error_message = 'Dimensions mismatch when adding a new library to database'
            raise RuntimeError(error_message)

        # Convert normalized to correct unit value
        new_inputs = copy.copy(self.basecase.inputs.xsgen)
        for i, key in enumerate(sorted(self.varied_ips)):
            new_inputs[key] = norm_coords[i] #* self.basecase.inputs.xsgen[key]

        # Check if input exists
        for lib in (self.slibs if screening else self.flibs):
            if lib.inputs.xsgen == new_inputs:
                print('Input already exists')
                return

        # Create paths
        lib_number = len(self.flibs)
        if screening:
            lib_number = len(self.slibs)
        if screening:
            ip_path = self.paths.database_path + self.paths.SR_Input_folder + self.paths.slash + str(lib_number) + '.py'
            op_path = self.paths.database_path + self.paths.SR_Output_folder + self.paths.slash + str(lib_number)+'.py'
        else:
            ip_path = self.paths.database_path + self.paths.FR_Input_folder + self.paths.slash + str(lib_number) + '.py'
            op_path = self.paths.database_path + self.paths.FR_Output_folder + self.paths.slash + str(lib_number)+'.py'

        # Make input file from basecase file
        ipfile = self.base_file

        # Name the input file reactor by string replacement
        base_line = [line for line in ipfile.split('\n') if 'reactor = ' in line]
        if len(base_line) != 1:
            error_message = 'Input creation error when building from basecase. ' +\
                            str(len(base_line)) + ' instances of: ' + 'reactor = '
            raise RuntimeError(error_message)
        ipfile = ipfile.replace(base_line[0], 'reactor = ' + str(lib_number))

        # Change input variables by string replacement
        for key, value in new_inputs.items():
            base_line = [line for line in ipfile.split('\n') if key in line]
            if len(base_line) != 1:
                error_message = 'Input creation error when building from basecase. ' +\
                                str(len(base_line)) + ' instances of: ' + key
                raise RuntimeError(error_message)
            ipfile = ipfile.replace(base_line[0], key + ' = ' + str(value))

        #TODO: make necessary change if screening

        # Write out the file
        with open(ip_path, 'w') as openfile: # bad naming here
            openfile.write(ipfile)

        # Read it in database
        new_lib = Library(self.paths.database_path, op_path, ip_path, lib_number, screening)
        self.slibs.append(new_lib) if screening else self.flibs.append(new_lib)
        if screening:
            self.update_metrics(screening=screening)
        else:
            self.update_metrics(screening=screening)
        return

    # Runs exploration and exploitation to build the database
    def build(self, exploration_count, exploitation_count, print_progress=False):
        self.update_metrics()
        self.print()

        # Add some initial points
        self.initial_exploration(False)
        self.run_pxsgen(False)

        # Perform exploration
        if print_progress:
            print('\n_____________________________________\n-- Exploration Step. Total points:', exploration_count)
        for i in range(exploration_count):
            if print_progress:
                print('Generating point', len(self.flibs))
            self.exploration(False)
            self.run_pxsgen(False)
            if print_progress:
                print('  Estimating error of point')
            self.estimate_error()
            #self.find_error(method='cubic')

        # Perform exploitation
        if print_progress:
            print('\n_____________________________________\n-- Exploitation Step. Total points:', exploitation_count)
        for i in range(exploitation_count):
            if print_progress:
                print('Generating point (exploitation)', len(self.flibs))
            self.exploitation()
            self.run_pxsgen(False)
            if print_progress:
                print('  Estimating error of point')
            self.estimate_error()
            #self.find_error(method='cubic')
        self.find_error(print_result=True)
        self.find_error(method='cubic', print_result=True)

        # Write errors
        print(self.paths.database_path, 'complete')
        ip_path = self.paths.database_path + 'errors.txt'
        with open(ip_path, 'w') as openfile:  # bad naming here
            openfile.write('max errors\n' + str(self.est_error_max).replace(',', ''))
            openfile.write('\nmin errors\n' + str(self.est_error_min).replace(',', ''))
            openfile.write('\nmean errors\n' + str(self.est_error_mean).replace(',', ''))
            openfile.write('\nreal errors\n' + str(self.database_error).replace(',', ''))

    # Estimates the error of the database using leave-1-out method
    def estimate_error(self, method='linear', save_result=True, print_result=False):
        # Skip if points are too few
        if len(self.flibs) < self.dimensions * 2:
            return

        # TODO: screening check
        lib_errors = []
        for i in range(len(self.flibs)):
            interpolated = self.interpolate(self.flibs[i].coordinate, method=method, exclude=self.flibs[i].number)
            real = main('', self.flibs[i].coordinate)
            try:
                self.flibs[i].excluded_error = 100 * abs(real - interpolated) / real
                lib_errors.append(self.flibs[i].excluded_error)
            except ZeroDivisionError:
                return
        if save_result:
            self.est_error_mean.append(round(sum(lib_errors)/max(len(lib_errors), 1), 2))
            self.est_error_max.append(round(max(lib_errors), 2))
            self.est_error_min.append(round(min(lib_errors), 2))
        if print_result:
            print('Estimated error:', round(sum(lib_errors)/max(len(lib_errors), 1), 2))

    # Exploitation loop, generates next point based on outputs
    def exploitation(self, print_output=False):
        # Check if there are enough libraries
        if len(self.flibs) < self.dimensions * 2 + 2:
            # Not enough libs to construct neighborhoods
            print('Not enough flibs in database for exploitation')
            return

        """"
        # Check if neighbors are found
        try:
            self.flibs[0].neighborhood
        except AttributeError:
            print('Building initial neighborhoods (this may take a while)')
            self.generate_neighbors()

        # Update neighbors, start by newest point
        try:
            self.flibs[-1].neighborhood
        except AttributeError:
            print('  Building neighborhood of new point')
            considered_libs = self.flibs[-1].proximity_order[:self.dimensions * 4]
            self.update_neighbors(self.flibs[-1], considered_libs=considered_libs)
            # Update the neighbors of new points neighbors
            print('  Updating neighbors of new points neighbors:')
            for i in self.flibs[-1].neighborhood.lib_numbers:
                print('    Updating neighbors of lib', self.flibs[i].number)
                considered_libs = self.flibs[i].proximity_order[:self.dimensions * 4]
                self.update_neighbors(self.flibs[i], considered_libs=considered_libs)
        """
        # Find ranks
        if print_output:
            print('  Finding ranks of database')
        self.generate_ranks()
        if print_output:
            print('    Ranks:')
            for lib in self.flibs:
                print(lib.number, lib.rank)

        # Find the point with highest rank and add it
        max_rank_i = [lib.rank for lib in self.flibs].index(max(lib.rank for lib in self.flibs))
        rounded_point = [round(i, 2) for i in self.flibs[max_rank_i].furthest_point]
        if print_output:
            print('Selected lib', self.flibs[max_rank_i].number, 'point:', rounded_point)
        self.add_lib(self.flibs[max_rank_i].furthest_point, False)

    # Finds the coordinates of next point to sample
    def exploration(self, screening=False):
        # Get the coordinates of the current database
        if screening:
            coords = [i.coordinates(self.varied_ips) for i in self.slibs]
        else:
            coords = [i.coordinates(self.varied_ips) for i in self.flibs]
        if len(coords) == 0:
            print('Exploration - no points in database error')
            return

        # Iterate through all random points
        rand_count = len(coords) * 5000 #TODO: make this better, should probably depend on dimensions too
        rand_points = [[random.random() for i in range(self.dimensions)] for i in range(rand_count)]
        #print('first rand point:', rand_points[0], 'len of rands:', len(rand_points))

        p_cand = [3,3]		# Candidate point to be selected next
        maximin = 0			# The maximin distance of the selected point (higher better)
        fail_count = 0		# Number of rejected rand points
        for counter, rand in enumerate(rand_points):
            #print('checking point', rand)
            projection_fail = False		# I know, this is a n00b way to iterate...FINE #TODO: have better flow control
            min_tot = 10
            for p in coords:
                tot_dist = 0

                for d in range(self.dimensions):
                    dist = (rand[d] - p[d])**2	# Cartesian distance will be calculated with this
                    if dist < self.proj_threshold:	# Projection check
                        projection_fail = True
                        #print('  failed point, dist:', dist)
                    else:
                        tot_dist += dist
                if projection_fail:
                    fail_count += 1
                    break
                tot_dist = (tot_dist)**0.5		# Total cartesian distance of p from rand point
                #print('  total distance:', tot_dist, ' min distance:', min_tot)
                if tot_dist < min_tot:			# Finds the closest distance (in coords) to rand point
                    #print('  assigned new min_tot')
                    min_tot = tot_dist
            if min_tot > maximin and not projection_fail:
                #print('assigned new maximin')
                maximin = min_tot
                p_cand = rand

            # Check if projection threshold is too low
            if (counter+1) % 100 == 0 and fail_count/(counter+1) > 0.50:
                self.proj_threshold /= 2
                # print('Increasing projection threshold to', self.proj_threshold, 'fail rate was at',
                #      round(fail_count/(counter+1), 3))

        # Create a new library with the selected next point and add it to flibs/slibs
        self.add_lib(p_cand, screening)    # Also updates metrics

    # Generates new points for the purpose of finding database error
    def find_error(self, method='linear', save_result=True, print_result=False, multiplier=15000):
        # Skip if points are too few
        if len(self.flibs) < 6:
            return

        # Generate random points for database
        rand_count = multiplier * self.dimensions
        values = copy.deepcopy(self.basecase.inputs.xsgen)
        rand_points = [copy.copy(values) for i in range(rand_count)]

        for i in range(rand_count):
            for key in self.varied_ips:
                rand_points[i][key] = random.random()

        # Iterate through points to find error of each
        tot_error = 0
        max_error = 0
        min_error = 100
        for rand in rand_points:
            rand_varied = []
            real = main('', [x, y])
            point_error = 100 * abs(real - interpolated) / real
            if point_error > max_error:
                max_error = point_error
            if point_error < min_error:
                min_error = point_error
            tot_error += point_error
        tot_error /= rand_count
        if save_result:
            self.database_error.append(round(tot_error, 2))
        if print_result:
            print('Real error:', round(tot_error, 2), '%. Max error:', round(max_error, 2), '%. Min error:',
                  round(min_error, 2), '%')

    # Generates all neighborhoods after exploration
    def generate_neighbors(self):
        self.update_proximity()
        # Go through all full libraries and generate initial neighborhood
        for i, lib in enumerate(self.flibs):
            #TODO: will need to update the initial neighborhood guess

            # Check if there are enough libraries
            if len(self.flibs) < self.dimensions * 2 + 2:
                # Not enough libs to construct neighborhoods
                print('Not enough flibs in database for generate_neighbors')
                return

            initial_neighbors = list(range(self.dimensions*2))
            # Avoid passing the lib in its own neighborhood
            if i in initial_neighbors:
                initial_neighbors[i] = self.dimensions*2

            neighbor_libs = [self.flibs[i] for i in initial_neighbors]
            neighbor_coordinates = [lib.coordinates(self.varied_ips) for lib in neighbor_libs]
            lib.neighborhood = Neighborhood(lib.coordinates(self.varied_ips), initial_neighbors, neighbor_coordinates)

            # Go through all combinations and pick best neighborhood for each point
            print(' Building neighbors of point', i, 'of', len(self.flibs)-1)
            considered_libs = lib.proximity_order[:self.dimensions * 3]
            self.update_neighbors(lib, considered_libs=considered_libs)

    # Finds the gradient of all flibs
    def generate_ranks(self):
        # Estimate voronoi cell sizes
        self.voronoi()      #TODO: in the future this will be optimized

        """"
        # Go through all flibs
        for lib in self.flibs:
            # Update outputs
            self.neighbor_outputs(lib)

            # Find nonlinearity score
            lib.neighborhood.calculate_nonlinearity()
        """
        self.estimate_error(save_result=False)

        # Go through all libs again now that nonlinearity scores are found
        total_nonlinearity = sum([lib.excluded_error for lib in self.flibs])
        for lib in self.flibs:
            # Calculate rank
            lib.rank = 3 * lib.voronoi_size + lib.excluded_error / total_nonlinearity
            # print(lib.number, lib.voronoi_size, lib.excluded_error / total_nonlinearity, lib.rank)

    # Creates the initial set of inputs before exploration begins
    def initial_exploration(self, screening):
        #TODO: this will cause the code to fail if there are only 1 or 2 points and the rest get skipped
        # Make sure this is really initial
        points = len(self.slibs) if screening else len(self.flibs)
        if points > 0:
            # print('Initial exploration method called, but there are already libraries in database')
            return

        # Assign all dimensions 0, 0.5, and 1 to create 3 initial input libs
        for val in [0, 0.5, 1]:
            coords = [val for i in range(self.dimensions)]
            self.add_lib(coords, screening)

        # Update metrics
        self.update_metrics()

    # Finds the interpolated output value at location using method
    def interpolate(self, location, method='linear', exclude=None):
        # Available methods: ‘linear’ or ‘cubic’
        # Database metrics should be updated before running

        data_matrix = copy.copy(self.lib_inputs)
        outputs = copy.copy(self.lib_outputs)

        if exclude is not None:
            data_matrix = data_matrix[:exclude] + data_matrix[exclude + 1:]
            outputs = outputs[:exclude] + outputs[exclude + 1:]

        data_matrix = np.asarray(data_matrix)
        outputs = np.asarray(outputs)
        location = np.asarray(location)

        interpolated = griddata(data_matrix, outputs, location, method=method).tolist()[0]

        if np.isnan(interpolated):
            distances = []
            for point in data_matrix:
                distances.append(distance.euclidean(point, location))
            interpolated = outputs[distances.index(min(distances))]

        return interpolated

    # Places the output data to the neighborhood of lib
    def neighbor_outputs(self, lib):
        outputs = []
        for i in lib.neighborhood.lib_numbers:
            #TODO: in the future this will access combined output
            outputs.append(self.flibs[i].max_BU + self.flibs[i].max_prod + self.flibs[i].max_dest)
        lib.neighborhood.outputs = outputs
        lib.neighborhood.p_output = lib.max_BU + lib.max_prod + lib.max_dest

    #TODO: PCA method need update
    def PCA(self):
        if len(self.max_prods) > 0:
            self.np_prods = np.asarray(self.max_prods)
            self.np_dests = np.asarray(self.max_dests)
            self.np_BUs = np.asarray(self.max_BUs)

            self.data_mat = np.column_stack((self.np_prods, self.np_dests, self.np_BUs))
            self.pca_mat = mlabPCA(self.data_mat)   # PCA matrix

    # Plots data
    def plot(self):
        # Plot input points
        x = [i[0] for i in self.lib_inputs]
        y = [i[1] for i in self.lib_inputs]
        s = [10 for i in self.lib_inputs]   # Sizes
        s[-1] = 50

        fig, ax = plt.subplots()
        ax.scatter(x, y, s=s)
        ax.set_xlim([-0.1, 1.1])
        ax.set_ylim([-0.1, 1.1])
        for i in range(len(x)):
            ax.annotate(str(i), (x[i], y[i]))

        # Plot underlying function
        grid_x, grid_y = np.mgrid[0:1:100j, 0:1:100j]
        plt.imshow(main('', [grid_x, grid_y]), extent=(0, 1, 0, 1), origin='lower')

        plt.show()

    # Prints information about database
    def print(self):
        print('Database ', self.paths.database_path)
        print(' Screening ips:', len(self.slibs), ' Full ips:', len(self.flibs), '  Dimensions:', self.dimensions)

    # Randomly selects and adds the next point
    def random_next(self, screening=False):
        self.add_lib([random.random() for i in range(self.dimensions)], screening)

    #TODO: basecase output and general numbering will need to be thought-out
    def read_base(self, path):
        database_path = self.paths.database_path
        op_path = database_path + self.paths.FR_Output_folder + self.paths.base_output
        self.basecase = Library(database_path, op_path, path, 0, False)

        # Also save this file to write new inputs later
        self.base_file = None
        with open(path, 'r') as bfile:
            self.base_file = bfile.read()

    def read_input(self, ip_path):
        # Make sure it's there
        if not os.path.exists(ip_path):
            error_message = 'The database input file does not exist. Looking for: ' \
                            + ip_path
            raise RuntimeError(error_message)

        doc = open(ip_path, "r")

        for line in doc.readlines():
            items = line.split()

            # Check empty line
            if len(items) == 0:
                continue

            # Check for comments
            if items[0][0] == '#' or items[0][0] == '/':
                continue

            # Check for a 6 digit nucid and that a range is given
            try:
                if len(str(int(items[0]))) == 6 and len(items) == 2:
                    self.ip_ranges.initial_heavy_metal[int(items[0])] = float(items[1])
                    self.varied_ips.append(int(items[0]))
            except ValueError:
                pass

            # Check the rest of (most of) defined inputs
            if len(items) == 2 and float(items[1]) > 0:
                if items[0] in self.ip_ranges.xsgen:
                    self.ip_ranges.xsgen[items[0]] = float(items[1])
                    self.varied_ips.append(items[0])
                if items[0] in self.inputs:	# Check database inputs
                    self.inputs[items[0]] = float(items[1])

        self.dimensions = len(self.varied_ips)

    # Runs pxsgen on all waiting inputs and adds result to database
    def run_pxsgen(self, screening):
        if screening:
            for i in range(len(self.slibs)):
                shell_arg = self.paths.pxsgen_command + ' ' + self.slibs[i].ip_path + ' ' + self.slibs[i].op_path
                if not os.path.exists(self.slibs[i].op_path):
                    subprocess.run(shell_arg, shell=True)
                    self.slibs[i].read_output(self.slibs[i].op_path, 1)
        else:
            for i in range(len(self.flibs)):
                shell_arg = self.paths.pxsgen_command + ' ' + self.flibs[i].ip_path + ' ' + self.flibs[i].op_path
                if not os.path.exists(self.flibs[i].op_path):
                    subprocess.run(shell_arg, shell=True)
                    self.flibs[i].read_output(self.flibs[i].op_path, 1)
        self.update_metrics()

    # A catch-all to update data
    def update_metrics(self, screening=None, libs=None):
        if len(self.slibs) == 0 and len(self.flibs) == 0:
            return

        # Figure out what to update
        if libs is not None:
            pass
        else:
            if screening is None:
                libs = self.slibs + self.flibs
            else:
                libs = self.slibs if screening else self.flibs

        # Update normalized values
        for lib in libs:
            for ip in self.basecase.inputs.xsgen.keys():
                #TODO: assigning normalized values will need to be fixed for xsgen
                lib.normalized[ip] = lib.inputs.xsgen[ip]
                lib.coordinate = lib.coordinates(self.varied_ips)

        # Update library proximity values
        self.update_proximity()

        # Update database coordinate matrix and output vector
        if not screening:
            self.lib_inputs = []
            self.lib_outputs = []
            for flib in self.flibs:
                self.lib_inputs.append(flib.coordinate)
                self.lib_outputs.append(flib.max_BU + flib.max_prod + flib.max_dest)

    # Updates the neighborhoods of flib in database
    def update_neighbors(self, flib, considered_libs=None):
        if len(self.flibs) < self.dimensions * 2 + 2:
            # Not enough libs to construct neighborhoods
            print('Not enough flibs in database for update_neighbors')
            return

        try:
            flib.neighborhood
        except AttributeError:  #TODO: this may not be necessary since it's done in generate_neighbors
            # Stupidly guess initial neighbors
            initial_neighbors = list(range(self.dimensions * 2))
            # Avoid passing the lib in its own neighborhood
            if flib.number in initial_neighbors:
                initial_neighbors[flib.number] = self.dimensions * 2
            neighbor_libs = [self.flibs[i] for i in initial_neighbors]
            neighbor_coordinates = [lib.coordinates(self.varied_ips) for lib in neighbor_libs]
            flib.neighborhood = Neighborhood(flib.coordinates(self.varied_ips), initial_neighbors, neighbor_coordinates)

        # Save current neighborhood information
        current_score = 0
        current_coords = flib.neighborhood.p_coords
        current_neighborhood = flib.neighborhood.lib_numbers

        # Create a list of all candidate libraries
        if considered_libs is None:
            cand_list = list(range(len(self.flibs)))
            cand_list.remove(flib.number)
        else:
            cand_list = considered_libs

        # Iterate through each combination
        for subset in itertools.combinations(cand_list, self.dimensions*2):
            n_coordinates = []
            for i in subset:
                n_coordinates.append(self.flibs[i].coordinates(self.varied_ips))
            new_neighborhood = Neighborhood(current_coords, subset, n_coordinates)

            # Update the current data if new has better score
            if current_score < new_neighborhood.neighbor_score:
                #TODO: this could be improved, information already saved in new_neighborhood
                current_score = new_neighborhood.neighbor_score
                # print(current_score, subset, round(100*current_iter/total_iters),'%')
                current_coords = new_neighborhood.p_coords
                current_neighborhood = copy.deepcopy(new_neighborhood)

        flib.neighborhood = copy.deepcopy(current_neighborhood)

    # Updates the proximity among of points in the database
    def update_proximity(self):
        # For each target point in the database find the distance from other points
        for lib_target in self.flibs:
            distances_target = {}
            for lib_dbase in self.flibs:
                lib_distance = distance.euclidean(lib_target.coordinate, lib_dbase.coordinate)
                distances_target[lib_dbase.number] = lib_distance
            # Order the distances, removing first one since it's the library itself
            lib_target.proximity_order = sorted(distances_target, key=distances_target.get)[1:]

    # Finds the estimate of voronoi cell sizes in database
    def voronoi(self, s_mult=500):
        # For the set of input points in d dimensional space,
        # generates samples number of random points. For each random
        # point, finds which point in p_coords is closest to it for
        # Voronoi cell volume approximation.
        samples = len(self.flibs) * s_mult
        p_coords = [i.coordinates(self.varied_ips) for i in self.flibs]
        p_vol = [0] * len(self.flibs)
        dimensions = len(p_coords[0])

        # Reset furthest points within Voronoi cells for each lib
        for lib in self.flibs:
            # Setting the furthest point to the point itself effectively resets it: any point will be further
            lib.furthest_point = lib.coordinates(self.varied_ips)
            lib.furthest_point_dist = 0

        # Create samples number of random coordinates
        s_coords = [[random.random() for i in range(dimensions)] for i in range(samples)]

        for s in s_coords:      # Random point, s
            min_dist = 9999
            p_closest = 0       # The point in the database closest to the given random point
            for p in p_coords:  # Point in database, p
                # Save the index and its distance if its closest
                distance_s = distance.euclidean(p, s)
                if min_dist > distance_s:
                    min_dist = distance_s
                    p_closest = p_coords.index(p)
            p_vol[p_closest] += 1.0
            # Update furthest point
            if min_dist > self.flibs[p_closest].furthest_point_dist:
                self.flibs[p_closest].furthest_point_dist = min_dist
                self.flibs[p_closest].furthest_point = s
                #print('-updated', self.flibs[p_closest].number, 'to', round(s[0], 2), round(s[1], 2), 'dist:', round(min_dist,2))

        p_vol = [i/samples for i in p_vol]

        self.voronoi_sizes = p_vol
        for i, lib in enumerate(self.flibs):
            lib.voronoi_size = p_vol[i]
