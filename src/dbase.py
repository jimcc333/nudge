import os
import copy
import itertools
import random
import subprocess

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.mlab import PCA as mlabPCA
from scipy.spatial import distance
from scipy.interpolate import griddata

from library import Library
from objects import xsgenParams, Neighborhood, PathNaming
from pxsgen import main

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

    def __init__(self, database_path):
        # Read database, assuming software may have been interrupted and
        # the folder may have some inputs and outputs

        # Database constants
        self.paths = PathNaming(os.name, database_path=database_path)
        self.varied_ips = []            # the varied inputs
        self.ip_ranges = xsgenParams()  # Ranges of varied inputs
        self.dimensions = None          # Number of varied inputs for database (assigned when inputs are read)
        self.basecase = None            # The basecase Library object (assigned when basecase is read)
        self.base_file = None           # The basecase as string (assigned when basecase is read)

        # Database inputs
        self.inputs = {
            'max_exploration': 0,       # Total number of new exploration points to add to the database
            'max_exploitation': 0,      # Total number of new exploitation points to add to the database
            'max_samples': None,        # Maximum database size (flibs)
            'max_error': None,          # In [%]
            'max_time': 100,            # In [hour]
            'scout_frac': 10,           # Weight of screening time allocation
            'explore_frac': 40,         # Weight of exploration time allocation
            'exploit_frac': 50,         # Weight of exploitation time allocation
            'explore_mult': 500,        # Exploration method Monte Carlo multiplier
            'max_projection': 0.0001,   # Projection check (exploration) threshold
            'voronoi_mult': 200,        # Voronoi method Monte Carlo multiplier
            'rank_factor': 1,           # The factor that multiplies error when finding rank
            'voronoi_adjuster': 0.8,    # The maximum ratio of voronoi cell adjustment (guided method) [0,1]
            'guide_increment': 0.0001,  # The increment to bring back selected guided sample back to original V cell
        }

        # Database libraries
        self.slibs = []                 # Screening libraries (with and without outputs)
        self.flibs = []                 # Full libraries (with and without outputs)

        # Database analysis
        self.lib_inputs = []            # (d,n) matrix of varied input normalized coordinates (n: number of flibs)
        self.lib_outputs = []           # (n) vector of library outputs
        self.voronoi_sizes = []         # Voronoi cell sizes of points in the database
        self.distance_factors = []      # Distance adjustment factors for Voronoi cell calculation
        self.np_prods = None            # numpy array of neutron production values
        self.np_dests = None            # numpy array of neutron destruction values
        self.np_BUs = None              # numpy array of burnup values
        self.data_mat = []              # numpy data matrix of neutron prod/dest and BU
        self.pca_mat = None             # numpy PCA matrix
        self.est_error_max = []         # Estimated maximum database error vector
        self.est_error_min = []         # Estimated minimum database error vector
        self.est_error_mean = []        # Estimated mean database error vector
        self.database_error = []        # "True" database error vector calculated using generated test points
        self.database_error_max = []    # "True" database max error vector calculated using generated test points

        # Check that database path exists
        if not os.path.isdir(self.paths.database_path):
            error_message = 'The database path does not exist. Looking for: ' + self.paths.database_path
            raise NotADirectoryError(error_message)

        # Read database input file
        self.read_input(self.paths.database_path + self.paths.dbase_input)

        # Read basecase input
        self.read_base(self.paths.database_path + self.paths.base_input)

        # Check to see if there's a screening library input folder
        if os.path.exists(self.paths.database_path + self.paths.SR_Input_folder):
            tot_sfiles = len(os.listdir(self.paths.database_path + self.paths.SR_Input_folder))
        # If the input screening library folder doesn't exist create it
        else:
            os.mkdir(self.paths.database_path + self.paths.SR_Input_folder)
            tot_sfiles = 0

        # Check to see if there's a full library input folder
        if os.path.exists(self.paths.database_path + self.paths.FR_Input_folder):
            tot_ffiles = len(os.listdir(self.paths.database_path + self.paths.FR_Input_folder))
        # If the input full library folder doesn't exist create it
        else:
            os.mkdir(self.paths.database_path + self.paths.FR_Input_folder)
            tot_ffiles = 0

        tot_sr_libs = 0
        tot_fr_libs = 0

        # If there are files SR_Inputs, read them
        if tot_sfiles > 0:
            for ip_number in range(tot_sfiles):
                    ip_path = self.paths.database_path + self.paths.SR_Input_folder + self.paths.slash + str(ip_number) + '.py'
                    # TODO: op_path will need to be updated to work with xsgen
                    op_path = self.paths.database_path + self.paths.SR_Output_folder + self.paths.slash + str(ip_number) + '.py'
                    #+ self.paths.xsgen_prefix + self.paths.sr_prefix + str(ip_number) + '/' + self.paths.xsgen_op_folder

                    if os.path.exists(ip_path):
                        inputlib = Library(database_path=self.paths.database_path, op_path=op_path, ip_path=ip_path,
                                           number=ip_number, scout=True)
                        self.slibs.append(copy.deepcopy(inputlib))
                        tot_sr_libs += 1
                    else:
                        # Stops reading inputs if they're not sequential. This will overwrite inputs
                        break

        # If there are files FR_Inputs, read them
        if tot_ffiles > 0:
            for ip_number in range(tot_ffiles):
                ip_path = self.paths.database_path + self.paths.FR_Input_folder + self.paths.slash + str(ip_number) + '.py'
                # TODO: op_path will need to be updated to work with xsgen
                op_path = self.paths.database_path + self.paths.FR_Output_folder + self.paths.slash + str(ip_number) + '.py'
                # + self.paths.xsgen_prefix + self.paths.sr_prefix + str(ip_number) + '/' + self.paths.xsgen_op_folder

                if os.path.exists(ip_path):
                    inputlib = Library(database_path=self.paths.database_path, op_path=op_path, ip_path=ip_path,
                                       number=ip_number, scout=True)
                    self.flibs.append(copy.deepcopy(inputlib))
                    tot_fr_libs += 1
                else:
                    # Could continue here instead of breaking, but at this point this is better
                    break

        # If there's no output folder, create it
        if not os.path.exists(self.paths.database_path + self.paths.SR_Output_folder):
            os.mkdir(self.paths.database_path + self.paths.SR_Output_folder)
        if not os.path.exists(self.paths.database_path + self.paths.FR_Output_folder):
            os.mkdir(self.paths.database_path + self.paths.FR_Output_folder)

        self.update_coordinates()

    # Once normalized coords are generated, pass here to add
    # TODO: outputting in correct units
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

        # TODO: make necessary change if screening

        # Write out the file
        with open(ip_path, 'w') as openfile: # bad naming here
            openfile.write(ipfile)
        openfile.close()

        # Read it in database
        new_lib = Library(self.paths.database_path, op_path, ip_path, lib_number, screening)
        self.slibs.append(new_lib) if screening else self.flibs.append(new_lib)
        if screening:
            self.update_metrics()
        else:
            self.update_metrics()

    # Runs exploration and exploitation to build the database from input file
    def build(self, exploration_to_add=0, exploitation_to_add=0, print_progress=False, record_errors=True,
              exploit_method='furthest'):
        self.update_metrics()
        self.print()

        # Decide how many of each
        lib_count = len(self.flibs)
        exploration_to_add = \
            max(self.inputs['max_exploration'] - lib_count, 0) if exploration_to_add == 0 else exploration_to_add
        exploitation_to_add = \
            max(self.inputs['max_exploitation'] - lib_count, 0) if exploitation_to_add == 0 else exploitation_to_add

        print('Adding', exploration_to_add, 'exploration and', exploitation_to_add, 'exploitation samples')

        # Add some initial points
        if lib_count < 3:
            print('\n_____________________________________')
            print('-- Initial Building Step. Samples to add: 3')
            self.initial_exploration(False)
            self.run_pxsgen(False)
            lib_count = len(self.flibs)
            exploration_to_add -= 3

        # Perform exploration
        if print_progress:
            print('\n_____________________________________')
            print('-- Exploration Step. Samples to add:', exploration_to_add)
        for i in range(exploration_to_add):
            if print_progress:
                print('Generating exploration sample', len(self.flibs))
            self.exploration(False)
            self.run_pxsgen(False)
            if record_errors:
                if print_progress:
                    print('  Estimating errors')
                self.estimate_error()
                self.find_error()

        # Perform exploitation
        if print_progress:
            print('\n_____________________________________')
            print('-- Exploitation Step. Total points:', exploitation_to_add)
        for i in range(exploitation_to_add):
            if print_progress:
                print('Generating exploitation sample', len(self.flibs), 'method:', exploit_method)
            self.exploitation(method=exploit_method)
            self.run_pxsgen(False)
            if record_errors:
                if print_progress:
                    print('  Estimating error of point')
                self.estimate_error()
                self.find_error()

        # Write errors
        print(self.paths.database_path, 'complete')
        if record_errors:
            self.write_errors()

    # Calculates the distance weighing factors for Voronoi cell calculation
    def calculate_factors(self, base_point_i):
        selected_error = self.flibs[base_point_i].excluded_error
        zeroed_errors = [lib.excluded_error / selected_error - 1 for lib in self.flibs]
        zeroed_errors[:] = [-0.5 if i < -0.5 else i for i in zeroed_errors]
        adjuster = max([abs(i) for i in zeroed_errors])
        if adjuster > self.inputs['voronoi_adjuster']:
            adjuster = self.inputs['voronoi_adjuster'] / adjuster
        self.distance_factors = [1 + i * adjuster for i in zeroed_errors]

    # Estimates the error of the database using leave-1-out method
    def estimate_error(self, method='linear', save_result=True, print_result=False, exclude_after=None, plot=False):
        # Skip if points are too few
        if len(self.flibs) < self.dimensions * 2 or len(self.flibs) < 10:
            return

        # Handle database exclusions
        exclude_list = []
        if exclude_after is not None:
            exclude_list = list(range(len(self.flibs)))[exclude_after:]

        # TODO: screening check
        lib_errors = []
        for i in range(len(self.flibs)):
            all_excluded = [self.flibs[i].number] + exclude_list
            interpolated = self.interpolate(self.flibs[i].coordinate, method=method, exclude=all_excluded)
            real = main('', self.flibs[i].inputs.xsgen)
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
            print('Estimated errors:', lib_errors)

        if plot:
            errors = [flib.excluded_error for flib in self.flibs]
            x = [i[0] for i in self.lib_inputs]
            y = [i[1] for i in self.lib_inputs]
            fig, ax = plt.subplots()
            ax.set_ylim([-0.01, 1.01])
            ax.set_xlim([-0.01, 1.01])
            ax.scatter(x, y, s=160, c=errors)
            plt.show()

    # Exploitation loop, generates next point based on outputs
    def exploitation(self, print_output=False, method='furthest'):
        # methods: furthest, guided(computationally more expensive)
        # Check if there are enough libraries
        if len(self.flibs) < self.dimensions * 2 + 2:
            # Not enough libs to construct neighborhoods
            print('Not enough flibs in database for exploitation')
            return

        if method not in ['furthest', 'guided']:
            raise RuntimeError('The exploitation method ', method, ' is not defined')

        # Find ranks
        if print_output:
            print('  Finding ranks of database')
        self.generate_ranks()
        if print_output:
            print('    Ranks:')
            for lib in self.flibs:
                print(lib.number, lib.rank)

        # Find the next point
        ranks = [lib.rank for lib in self.flibs]
        max_rank_i = ranks.index(max(ranks))    # The next point is selected near this point (both methods)
        # Find the point with highest rank and add it
        selected_point = self.flibs[max_rank_i].furthest_point
        # print('ranks:', [round(i.rank, 2) for i in self.flibs])

        if method == 'guided':
            # Find adjustment factors
            self.calculate_factors(max_rank_i)
            # Find adjusted voronoi cells
            self.voronoi(factors=self.distance_factors)
            base = self.flibs[max_rank_i].coordinate
            furthest = copy.copy(self.flibs[max_rank_i].furthest_point)
            adjusted_point = furthest

            # Adjust coordinates of selected point so that its in the original voronoi cell
            #   Find normalized vector from furthest to base point
            normal_vector = []
            for d in range(self.dimensions):
                normal_vector.append(furthest[d] - base[d])
            total = sum([abs(i) for i in normal_vector])
            total = 1 if total == 0 else total
            normal_vector[:] = [value/total for value in normal_vector]
            closest_to_base = False

            # print('Picked point', max_rank_i)
            # print('Furthest:', [round(i, 3) for i in furthest])
            # print('  closest to:', self.find_closest(furthest))
            # If it's closest to base, then it's too close: move it away by flipping normal_vector and condition
            exit_condition = True
            if self.find_closest(adjusted_point) == max_rank_i:
                normal_vector[:] = [i for i in normal_vector]
                exit_condition = False
                closest_to_base = True
            # print('exit condition:', exit_condition, '   closest to base:', closest_to_base)
            # print('normal_vector')
            # print(normal_vector)
            while closest_to_base is not exit_condition:
                # Move the point closer to the base point (max_rank_i point)
                adjusted_point = [adjusted_point[i] + normal_vector[i] * self.inputs['guide_increment'] for i in
                                  range(self.dimensions)]
                # Check if adjusting hits the edge
                at_edge = False
                for value in adjusted_point:
                    if value > 1:
                        value = 1
                        at_edge = True
                    if value < 0:
                        value = 0
                        at_edge = True
                if at_edge:
                    break
                # Find which sample the adjusted point is closest
                closest_to_base = True if self.find_closest(adjusted_point) == max_rank_i else False

            # Fix edge issues
            for value in adjusted_point:
                if value > 1:
                    value = 1
                if value < 0:
                    value = 0
            # print('Initial furthest, adjusted:', [round(i, 3) for i in furthest], [round(i, 3) for i in adjusted_point])
            selected_point = adjusted_point

        rounded_point = [round(i, 2) for i in selected_point]
        if print_output:
            print('Selected lib:', self.flibs[max_rank_i].number, ' Coordinates:', rounded_point)
        self.add_lib(selected_point, False)

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
        rand_count = int((len(coords)) ** 0.5) * self.inputs['explore_mult']  # TODO: could depend on dimensions too
        rand_points = [[random.random() for i in range(self.dimensions)] for i in range(rand_count)]
        # print('first rand point:', rand_points[0], 'len of rands:', len(rand_points))

        p_cand = [3,3]		# Candidate point to be selected next
        maximin = 0			# The maximin distance of the selected point (higher better)
        fail_count = 0		# Number of rejected rand points
        for counter, rand in enumerate(rand_points):
            # print('checking point', rand)
            projection_fail = False		# I know, this is a n00b way to iterate...FINE # TODO: have better flow control
            min_tot = 10
            for p in coords:
                tot_dist = 0

                for d in range(self.dimensions):
                    dist = (rand[d] - p[d])**2	# Cartesian distance will be calculated with this
                    if dist < self.inputs['max_projection']:  # Projection check
                        projection_fail = True
                        # print('  failed point, dist:', dist)
                    else:
                        tot_dist += dist
                if projection_fail:
                    fail_count += 1
                    break
                tot_dist = (tot_dist)**0.5		# Total cartesian distance of p from rand point
                # print('  total distance:', tot_dist, ' min distance:', min_tot)
                if tot_dist < min_tot:			# Finds the closest distance (in coords) to rand point
                    # print('  assigned new min_tot')
                    min_tot = tot_dist
            if min_tot > maximin and not projection_fail:
                # print('assigned new maximin')
                maximin = min_tot
                p_cand = rand

            # Check if projection threshold is too low
            if (counter+1) % 100 == 0 and fail_count/(counter+1) > 0.50:
                self.inputs['max_projection'] /= 2

        # Create a new library with the selected next point and add it to flibs/slibs
        self.add_lib(p_cand, screening)    # Also updates metrics

    # Finds the closest point in self.flibs to the given point and returns the index
    def find_closest(self, point):
        closest_dist = 10
        closest_lib = -1
        for i, lib in enumerate(self.flibs):
            lib_dist = distance.euclidean(point, lib.coordinate)
            if lib_dist < closest_dist:
                closest_dist = lib_dist
                closest_lib = i
        return closest_lib

    # Generates new points for the purpose of finding database error
    def find_error(self, method='linear', save_result=True, print_result=False, multiplier=5000):
        # Skip if points are too few
        if len(self.flibs) < self.dimensions * 2 or len(self.flibs) < 10:
            return

        # Generate random points for database
        rand_count = multiplier  # This prevents scaling issues in very high dimensions
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
            for key in sorted(self.varied_ips):
                rand_varied.append(rand[key])
            interpolated = self.interpolate(rand_varied, method=method)
            real = main('', rand)
            point_error = 100 * abs(real - interpolated) / real
            if point_error > max_error:
                max_error = point_error
            if point_error < min_error:
                min_error = point_error
            tot_error += point_error
        tot_error /= rand_count
        if save_result:
            self.database_error.append(round(tot_error, 2))
            self.database_error_max.append(round(max_error, 2))
        if print_result:
            print('Real error:', round(tot_error, 2), '%. Max error:', round(max_error, 2), '%. Min error:',
                  round(min_error, 2), '%')

    # Finds the gradient of all flibs
    def generate_ranks(self):
        # Estimate voronoi cell sizes
        self.voronoi()      # TODO: in the future this will be optimized
        self.estimate_error(save_result=False)

        # Go through all libs again now that nonlinearity scores are found
        total_nonlinearity = sum([lib.excluded_error for lib in self.flibs])

        for lib in self.flibs:
            # Calculate rank
            lib.rank = lib.voronoi_size + lib.excluded_error / total_nonlinearity * self.inputs['rank_factor']
            # print(lib.number, lib.voronoi_size, lib.excluded_error / total_nonlinearity, lib.rank)

    # Creates the initial set of inputs before exploration begins
    def initial_exploration(self, screening):
        # TODO: this will cause the code to fail if there are only 1 or 2 points and the rest get skipped
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

        # Handle exclusions
        if exclude is not None:
            data_matrix = []
            outputs = []
            for i in range(len(self.lib_outputs)):
                if i not in exclude:
                    data_matrix.append(self.lib_inputs[i])
                    outputs.append(self.lib_outputs[i])
        else:
            data_matrix = copy.copy(self.lib_inputs)
            outputs = copy.copy(self.lib_outputs)

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
            # TODO: in the future this will access combined output
            outputs.append(self.flibs[i].max_BU + self.flibs[i].max_prod + self.flibs[i].max_dest)
        lib.neighborhood.outputs = outputs
        lib.neighborhood.p_output = lib.max_BU + lib.max_prod + lib.max_dest

    # TODO: PCA method need update
    def PCA(self):
        if len(self.max_prods) > 0:
            self.np_prods = np.asarray(self.max_prods)
            self.np_dests = np.asarray(self.max_dests)
            self.np_BUs = np.asarray(self.max_BUs)

            self.data_mat = np.column_stack((self.np_prods, self.np_dests, self.np_BUs))
            self.pca_mat = mlabPCA(self.data_mat)   # PCA matrix

    # Plots data
    def plot(self, numbers=False, points=True, est_errors=False, mark_last=False):
        # Plot input points
        x = [i[0] for i in self.lib_inputs]
        y = [i[1] for i in self.lib_inputs]
        s = [10 for i in self.lib_inputs]   # Sizes

        fig, ax = plt.subplots()
        ax.set_ylim([-0.01, 1.01])
        ax.set_xlim([-0.01, 1.01])
        if points:
            errors = [flib.excluded_error if est_errors else 'g' for flib in self.flibs]
            if mark_last is True:
                errors[-1] = 1 if est_errors else 'y'
            ax.scatter(x, y, s=200, c=errors)
            if numbers:
                for i in range(len(x)):
                    ax.annotate(str(i), (x[i], y[i]))

        # Plot underlying function
        grid_x, grid_y = np.mgrid[0:1:100j, 0:1:100j]
        values = copy.deepcopy(self.basecase.inputs.xsgen)
        outputs = copy.deepcopy(grid_x)
        for x in range(len(grid_x)):
            for y in range(len(grid_x[0])):
                values[self.varied_ips[0]] = grid_x[x, y]
                values[self.varied_ips[1]] = grid_y[x, y]
                outputs[x, y] = main('', values)
        plt.imshow(outputs, extent=(0, 1, 0, 1), origin='lower')

        plt.show()

    # Plots database estimate of blackbox output
    def plot_estimate(self, exclude_after=None, diff=False, abs_max=None, abs_min=None, est_errors=False,
                      mark_last=False, print_range=False):
        # Handle database exclusion
        exclude = None
        if exclude_after is not None:
            exclude = list(range(len(self.flibs)))[exclude_after:]
        # Get used points
        data_x = [i[0] for i in self.lib_inputs[:exclude_after]]
        data_y = [i[1] for i in self.lib_inputs[:exclude_after]]

        # Get estimate data
        grid_x, grid_y = np.mgrid[0:1:80j, 0:1:80j]
        values = copy.deepcopy(grid_x)
        for x in range(len(grid_x)):
            for y in range(len(grid_x[0])):
                values[x, y] = self.interpolate([grid_y[x, y], grid_x[x, y]], exclude=exclude)

        # Check difference plotting
        if diff:
            inputs = copy.deepcopy(self.basecase.inputs.xsgen)
            for x in range(len(values)):
                for y in range(len(values[0])):
                    inputs[self.varied_ips[0]] = grid_x[x, y]
                    inputs[self.varied_ips[1]] = grid_y[x, y]
                    values[x, y] = abs(values[x, y] - main('', inputs))
            if print_range:
                print('Abs max:', values.max())
                print('Abs min:', values.min())

        # Check absolute values for coloring limits
        if abs_max is None:
            abs_max = np.amax(values)
        if abs_min is None:
            abs_min = np.amin(values)

        # Plot data
        fig, ax = plt.subplots()
        errors = [flib.excluded_error if est_errors else 'g' for flib in self.flibs]
        if mark_last:
            errors[-1] = max(errors)*2 if est_errors else 'y'
        ax.scatter(data_x, data_y, s=200, c=errors)
        ax.set_xlim([-0.01, 1.01])
        ax.set_ylim([-0.01, 1.01])
        plt.imshow(values, extent=(0, 1, 0, 1), origin='lower', vmin=abs_min, vmax=abs_max)
        plt.show()
        return

    def plot_voronoi(self, resolution=100, base_point_i=None):
        print('Plotting 2D voronoi cells of the database')
        # Generate a grid and get coords of samples
        grid_x, grid_y = np.mgrid[0:1:(resolution*1j), 0:1:(resolution*1j)]
        colors = np.zeros((resolution, resolution))

        samples_x = [i[0] for i in self.lib_inputs]
        samples_y = [i[1] for i in self.lib_inputs]

        if base_point_i is not None:
            self.calculate_factors(base_point_i)

        for x in range(resolution):
            for y in range(resolution):
                min_dist = 9999
                closest_s = None
                for i in range(len(samples_x)):
                    sample_dist = distance.euclidean((grid_x[x, y], grid_y[x, y]), (samples_x[i], samples_y[i]))
                    if base_point_i is not None:
                        sample_dist *= self.distance_factors[i]
                    if min_dist > sample_dist:
                        min_dist = sample_dist
                        closest_s = i
                colors[y][x] = closest_s

        colors = np.divide(colors, colors.max())
        errors = [flib.excluded_error for flib in self.flibs]
        fig, ax = plt.subplots()
        ax.set_xlim([-0.01, 1.01])
        ax.set_ylim([-0.01, 1.01])
        ax.scatter(samples_x, samples_y, s=100, c=errors)
        plt.imshow(colors, extent=(0, 1, 0, 1), origin='lower', interpolation='hermite')
        plt.show()

    # Prints information about database
    def print(self):
        print('Database ', self.paths.database_path)
        # print(' Screening ips:', len(self.slibs), ' Full ips:', len(self.flibs), '  Dimensions:', self.dimensions)

    # Randomly selects and adds the next point
    def random_selection(self, count=1, screening=False, print_progress=False):
        self.print()

        for i in range(count):
            self.add_lib([random.random() for i in range(self.dimensions)], screening)
            self.run_pxsgen(screening)
            if print_progress:
                print('  Estimating error of database')
            self.estimate_error()
            self.find_error()

        # Write errors
        print(self.paths.database_path, 'complete')
        ip_path = self.paths.database_path + 'errors.txt'
        with open(ip_path, 'w') as openfile:  # bad naming here
            openfile.write('max errors\n' + str(self.est_error_max).replace(',', ''))
            openfile.write('\nmin errors\n' + str(self.est_error_min).replace(',', ''))
            openfile.write('\nmean errors\n' + str(self.est_error_mean).replace(',', ''))
            openfile.write('\nreal max errors\n' + str(self.database_error_max).replace(',', ''))
            openfile.write('\nreal errors\n' + str(self.database_error).replace(',', ''))

    # Reads basecase input
    def read_base(self, path):
        # Check that basecase input exists
        if not os.path.exists(path):
            error_message = 'The database base-case input file does not exist. Looking for: ' + path
            raise FileNotFoundError(error_message)

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
            if len(items) == 2:
                if items[0] in self.ip_ranges.xsgen:
                    self.ip_ranges.xsgen[items[0]] = float(items[1])
                    self.varied_ips.append(items[0])
                # Check database inputs
                if items[0] in ['max_exploration', 'max_exploitation', 'max_samples']:
                    self.inputs[items[0]] = int(items[1])
                    continue
                if items[0] in self.inputs:
                    self.inputs[items[0]] = float(items[1])
        doc.close()
        self.dimensions = len(self.varied_ips)

    # Runs pxsgen on all waiting inputs and adds result to database
    def run_pxsgen(self, screening):
        devnull = open(os.devnull, 'w')
        if screening:
            for i in range(len(self.slibs)):
                shell_arg = self.paths.pxsgen_command + ' ' + self.slibs[i].ip_path + ' ' + self.slibs[i].op_path
                if not os.path.exists(self.slibs[i].op_path):
                    subprocess.run(shell_arg, shell=True, stdout=devnull, stderr=devnull)
                    self.slibs[i].read_output(self.slibs[i].op_path, 1)
        else:
            for i in range(len(self.flibs)):
                shell_arg = self.paths.pxsgen_command + ' ' + self.flibs[i].ip_path + ' ' + self.flibs[i].op_path
                if not os.path.exists(self.flibs[i].op_path):
                    subprocess.run(shell_arg, shell=True, stdout=devnull, stderr=devnull)
                    self.flibs[i].read_output(self.flibs[i].op_path, 1)

        self.update_metrics()

    # Updates flib coordinates
    def update_coordinates(self):
        # Update normalized values
        libs = self.slibs + self.flibs
        for lib in libs:
            for ip in self.basecase.inputs.xsgen.keys():
                # TODO: assigning normalized values will need to be fixed for xsgen
                lib.normalized[ip] = lib.inputs.xsgen[ip]
                lib.coordinate = lib.coordinates(self.varied_ips)

        self.lib_inputs = []
        self.lib_outputs = []
        for flib in self.flibs:
            self.lib_inputs.append(flib.coordinate)
            self.lib_outputs.append(flib.max_BU + flib.max_prod + flib.max_dest)

    # A catch-all to update data
    def update_metrics(self):
        self.update_coordinates()
        self.update_proximity()
        self.estimate_error(save_result=False)

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
    def voronoi(self, factors=None):
        # For the set of input points in d dimensional space,
        # generates samples number of random points. For each random
        # point, finds which point in p_coords is closest to it for
        # Voronoi cell volume approximation.
        samples = int((len(self.flibs)) ** 0.5) * self.inputs['voronoi_mult']
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

        if factors is None:
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
        else:
            if len(factors) != len(self.flibs):
                raise RuntimeError('Size mismatch of factor and library vectors in voronoi cell calculation')
            for s_i, s in enumerate(s_coords):      # Random point, s
                min_dist = 9999
                p_closest = 0       # The point in the database closest to the given random point
                for p_i, p in enumerate(p_coords):  # Point in database, p
                    # Save the index and its distance if its closest
                    distance_s = distance.euclidean(p, s) * factors[p_i]
                    if min_dist > distance_s:
                        min_dist = distance_s
                        p_closest = p_coords.index(p)
                p_vol[p_closest] += 1.0
                # Update furthest point
                if min_dist > self.flibs[p_closest].furthest_point_dist:
                    self.flibs[p_closest].furthest_point_dist = min_dist
                    self.flibs[p_closest].furthest_point = s

        p_vol = [i/samples for i in p_vol]

        self.voronoi_sizes = p_vol
        for i, lib in enumerate(self.flibs):
            lib.voronoi_size = p_vol[i]

    def write_errors(self):
        # This should check if the file exists to prevent loss of data
        ip_path = self.paths.database_path + 'errors.txt'
        with open(ip_path, 'w') as openfile:  # bad naming here
            openfile.write('max errors\n' + str(self.est_error_max).replace(',', ''))
            openfile.write('\nmin errors\n' + str(self.est_error_min).replace(',', ''))
            openfile.write('\nmean errors\n' + str(self.est_error_mean).replace(',', ''))
            openfile.write('\nreal max errors\n' + str(self.database_error_max).replace(',', ''))
            openfile.write('\nreal errors\n' + str(self.database_error).replace(',', ''))
