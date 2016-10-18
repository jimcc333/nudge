import os
import copy
import itertools
import random
import subprocess

import numpy as np
from matplotlib.mlab import PCA as mlabPCA
from scipy.spatial import distance

from library import Library
from objects import xsgenParams

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
        self.voronoi_sizes = []         # Voronoi cell sizes of points in the database
        self.np_prods = None            # numpy array of neutron production values
        self.np_dests = None            # numpy array of neutron destruction values
        self.np_BUs = None              # numpy array of burnup values
        self.data_mat = []              # numpy data matrix of neutron prod/dest and BU
        self.pca_mat = None             # numpy PCA matrix
        self.est_error_max = []         # Estimated maximum database error
        self.est_error_min = []         # Estimated minimum database error
        self.est_error_mean = []        # Estimated mean database error
        self.database_error = []        # "True" database error calculated using generated test points

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
        new_inputs = self.basecase.inputs.xsgen
        for key, value in norm_coords.items():
            new_inputs[key] = value #* self.basecase.inputs.xsgen[key]

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
            self.update_metrics(screening=screening, libs=[self.slibs[-1]])
        else:
            self.update_metrics(screening=screening, libs=[self.flibs[-1]])
        return

    # Estimates the error of the database using leave-1-out method
    def estimate_error(self):
        # TODO: screening check
        lib_errors = []
        for i in range(len(self.flibs)):
            int_lib = self.interpolate_lib(self.flibs[:i] + self.flibs[i + 1:], self.flibs[i].normalized)
            try:
                lib_errors.append(100 * abs(self.flibs[i].max_BU - int_lib.max_BU) / self.flibs[i].max_BU)
            except ZeroDivisionError:
                return
        print('Estimated average error:', round(sum(lib_errors)/max(len(lib_errors), 1), 2), '%')
        print('Estimated maximum error:', round(max(lib_errors), 2), '%')
        print('Estimated minimum error:', round(min(lib_errors), 2), '%')
        self.est_error_mean.append(round(sum(lib_errors)/max(len(lib_errors), 1), 2))
        self.est_error_max.append(round(max(lib_errors), 2))
        self.est_error_min.append(round(min(lib_errors), 2))

    # Exploitation loop, generates next point based on outputs
    def exploitation(self):
        # Check if neighbors are found
        try:
            self.flibs[0].neighborhood
        except AttributeError:
            print('Building initial neighborhoods (this may take a while)')
            self.generate_neighbors()

            print('y:', self.flibs[0].neighborhood)

        print('y:', self.flibs[0].neighborhood)
        # Update neighbors, start by newest point
        try:
            self.flibs[-1].neighborhood
        except AttributeError:
            print(' Building neighborhood of new point')
            self.update_neighbors(self.flibs[-1])
            # Update the neighbors of new points neighbors
            for i in self.flibs[-1].neighborhood.lib_numbers:
                print(' Updating neighbors of new points neighbors')
                self.update_neighbors(self.flibs[i])

        # Find ranks
        print(' generating ranks')
        self.generate_ranks()
        # Find the point with highest rank and add it
        max_rank_i = [lib.rank for lib in self.flibs].index(max(lib.rank for lib in self.flibs))
        next_point = {}
        for i, key in enumerate(self.varied_ips):
            next_point[key] = self.flibs[max_rank_i].furthest_point[i]
        print('adding next point')
        self.add_lib(next_point, False)

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
        rand_count = len(coords) * 100 #TODO: make this better, should probably depend on dimensions too
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
                p_list = list(p.values())		# Python3 guarantees the order will be consistent if dict not altered

                for d in range(self.dimensions):
                    dist = (rand[d] - p_list[d])**2	# Cartesian distance will be calculated with this
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
                print('Increasing projection threshold to', self.proj_threshold, 'fail rate was at',
                      round(fail_count/(counter+1), 3))

        #print('Selected:', p_cand, ' rejected:', 100*round(fail_count/len(rand_points),3),'%')

        # Create a new library with the selected next point and add it to flibs/slibs
        # Convert p_cand to dict
        cand_coords = {}
        for counter, key in enumerate(coords[0].keys()):
            cand_coords[key] = p_cand[counter]

        self.add_lib(cand_coords, screening)

    # Generates new points for the purpose of finding database error
    def find_error(self):
        # Generate random points for database
        rand_count = 3000 * self.dimensions
        values = copy.deepcopy(self.flibs[0].inputs.xsgen)
        rand_points = [copy.copy(values) for i in range(rand_count)]

        for i in range(rand_count):
            for key in self.varied_ips:
                rand_points[i][key] = random.random()

        # Iterate through points to find error of each
        tot_error = 0
        for rand in rand_points:
            rand_varied = {}
            for key in self.varied_ips:
                rand_varied[key] = rand[key]
            int_lib = self.interpolate_lib(self.flibs, rand_varied)
            lib_BU = int_lib.max_BU
            real_BU = burnup_maker(rand)
            tot_error += 100 * abs(real_BU - lib_BU) / real_BU
        tot_error /= rand_count
        print('Real error:', round(tot_error, 2), '%')
        self.database_error.append(round(tot_error, 2))

    # Generates all neighborhoods after exploration
    def generate_neighbors(self):
        # Go through all full libraries and generate initial neighborhood
        for i, lib in enumerate(self.flibs):
            #TODO: will need to update the initial neighborhood guess
            initial_neighbors = list(range(self.dimensions*2))
            # Avoid passing the lib in its own neighborhood
            if i in initial_neighbors:
                initial_neighbors[i] = self.dimensions*2
            neighbor_libs = [self.flibs[i] for i in initial_neighbors]
            neighbor_coordinates = [lib.coordinates(self.varied_ips) for lib in neighbor_libs]
            lib.neighborhood = Neighborhood(lib.coordinates(self.varied_ips), initial_neighbors, neighbor_coordinates)

            print('z:', lib.neighborhood)
            # Go through all combinations and pick best neighborhood for each point
            print(' Building neighbors of point', i, 'of', len(self.flibs)-1)
            self.update_neighbors(lib)

        print('y:', self.flibs[0].neighborhood)

    # Finds the gradient of all flibs
    def generate_ranks(self):
        print('x:', self.flibs[0].neighborhood)
        # Estimate voronoi cell sizes
        self.voronoi()      #TODO: in the future this will be optimized
        print('end voronoi')
        # Go through all flibs
        for lib in self.flibs:
            # Update outputs
            print(lib.neighborhood)
            self.neighbor_outputs(lib)

            # Find nonlinearity score
            lib.neighborhood.calculate_nonlinearity()

        # Go through all libs again now that nonlinearity scores are found
        total_nonlinearity = sum([lib.neighborhood.nonlinearity for lib in self.flibs])
        for lib in self.flibs:
            # Calculate rank
            lib.rank = lib.voronoi_size + lib.neighborhood.nonlinearity / total_nonlinearity

    # Creates the initial set of inputs before exploration begins
    def initial_exploration(self, screening):
        # Make sure this is really initial
        points = len(self.slibs) if screening else len(self.flibs)
        if points > 0:
            return

        # Assign all dimensions 0, 0.5, and 1 to create 3 input libs
        coords = {}
        for val in [0, 0.5, 1]:
            for ip in self.varied_ips:
                coords[ip] = val
            self.add_lib(coords, screening)

        # Update metrics
        self.update_metrics()

    # Uses inverse distance weighing to find library at target (t_) metrics
    def interpolate_lib(self, neighbor_libs, norm_coords, alpha=500):
        # Notes:
        # - The passed variables need to be normalized [0,1]
        if len(neighbor_libs) == 1:
            print('Warning! Only 1 library passed for interpolation')
            return neighbor_libs[0]

        # Read parameters and store them in metrics lists
        t_coords = {}
        for key, value in norm_coords.items():
            if value is not None:
                t_coords[key] = value

        # Calculate distances from each neighbor_lib
        lib_distances = []

        for lib in neighbor_libs:
            distance_lib = 1
            for key, value in t_coords.items():
                dim_dist = 1 - (value - lib.normalized[key]) ** 2
                if dim_dist == 0:  # TODO: fix this to make it global var
                    dim_dist = 0.0001 ** 2
                distance_lib *= dim_dist
            distance_lib **= (alpha / 2)
            lib_distances.append(distance_lib)
        tot_dist = sum(lib_distances)

        # Generate interpolated library object
        interpolated_lib = copy.deepcopy(neighbor_libs[0])
        interpolated_lib.reset()
        for key, value in t_coords.items():
            interpolated_lib.normalized[key] = value

        for lib_i, lib in enumerate(neighbor_libs):
            interpolated_lib.max_prod += lib.max_prod * lib_distances[lib_i] / tot_dist
            interpolated_lib.max_dest += lib.max_dest * lib_distances[lib_i] / tot_dist
            interpolated_lib.max_BU += lib.max_BU * lib_distances[lib_i] / tot_dist

        return interpolated_lib

    # Places the output data to the neighborhood of lib
    def neighbor_outputs(self, lib):
        outputs = []
        print('lib')
        print(lib.neighborhood)
        print('lib numbers:', lib.neighborhood.lib_numbers)
        for i in lib.neighborhood.lib_numbers:
            #TODO: in the future this will access combined output
            outputs.append(self.flibs[i].max_BU)
        lib.neighborhood.outputs = outputs
        lib.neighborhood.p_output = lib.max_BU

    def PCA(self):
        if len(self.max_prods) > 0:
            self.np_prods = np.asarray(self.max_prods)
            self.np_dests = np.asarray(self.max_dests)
            self.np_BUs = np.asarray(self.max_BUs)

            self.data_mat = np.column_stack((self.np_prods, self.np_dests, self.np_BUs))
            self.pca_mat = mlabPCA(self.data_mat)   # PCA matrix

    # Prints information about database
    def print(self):
        print('Database screening ips:', len(self.slibs), ' full ips:', len(self.flibs))
        print('  Dimensions:', self.dimensions)

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

            # Check the rest of defined inputs
            if len(items) == 2 and float(items[1]) > 0:
                if items[0] in self.ip_ranges.xsgen:
                    self.ip_ranges.xsgen[items[0]] = float(items[1])
                    self.varied_ips.append(items[0])
                if items[0] in self.inputs:	# Check database inputs
                    self.inputs[items[0]] = float(items[1])

        self.dimensions = len(self.varied_ips)

    # Runs pxsgen on all waiting inputs and adds result to database
    def run_pxsgen(self, exploitation):
        if exploitation:
            for i in range(len(self.slibs)):
                shell_arg = self.paths.pxsgen_command + ' ' + self.slibs[i].ip_path + ' ' + self.slibs[i].op_path
                if not os.path.exists(self.slibs[i].op_path):
                    subprocess.run(shell_arg, shell=True)
                    self.slibs[i].read_output(0, self.slibs[i].op_path, 1)
        else:
            for i in range(len(self.flibs)):
                shell_arg = self.paths.pxsgen_command + ' ' + self.flibs[i].ip_path + ' ' + self.flibs[i].op_path
                if not os.path.exists(self.flibs[i].op_path):
                    subprocess.run(shell_arg, shell=True)
                    self.flibs[i].read_output(0, self.flibs[i].op_path, 1)

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
            for ip in self.varied_ips:
                lib.normalized[ip] = lib.inputs.xsgen[ip]

    # Updates the neighborhoods of flib in database
    def update_neighbors(self, flib):
        print('update neighbors')

        print('a:', flib.neighborhood)
        if len(self.flibs) < self.dimensions * 2 + 2:
            # Not enough libs to construct neighborhoods
            print('not enough flibs in database for update_neighbors')
            return

        try:
            flib.neighborhood
        except AttributeError:
            initial_neighbors = list(range(self.dimensions * 2))
            # Avoid passing the lib in its own neighborhood
            if flib.number in initial_neighbors:
                initial_neighbors[flib.number] = self.dimensions * 2
            neighbor_libs = [self.flibs[i] for i in initial_neighbors]
            neighbor_coordinates = [lib.coordinates(self.varied_ips) for lib in neighbor_libs]
            print('aa:', flib.neighborhood)
            flib.neighborhood = Neighborhood(flib.coordinates(self.varied_ips), initial_neighbors, neighbor_coordinates)

        print('b:', self.flibs[0].neighborhood)
        # Save current neighborhood information
        current_score = flib.neighborhood.neighbor_score
        current_coords = flib.neighborhood.p_coords
        current_neighborhood = flib.neighborhood.lib_numbers

        # Create a list of all candidate libraries
        cand_list = list(range(len(self.flibs)))
        cand_list.remove(flib.number)

        # Iterate through each combination
        current_iter = 0
        for subset in itertools.combinations(cand_list, self.dimensions*2):
            current_iter += 1
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
                current_neighborhood = new_neighborhood
        print('cur neigh',current_neighborhood)
        flib.neighborhood = current_neighborhood

    # Finds the estimate of voronoi cell sizes in database
    def voronoi(self, s_mult = 500):
        """
        The algorithm should work in the following way:

        ep_vols = list(0, len(existing_points))			# Voronoi cell vol of each ep (existing point)
        ep_cands = list(None, len(existing_points))		# The best candidate for next input point near ep

        for rp in random_points:
            min_dist = 1
            for ep in existing_points:
                if dist(rp,ep) < min_dist:
                    min_dist = dist(rp,ep)
                    rp_cell = ep						# rp is in the Voronoi cell of ep
            ep_vols[rp_cell] += 1/len(random_points)	# this rp contributes to the volume of rp_cell
            if dist(rp_cell,rp) > dist(rp_cell, ep_cands[rp_cell]) or ep_cands[rp_cell] == None:
                if passes_other_crit(rp):				# Catch-all, but mainly projective property check
                    ep_cands[rp_cell] = rp
        """
        # For the set of input points in d dimensional space,
        # generates samples number of random points. For each random
        # point, finds which point in p_coords is closest to it for
        # Voronoi cell volume approximation.
        # Uses the info to estimate voronoi cell volume of points in p_coords.
        samples = len(self.flibs) * s_mult
        p_coords = [i.coordinates(self.varied_ips) for i in self.flibs]
        p_vol = [0] * len(self.flibs)
        dimensions = len(p_coords[0])

        # Reset furthest points within Voronoi cells for each lib
        for lib in self.flibs:
            # Setting the furthest point to the point itself effectively resets it: any point will be further
            lib.furthest_point = lib.coordinates(self.varied_ips)

        # Create samples number of random coordinates
        #random.seed(1) #TODO: remove this line eventually
        s_coords = [[random.random() for i in range(dimensions)] for i in range(samples)]

        for s in s_coords:      # Random point, s
            min_dist = 9999
            p_closest = 0       # The point in the database closest to the given random point
            for p in p_coords:  # Point in database, p
                # Save the index and its distance if its closest
                x1 = [value for value in p.values()]
                distance_s = distance.euclidean(x1, s)
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
