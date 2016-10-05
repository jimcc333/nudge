import os
import copy
import math
import itertools
import random
import subprocess

import numpy as np
from matplotlib.mlab import PCA as mlabPCA
from scipy.spatial import distance

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class PathNaming:
    # A class that holds all info about the naming of system file paths
    def __init__(self, os_name, database_path="/home/cem/nudge/db_dbtest1/"):
        if os_name == 'nt': self.win = True
        else: self.win = False

        self.slash = '/'
        if self.win: self.slash = '\\'

        self.database_path = database_path
        self.xsgen_command = 'xsgen --rc'
        self.pxsgen_command = 'python3 src/pxsgen.py'
        if self.win:
            self.pxsgen_command = 'python src\\pxsgen.py'

        self.base_input = 'basecase.py'		# This file is used as the base case for the database
        self.base_output = self.slash + 'basecase.py'
        self.dbase_input = 'inputs.txt'		# Used for database creation inputs

        self.SR_Input_folder = 'SR_Inputs'	# Where scout run inputs are stored
        self.SR_Output_folder = 'SR_Outputs'# Where scout run outputs are stored
        self.FR_Input_folder = 'FR_Inputs'	# Where full run inputs are stored
        self.FR_Output_folder = 'FR_Outputs'# Where full run outputs are stored

        self.xsgen_prefix = 'build-'		# The prefix xsgen assigns to output folders
        self.xsgen_op_folder = 'brightlite0'# The folder xsgen places bright-lite formatted outputs

        self.sr_prefix = 'sr'				# The prefix given to scout libs
        self.fr_prefix = 'fr'				# The prefix given to full libs

        self.database_name = 'database1'

class xsgenParams:
    # A class that holds all parameters also used in xsgeninputs
    def __init__(self):
        # Initial heavy metal mass fraction distribution
        self.initial_heavy_metal = {}

        self.xsgen = {
            # Geometry inputs
            'fuel_cell_radius': None,	# [cm]
            'void_cell_radius': None,	# [cm]
            'clad_cell_radius': None,	# [cm]
            'unit_cell_pitch': None,	# [cm]
            'unit_cell_height': None,	# [cm]
            # Density inputs
            'fuel_density': None,		# Fuel density [g/cc]
            'clad_density': None,		# Cladding Density [g/cc]
            'cool_density': None,		# Coolant Density [g/cc]
            # Other inputs
            'enrichment': None,		# Fuel enrichment (Uranium fuel only) as atom fraction
            'flux': None,	  		# Average reactor flux [n/cm2/s]
            'k_particles': None,	# Number of particles to run per kcode cycle
        }

    # Returns the number of defined inputs
    def DefinedCount(self):
        defined_inputs = 0
        defined_inputs += sum(1 for i in self.initial_heavy_metal.values() if i != None)
        defined_inputs += sum(1 for i in self.xsgen.values() if i != None)
        return defined_inputs

    def PrintDefined(self):
        for key, value in self.initial_heavy_metal.items():
            if value != None:
                print(key, value)
        for key, value in self.xsgen.items():
            if value != None:
                print(key, value)


class Neighborhood:
    # A class that holds information about the sample point neighborhood
    # Neighborhoods are used to calculate gradients of points
    # A library does not need outputs to have a fully defined neighborhood

    def __init__(self, p_coords, lib_numbers, coordinates):
        self.p_coords = p_coords			# Coordinates of the center point
        self.lib_numbers = lib_numbers		# flib numbers of libraries forming the neighbors
        self.coordinates = coordinates		# (normalized) coordinates of the neighbors
        self.cohesion = 1					# C=1 implies all points are as far away as possible
        self.adhesion = 0					# A=1 implies all points are on the same spot
        self.neighbor_score = 0				# The neighborhood score
        self.p_output = None                # The output of the center point
        self.outputs = []                   # The outputs of the neighbors (assigned after neighborhood is built)
        self.nonlinearity = 0               # The nonlinearity score of the neighborhood (needs outputs)
        self.calculate_score()

    def calculate_score(self):
        # Check if center is in its own neighbors
        if self.p_coords in self.coordinates:
            error_message = 'Point given in its own neighborhood for point ' + str(self.p_coords)
            raise RuntimeError(error_message)

        # Cohesion: avrg(distance(coord, center))
        distances = []
        p = [value for value in self.p_coords.values()]
        for coord in self.coordinates:
            x = [value for value in coord.values()]
            distances.append(distance.euclidean(p,x))
        self.cohesion = np.mean(distances)

        # Adhesion: avrg(min(distance(coord1,coord2)))
        tot_libs = len(self.lib_numbers)
        distances.clear()
        lib_distances = []
        for lib1 in range(tot_libs):
            for lib2 in range(tot_libs):
                if lib1 != lib2:
                    x1 = [value for value in self.coordinates[lib1].values()]
                    x2 = [value for value in self.coordinates[lib2].values()]
                    lib_distances.append(distance.euclidean(x1,x2))
            distances.append(min(lib_distances))
            lib_distances.clear()
        self.adhesion = np.mean(distances)

        # Neighborhood score: A/(sqrt(2)*C^2)
        self.neighbor_score = self.adhesion / ((self.cohesion**2)*math.sqrt(2))

    # Calculates the gradient based on neighbor outputs
    def calculate_nonlinearity(self):
        # Check if outputs are there (lazily)
        if len(self.outputs) == 0:
            print('No outputs of neighborhood')
            return

        # Generate matrix A
        A = []      # This is the matrix so that the center point is in the origin
        for lib_dict in self.coordinates:
            row = []
            for key, value in lib_dict.items():
                row.append(value-self.p_coords[key])    # Coordinates of center subtracted to make it in origin
            A.append(row)

        # Solve the linear equation
        matrix_A = np.array(A)
        vector_b = np.array(self.outputs)
        gradient = np.linalg.lstsq(matrix_A, vector_b.transpose())[0]

        # Find nonlinearity score using gradient
        self.nonlinearity = 0
        for i, point in enumerate(A):
            self.nonlinearity += abs(self.p_output - vector_b[i] + np.dot(gradient, point))


class Library:
    """ A class that holds library information """
    #TODO "run" routine to talk to xsgen, needs to be parallizable

    def __init__(self, database_path, op_path, ip_path, number, scout):

        self.ip_path = ip_path	# Path of the input file
        self.op_path = op_path	# path to the library folder w/ Bright-lite formatted .txt files in it
        self.number = number	# Unique number of the library
        self.scout = scout
        self.inputs = xsgenParams()

        self.max_prod = 0
        self.max_dest = 0
        self.max_BU = 0

        self.voronoi_size = 0           # The Voronoi cell size of the library
        self.furthest_point = []        # The furthest found point within the Voronoi cell
        self.furthest_point_dist = 0    # The distance of the furthest point
        self.rank = 0                   # The rank of the library, used during exploitation

        # --- Normalized Values ---
        self.normalized = {
            'fuel_cell_radius': None,
            'void_cell_radius': None,
            'clad_cell_radius': None,
            'unit_cell_pitch': None,
            'unit_cell_height': None,
            'fuel_density': None,
            'clad_density': None,
            'cool_density': None,
            'enrichment': None,
            'flux': None
        }

        # Read input
        self.ReadInput(ip_path)

        # Read output if exists
        #TODO: pass combining fractions (frac) better
        if os.path.isdir(op_path):
            self.completed = True

            u235_file = op_path + "/922350.txt"
            self.ReadOutput("U235", u235_file, 0.04)

            u238_file = op_path + "/922380.txt"
            self.ReadOutput("U235", u238_file, 0.96)
        elif os.path.exists(op_path):
            self.completed = True
            self.read_output("pxsgen", op_path, 1)
        else:
            self.completed = False
            #TODO: read full run output
            #TODO: add library to queue

    def read_output(self, nuclide, file_path, frac):
        try:
            doc = open(file_path, "r")
        except IOError:
            print("Could not open ", file_path)
            return

        for line in doc.readlines():
            items = line.split()

            if items[0] == "NEUT_PROD":
                self.max_prod += float(items[len(items)-1]) * frac
            if items[0] == "NEUT_DEST":
                self.max_dest += float(items[len(items)-1]) * frac
            if items[0] == "BUd":
                self.max_BU += sum( [float(i) for i in items[2:]] ) * frac
        doc.close()

    def ReadInput(self, ip_path):
        if ip_path == "x":
            return

        try:
            doc = open(ip_path, "r")
        except IOError:
            print("Could not open ", file_path)
            return

        max_lines = 500
        ip = doc.readlines()

        for line_i in range(len(ip)):
            items = ip[line_i].split()

            if len(items) < 3:
                continue
            if items[0] in self.inputs.xsgen:
                self.inputs.xsgen[items[0]] = float(items[2])
            if items[0] == 'initial_heavy_metal':
                while line_i < max_lines:
                    line_i += 1
                    items = ip[line_i].replace(':',' ').replace(',',' ').split()
                    if len(items) > 2:
                        error_message = 'Input file in ' + self.ip_path + \
                                        ' has formatting error at initial_heavy_metal. Make sure each NUCID is in a new line and close bracket (}) at a new line.'
                        raise RuntimeError(error_message)
                    if '}' not in items[0]:
                        self.inputs.initial_heavy_metal[int(items[0])] = float(items[1])
                    else:
                        # End of initial heavy metal
                        break

        if self.inputs.xsgen['enrichment'] != None:
            if self.inputs.xsgen['enrichment'] - (self.inputs.initial_heavy_metal[922350] \
                /(self.inputs.initial_heavy_metal[922350] + self.inputs.initial_heavy_metal[922380])) > 0.001:
                    error_message = 'Input file in ' + self.ip_path + ' has inconsistency between enrichment and given mass compositions'
                    raise RuntimeError(error_message)

    def Print(self, detail=0):
        if self.ip_path == 'x':
            print('Interpolated library output information: ')
            print('  Max prod: ', self.max_prod)
            print('  Max dest: ', self.max_dest)
            print('  Max BU  : ', self.max_BU)
        else:
            print('Lib #', self.number, ' input information:')
            if(detail):
                print(self.inputs.xsgen)
            print(' Path:', self.ip_path)
            print(' Normalized values')
            print(self.normalized)

    def coordinates(self, varied_ips):
        coordinates = {}
        for key, value in self.normalized.items():
            if key in varied_ips:
                coordinates[key] = value
        return coordinates


class DBase:
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
    slibs = []		# Scout libs
    flibs = []		# Full libs

    complete_slibs = 0
    complete_flibs = 0

    # Library neighborhoods
    slib_neighbors = []
    flib_neighbors = []

    # Database inputs
    inputs = {
        'max_error': None,	# In [%]
        'max_time': 100,	# In [hour]
        'scout_frac': 10,	# Weight of scouting time allocation
        'explore_frac': 40,	# Weight of exploration time allocation
        'exploit_frac': 50,	# Weight of exploitation time allocation
    }

    # Scout lib output parameters
    max_prods = []
    max_dests = []
    max_BUs = []

    # Min and max values ("range") in database
    range_fuel_radius = [1,0]	# [min,max]

    range_fuel_density = [1,0]
    range_clad_density = [1,0]
    range_cool_density = [1,0]

    range_enrichment = [1,0]

    def __init__(self, paths):
        # Read database, assuming software may have been interrupted and
        #	the folder may have some inputs and outputs
        self.paths = paths
        self.voronoi_sizes = []         # Voronoi cell sizes of points in the database

        if not os.path.isdir(paths.database_path):
            error_message = 'The database path does not exist. Looking for: ' \
                             + paths.database_path
            raise RuntimeError(error_message)

        if not os.path.exists(paths.database_path + paths.dbase_input):
            error_message = 'The database input file does not exist. Looking for: ' \
                            + paths.database_path + paths.dbase_input
            raise RuntimeError(error_message)

        if not os.path.exists(paths.database_path + paths.base_input):
            error_message = 'The database base-case input file does not exist. Looking for: ' \
                            + paths.database_path + paths.base_input
            raise RuntimeError(error_message)

        self.name = paths.database_name

        # Read database inputs
        self.ReadInput(paths.database_path + paths.dbase_input)

        # Read basecase input
        self.ReadBase(paths.database_path + paths.base_input)

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
                    ip_path = paths.database_path + paths.SR_Input_folder + paths.slash + str(ip_number) +'.py'
                    #TODO: op_path will need to be updated to work with xsgen
                    op_path = paths.database_path + paths.SR_Output_folder + paths.slash + str(ip_number) + '.py'
                    #+ paths.xsgen_prefix + paths.sr_prefix + str(ip_number) + '/' + paths.xsgen_op_folder

                    if os.path.exists(ip_path):
                        inputlib = Library(database_path=paths.database_path, op_path=op_path, ip_path=ip_path, number=ip_number, scout=True)
                        #TODO: fix this!!!
                        self.slibs.append(inputlib)
                        tot_sr_libs += 1
                    else:
                        # Could continue here instead of breaking, but at this point this is better
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
                    # TODO: fix this!!!
                    self.flibs.append(inputlib)
                    tot_fr_libs += 1
                else:
                    # Could continue here instead of breaking, but at this point this is better
                    break

        # If there's no output folder, create it
        if not os.path.exists(paths.database_path + paths.SR_Output_folder):
            os.mkdir(paths.database_path + paths.SR_Output_folder)
        if not os.path.exists(paths.database_path + paths.FR_Output_folder):
            os.mkdir(paths.database_path + paths.FR_Output_folder)

    def ReadInput(self, ip_path):
        # Make sure it's there
        if not os.path.exists(ip_path):
            error_message = 'The database input file does not exist. Looking for: ' \
                            + ip_path
            raise RuntimeError(error_message)

        self.ip_ranges = xsgenParams()	# Ranges of varied inputs
        self.varied_ips = []			# the varied inputs
        self.proj_threshold = 0.0001	# Projection check (exploration) threshold

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

    #TODO: basecase output and general numbering will need to be thought-out
    def ReadBase(self, path):
        database_path = self.paths.database_path
        op_path = database_path + self.paths.FR_Output_folder + self.paths.base_output
        self.basecase = Library(database_path, op_path, path, 0, False)

        # Also save this file to write new inputs later
        self.basefile = None
        with open(path, 'r') as bfile:
            self.basefile = bfile.read()

    # Once normalized coords are generated, pass here to add
    #TODO: outputting in correct units
    def AddLib(self, norm_coords, screening):
        # Adding a library will:
        #	- Convert normalized coords to correct units
        #	- Generate the input file
        #	- Add to flib or slib list in this object

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
        ipfile = self.basefile

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
            self.UpdateMetrics(screening = screening, libs = [self.slibs[-1]])
        else:
            self.UpdateMetrics(screening = screening, libs = [self.flibs[-1]])
        return

    #TODO: this probably isnt used and aint even right
    def Exists(self, number):
        for lib in self.slibs:
            if lib.number == number:
                return True
        return False

    #TODO: divide to completed and queued

    def PCA(self):
        if len(self.max_prods) > 0:
            self.np_prods = np.asarray(self.max_prods)
            self.np_dests = np.asarray(self.max_dests)
            self.np_BUs = np.asarray(self.max_BUs)

            self.data_mat = np.column_stack((self.np_prods, self.np_dests, self.np_BUs))
            self.pca_mat = mlabPCA(self.data_mat)	# PCA matrix

    # A catch-all to update data
    def UpdateMetrics(self, screening = None, libs = None):
        if len(self.slibs) == 0 and len(self.flibs) == 0:
            return

        # Figure out what to update
        if libs != None:
            pass
        else:
            if screening == None:
                libs = self.slibs + self.flibs
            else:
                libs = self.slibs if screening else self.flibs

        # Update normalized values
        for lib in libs:
            for ip in self.varied_ips:
                lib.normalized[ip] = lib.inputs.xsgen[ip]

    # Uses inverse distance weighing to find library at target (t_) metrics
    def interpolate_lib(self, neighbor_libs, norm_coords, alpha=5):
        # Notes:
        # 	- The passed variables need to be normalized [0,1]
        if len(neighbor_libs) == 1:
            print('Warning! Only 1 library passed for interpolation')
            return neighbor_libs[0]

        # Read parameters and store them in metrics lists
        t_coords = {}
        for key, value in norm_coords.items():
            if value != None:
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
        interpolated_lib.Reset()
        for key, value in t_coords.items():
            interpolated_lib.normalized[key] = value

        for lib_i, lib in enumerate(neighbor_libs):
            interpolated_lib.max_prod += lib.max_prod * lib_distances[lib_i] / tot_dist
            interpolated_lib.max_dest += lib.max_dest * lib_distances[lib_i] / tot_dist
            interpolated_lib.max_BU += lib.max_BU * lib_distances[lib_i] / tot_dist

        return interpolated_lib

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
            self.AddLib(coords, screening)

        # Update metrics
        self.UpdateMetrics()

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
                print('Increasing projection threshold to', self.proj_threshold, \
                        'fail rate was at', round(fail_count/(counter+1), 3))

        #print('Selected:', p_cand, ' rejected:', 100*round(fail_count/len(rand_points),3),'%')

        # Create a new library with the selected next point and add it to flibs/slibs
        #	Convert p_cand to dict
        cand_coords = {}
        for counter, key in enumerate(coords[0].keys()):
            cand_coords[key] = p_cand[counter]

        self.AddLib(cand_coords, True)

        '''
        x = [i[0] for i in exp_coords]
        y = [i[1] for i in exp_coords]

        rx = [i[0] for i in rand_points]
        ry = [i[1] for i in rand_points]

        fig, ax = plt.subplots()
        ax.set_xlim([-0.1,1.1])
        ax.set_ylim([-0.1,1.1])
        sizes = 10
        labels = [i for i in range(len(x))]
        ax.scatter(x, y, s=200, c='b')
        #ax.scatter(rx, ry, s=20, c='g')
        #ax.scatter(p_cand[0], p_cand[1], s=200, c='r')
        ax.grid(True)
        #fig.tight_layout()
        for i, txt in enumerate(labels):
            ax.annotate(txt, (x[i],y[i]), xytext = (x[i]-0.03,y[i]+0.03))
        plt.show()
        '''

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
        # Also for each
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
        s_coords = 	[[random.random() for i in range(dimensions)] for i in range(samples)]

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
            # Go through all combinations and pick best neighborhood for each point
            self.update_neighbors(lib)

    # Updates the neighborhoods of flib in database
    def update_neighbors(self, flib):
        if len(self.flibs) < self.dimensions * 2 + 2:
            # Not enough libs to construct neighborhoods
            return

        # Save current neighborhood information
        current_score = flib.neighborhood.neighbor_score
        current_coords = flib.neighborhood.p_coords
        current_neighborhood = flib.neighborhood.lib_numbers

        # Create a list of all candidate libraries
        cand_list = list(range(len(self.flibs)))
        cand_list.remove(flib.number)

        # Iterate through each combination
        # total_iters = len(list(itertools.combinations(cand_list, self.dimensions*2)))
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
        flib.neighborhood = current_neighborhood
        return

    # Finds the gradient of all flibs
    def generate_ranks(self):
        # Estimate voronoi cell sizes
        self.voronoi()      #TODO: in the future this will be optimized

        # Go through all flibs
        for lib in self.flibs:
            # Update outputs
            self.neighbor_outputs(lib)

            # Find nonlinearity score
            lib.neighborhood.calculate_nonlinearity()

        # Go through all libs again now that nonlinearity scores are found
        total_nonlinearity = sum([lib.neighborhood.nonlinearity for lib in self.flibs])
        for lib in self.flibs:
            # Calculate rank
            lib.rank = lib.voronoi_size + lib.neighborhood.nonlinearity / total_nonlinearity

    # Places the output data to the neighborhood of lib
    def neighbor_outputs(self, lib):
        outputs = []
        for i in lib.neighborhood.lib_numbers:
            #TODO: in the future this will access combined output
            outputs.append(self.flibs[i].max_BU)
        lib.neighborhood.outputs = outputs
        lib.neighborhood.p_output = lib.max_BU

    # Exploitation loop, generates next point based on outputs
    def exploitation(self):
        # Find neighbors
        self.generate_neighbors()
        # Find ranks
        self.generate_ranks()
        # Find the point with highest rank and add it
        max_rank_i = [lib.rank for lib in self.flibs].index(max(lib.rank for lib in self.flibs))
        next_point = {}
        for i, key in enumerate(self.varied_ips):
            next_point[key] = self.flibs[max_rank_i].furthest_point[i]
        self.AddLib(next_point, False)

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

    # Estimates the error of the database using leave-1-out method
    # INCOMLETE
    def estimate_error(self):
        # TODO: screening check
        tot_error = 0
        for i in range(len(self.flibs)):
            int_lib = self.interpolate_lib(self.slibs[:i] + self.slibs[i + 1:], self.slibs[i].normalized)
            tot_error += 100 * abs(self.slibs[i].max_BU - int_lib.max_BU) / self.slibs[i].max_BU
        tot_error /= len(self.slibs)

    # Generates new points for the purpose of finding database error
    def find_error(self):
        # Generate random points for database
        rand_count = 3000 * self.dimensions
        values = copy.deepcopy(self.slibs[0].inputs.xsgen)
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
            int_lib = self.interpolate_lib(self.slibs, rand_varied)
            lib_BU = int_lib.max_BU
            real_BU = burnup_maker(rand)
            tot_error += 100 * abs(real_BU - lib_BU) / real_BU
        tot_error /= rand_count
        print('Real error:', round(tot_error, 2), '%')

    # Prints information about database
    def print(self):
        print('Database screening ips:', len(self.slibs), ' full ips:', len(self.flibs))
        print('  Dimensions:', self.dimensions)
        print("Database screening run ranges: ")
        print("  Fuel radius range:  ", self.range_fuel_radius)
        print("  Fuel density range: ", self.range_fuel_density)
        print("  Clad denisty range: ", self.range_clad_density)
        print("  Cool density range: ", self.range_cool_density)
        print("  Enrichment range:   ", self.range_enrichment)
        #print("  Max prods: ", self.max_prods)
        #print("  Max desds: ", self.max_dests)
        #print("  Max BUs  : ", self.max_BUs)

# copied from pxsgen, delete later
def burnup_maker(inputs):
    random.seed(inputs['enrichment'] * inputs['flux'] * 1245)

    x1 = (2 + inputs['enrichment']) ** 3
    x2 = (2.5635 - inputs['cool_density']) ** 2.12
    x3 = (3 - inputs['clad_density']) ** 2.1234
    x4 = (0.5 + inputs['fuel_density']) ** 3.01
    x5 = np.sin(inputs['fuel_cell_radius'] * 2) * 8 + 2
    x6 = 8 - ((inputs['flux'] - 0.71) * 2) ** 4
    x7 = 3 * inputs['void_cell_radius'] + 2 * inputs['clad_cell_radius'] + random.random()

    random.seed()

    return ((x1 + x2 + x3 + x4 + x5 + x6) ** 2.3) / 100 + x7





















