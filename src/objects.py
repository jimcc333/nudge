import math
import random

import numpy as np
from scipy.spatial import distance


class PathNaming:
    # A class that holds all info about the naming of system file paths
    def __init__(self, os_name, database_path="/home/cem/nudge/db_dbtest1/"):
        if os_name == 'nt': self.win = True
        else: self.win = False

        self.slash = '/'
        if self.win:
            self.slash = '\\'

        self.database_path = database_path
        self.xsgen_command = 'xsgen --rc'
        self.pxsgen_command = 'python3 src/pxsgen.py'
        if self.win:
            self.pxsgen_command = 'python src\\pxsgen.py'

        self.base_input = 'basecase.py'     # This file is used as the base case for the database
        self.base_output = self.slash + 'basecase.py'
        self.dbase_input = 'inputs.txt'     # Used for database creation inputs

        self.SR_Input_folder = 'SR_Inputs'      # Where scout run inputs are stored
        self.SR_Output_folder = 'SR_Outputs'    # Where scout run outputs are stored
        self.FR_Input_folder = 'FR_Inputs'      # Where full run inputs are stored
        self.FR_Output_folder = 'FR_Outputs'    # Where full run outputs are stored

        self.xsgen_prefix = 'build-'            # The prefix xsgen assigns to output folders
        self.xsgen_op_folder = 'brightlite0'    # The folder xsgen places bright-lite formatted outputs

        self.sr_prefix = 'sr'                   # The prefix given to scout libs
        self.fr_prefix = 'fr'                   # The prefix given to full libs

        self.database_name = 'database1'


class xsgenParams:
    # A class that holds all parameters also used in xsgen inputs
    def __init__(self):
        # Initial heavy metal mass fraction distribution
        self.initial_heavy_metal = {}

        self.xsgen = {
            # Geometry inputs
            'fuel_cell_radius': None,   # [cm]
            'void_cell_radius': None,   # [cm]
            'clad_cell_radius': None,   # [cm]
            'unit_cell_pitch': None,    # [cm]
            'unit_cell_height': None,   # [cm]
            # Density inputs
            'fuel_density': None,       # Fuel density [g/cc]
            'clad_density': None,       # Cladding Density [g/cc]
            'cool_density': None,       # Coolant Density [g/cc]
            # Other inputs
            'enrichment': None,         # Fuel enrichment (Uranium fuel only) as atom fraction
            'flux': None,               # Average reactor flux [n/cm2/s]
            'k_particles': None,        # Number of particles to run per kcode cycle
        }

    # Returns the number of defined inputs
    def defined_count(self):
        defined_inputs = 0
        defined_inputs += sum(1 for i in self.initial_heavy_metal.values() if i is not None)
        defined_inputs += sum(1 for i in self.xsgen.values() if i is not None)
        return defined_inputs

    def print_defined(self):
        for key, value in self.initial_heavy_metal.items():
            if value is not None:
                print(key, value)
        for key, value in self.xsgen.items():
            if value is not None:
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

    # Calculates the neighborhood score using adhesion and cohesion criteria
    def calculate_score(self):
        # Check if center is in its own neighbors
        if self.p_coords in self.coordinates:
            error_message = 'Point given in its own neighborhood for point ' + str(self.p_coords)
            raise RuntimeError(error_message)

        # Cohesion: average(distance(coord, center))
        distances = []
        for coord in self.coordinates:
            distances.append(distance.euclidean(self.p_coords, coord))
        self.cohesion = np.mean(distances)

        # Adhesion: average(min(distance(coord1,coord2)))
        tot_libs = len(self.lib_numbers)
        distances.clear()
        lib_distances = []
        for lib1 in range(tot_libs):
            for lib2 in range(tot_libs):
                if lib1 != lib2:
                    lib_distances.append(distance.euclidean(self.coordinates[lib1], self.coordinates[lib2]))
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
        for neighbor_coords in self.coordinates:
            row = []
            for i, value in enumerate(neighbor_coords):
                row.append(value - self.p_coords[i])    # Coordinates of center subtracted to make it in origin
            A.append(row)

        # Solve the linear equation
        matrix_A = np.array(A)
        vector_b = np.array(self.outputs)
        gradient = np.linalg.lstsq(matrix_A, vector_b.transpose())[0]

        # Find nonlinearity score using gradient
        self.nonlinearity = 0
        for i, point in enumerate(A):
            self.nonlinearity += abs(self.p_output - vector_b[i] + np.dot(gradient, point))
