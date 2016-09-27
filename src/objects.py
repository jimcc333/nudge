import os
import copy
import math
import itertools
import random

import numpy as np
from matplotlib.mlab import PCA as mlabPCA
from operator import attrgetter
from scipy.spatial import distance

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class PathNaming:
	# A class that holds all info about the naming of system file paths
	def __init__(self, database_path = "/home/cem/nudge/db_dbtest1/"):
		self.database_path = database_path
		self.xsgen_command = 'xsgen --rc'
		
		self.base_input = 'basecase.py'		# This file is used as the base case for the database
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
	#	Neighborhoods are used to calculate gradients of points
	#	A library does not need outputs to have a fully defined neighborhood
	
	def __init__(self, p_coords, lib_numbers, coordinates):
		self.p_coords = p_coords
		self.lib_numbers = lib_numbers
		self.cohesion = 1	# C=1 implies all points are as far away as possible
		self.adhesion = 0	# A=1 implies all points are on the same spot
		self.coordinates = coordinates
		self.CalculateScore()
	
	def CalculateScore(self):
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
		self.neighbor_score = self.adhesion / ( (self.cohesion**2)*math.sqrt(2))
	
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
		
		# --- Normalized Values ---
		self.normalized = {
			'fuel_cell_radius': None,
			'void_thickness': None,
			'clad_thickness': None,
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
		
		#TODO: pass combining fractions (frac) better
		if os.path.isdir(op_path):
			self.completed = True
			
			u235_file = op_path + "/922350.txt"
			self.ReadOutput("U235", u235_file, 0.04)
			
			u238_file = op_path + "/922380.txt"
			self.ReadOutput("U235", u238_file, 0.96)
			
			#print("Completed reading scout output #" + str(number))
		else:
			self.completed = False
			
			#TODO: read full run output
			#TODO: add library to queue
	
	
	def ReadOutput(self, nuclide, file_path, frac):	
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
				self.max_BU += sum( [float(i) for i in items[1:]] ) * frac
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

	def Coordinates(self, varied_ips):
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
		slib_neighbors, flib_neighbors
		
	
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
		print('reading basecase')
		
		# Read basecase input
		self.ReadBase(paths.database_path + paths.base_input)
		print('read basecase')
		
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
		
		# If there are files SR_Inputs, read them
		if tot_sfiles > 0:
			for ip_number in range(tot_sfiles):
					ip_path = paths.database_path + paths.SR_Input_folder + '/' + str(ip_number) +'.py'
					op_path = paths.database_path + paths.SR_Output_folder + '/' + paths.xsgen_prefix \
								+ paths.sr_prefix + str(ip_number) + '/' + paths.xsgen_op_folder
								
					if os.path.exists(ip_path):
						inputlib = Library(database_path=paths.database_path, op_path=op_path, ip_path=ip_path, number=ip_number, scout=True)
						#TODO: fix this!!! 
						self.slibs.append(inputlib)
						tot_sr_libs += 1
					else:
						# Could continue here instead of breaking, but at this point this is better
						break
		else:
			# If there are no inputs, the first one is basecase
			print(self.basefile)
			
			
		#TODO: also read full input libs
		
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
		op_path = database_path + self.paths.FR_Output_folder + '/basecase'
		self.basecase = Library(database_path, op_path, path, 0, False)
		
		# Also save this file to write new inputs later
		self.basefile = None
		with open(path, 'r') as bfile:
			self.basefile = bfile.read()	
	
	# Once normalized coords are generated, pass here to add
	#TODO: should check if exists first
	#TODO: should update relevant neighborhoods
	#TODO: outputting in correct units
	def AddLib(self, norm_coords, screening):
		# Adding a library will:
		#	- Convert normalized coords to correct units
		#	- Generate the input file
		#	- Add to flib or slib list in this object
		print('begin add lib')
		# Dimension consistency check (note that this is techically not sufficient)
		if len(norm_coords) != self.dimensions:
			error_message = 'Dimensions mismatch when adding a new library to database'
			raise RuntimeError(error_message)
		
		# Convert normalized to correct unit value
		new_inputs = self.basecase.inputs.xsgen
		print(new_inputs)
		for key, value in norm_coords.items():
			print(key, value)
			new_inputs[key] = value #* self.basecase.inputs.xsgen[key]
		
		# Check if input exists
		for lib in (self.slibs if screening else self.flibs):
			if lib.inputs.xsgen == new_inputs:
				print('Input already exists')
				return
		
		# Create paths
		lib_number = len(self.slibs)
		source_path = self.paths.database_path + self.paths.base_input
		if screening:
			ip_path = self.paths.database_path + self.paths.SR_Input_folder +\
					'/' + str(lib_number) + '.py'
			op_path = self.paths.database_path + self.paths.SR_Output_folder +\
					'/' + str(lib_number) + '.py'
		else:
			ip_path = self.paths.database_path + self.paths.FR_Input_folder +\
					'/' + str(lib_number) + '.py'
			op_path = self.paths.database_path + self.paths.FR_Output_folder +\
					'/' + str(lib_number) + '.py'
		print(source_path, ip_path)
		
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
		
		print('end add lib')
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
		
	def UpdateMetrics(self):
		if len(self.slibs) == 0 and len(self.flibs) == 0:
			return
		# work in progress
		
		# Delete old data
		self.max_prods.clear()
		self.max_dests.clear()
		self.max_BUs.clear()
			
		# Rebuild lib values
		self.complete_slibs = 0
		for i in self.slibs:
			if i.completed:
				self.max_prods.append(i.max_prod)
				self.max_dests.append(i.max_dest)
				self.max_BUs.append(i.max_BU)
				
				self.complete_slibs += 1
		
		#TODO: have a class to handle metrics
		# Update the range of metrics
		self.range_fuel_radius[0] = min([i.inputs.xsgen['fuel_cell_radius'] for i in self.slibs])
		self.range_fuel_radius[1] = max([i.inputs.xsgen['fuel_cell_radius'] for i in self.slibs])
		
		self.range_fuel_density[0] = min([i.inputs.xsgen['fuel_density'] for i in self.slibs])
		self.range_fuel_density[1] = max([i.inputs.xsgen['fuel_density'] for i in self.slibs])
		
		self.range_clad_density[0] = min([i.inputs.xsgen['clad_density'] for i in self.slibs])
		self.range_clad_density[1] = max([i.inputs.xsgen['clad_density'] for i in self.slibs])
		
		self.range_cool_density[0] = min([i.inputs.xsgen['cool_density'] for i in self.slibs])
		self.range_cool_density[1] = max([i.inputs.xsgen['cool_density'] for i in self.slibs])
		
		self.range_enrichment[0] = min([i.inputs.xsgen['enrichment'] for i in self.slibs])
		self.range_enrichment[1] = max([i.inputs.xsgen['enrichment'] for i in self.slibs])
		
		# Update the metrics in libraries
		for lib in self.slibs:
			lib.normalized['fuel_cell_radius'] = (lib.inputs.xsgen['fuel_cell_radius'] - self.range_fuel_radius[0]) / \
									(self.range_fuel_radius[1] - self.range_fuel_radius[0])
									
			lib.normalized['fuel_density'] = (lib.inputs.xsgen['fuel_density'] - self.range_fuel_density[0]) / \
									(self.range_fuel_density[1] - self.range_fuel_density[0])
									
			lib.normalized['clad_density'] = (lib.inputs.xsgen['clad_density'] - self.range_clad_density[0]) / \
									(self.range_clad_density[1] - self.range_clad_density[0])
									
			lib.normalized['cool_density'] = (lib.inputs.xsgen['cool_density'] - self.range_cool_density[0]) / \
									(self.range_cool_density[1] - self.range_cool_density[0])
									
			lib.normalized['enrichment'] = (lib.inputs.xsgen['enrichment'] - self.range_enrichment[0]) / \
									(self.range_enrichment[1] - self.range_enrichment[0])
	
	""" 
	#TODO: needs to be rewritten for compatibility
	# Uses inverse distance weighing to find library at target (t_) metrics		
	def EstLib(self, neighbor_libs, alpha=0.5, t_fuel_radius=-1, t_fuel_density=-1, t_clad_density=-1, \
				t_cool_density=-1, t_enrichment=-1):
		# Notes: 
		# 	- The passed variables need to be normalized [0,1]
		print('Begining workflow to interpolate a library')
		
		# Read parameters and store them in metrics lists
		t_metrics = []
		lib_metrics =[]
		lib_metrics_names = []
		
		if t_fuel_radius >= 0 and t_fuel_radius <= 1:
			metrics = []
			for lib in neighbor_libs:
				metrics.append(lib.inputs.norm_fuel_radius)
			lib_metrics.append(metrics)
			lib_metrics_names.append('Norm Fuel Rad:')
			t_metrics.append(t_fuel_radius)
		
		if t_fuel_density >= 0 and t_fuel_density <= 1:
			metrics = []
			for lib in neighbor_libs:
				metrics.append(lib.inputs.norm_fuel_density)
			lib_metrics.append(metrics)
			lib_metrics_names.append('Norm Fuel Dens:')
			t_metrics.append(t_fuel_density)
		
		if t_clad_density >= 0 and t_clad_density <= 1:
			metrics = []
			for lib in neighbor_libs:
				metrics.append(lib.inputs.norm_clad_density)
			lib_metrics.append(metrics)
			lib_metrics_names.append('Norm Clad Dens:')
			t_metrics.append(t_clad_density)
		
		if t_cool_density >= 0 and t_cool_density <= 1:
			metrics = []
			for lib in neighbor_libs:
				metrics.append(lib.inputs.norm_cool_density)
			lib_metrics.append(metrics)
			lib_metrics_names.append('Norm Cool Dens:')
			t_metrics.append(t_cool_density)
				
		if t_enrichment >= 0 and t_enrichment <= 1:
			metrics = []
			for lib in neighbor_libs:
				metrics.append(lib.inputs.norm_enrichment)
			lib_metrics.append(metrics)
			lib_metrics_names.append('Norm Enrichment:')
			t_metrics.append(t_enrichment)
		
		if len(lib_metrics) < 0:
			print('Error, no parameters for interpolation.')
			return
		if len(lib_metrics[0]) < 1:
			print('Error, not enough libraries for interpolation')
			return
		
		print(' Parameters')
		for i in range(len(lib_metrics_names)):
			print('  ', lib_metrics_names[i], t_metrics[i])
		print(' Libraries for interpolation: ', len(lib_metrics[0]))
		
		# Calculate distances
		lib_distances = []
		
		for lib_i in range(len(lib_metrics[0])):
			distance = 1
			for met_i in range(len(lib_metrics)):
				met_dist = 1-(lib_metrics[met_i][lib_i] - t_metrics[met_i])**2
				#TODO: make this a global treshold variable or return the matching lib
				if met_dist == 0:
					met_dist = 0.0001**2
				distance *= met_dist
			distance = distance**(alpha/2)
			lib_distances.append(distance)		
		tot_dist = sum(lib_distances)
		
		# Interpolate library
		interpolated_lib = Library("x","x",-1)
		for lib_i, lib in enumerate(neighbor_libs):
			interpolated_lib.max_prod += lib.max_prod * lib_distances[lib_i] / tot_dist
			interpolated_lib.max_dest += lib.max_dest * lib_distances[lib_i] / tot_dist
			interpolated_lib.max_BU   += lib.max_BU * lib_distances[lib_i] / tot_dist
		
		interpolated_lib.Print()
		print(' Library interpolation complete')
		return interpolated_lib
		"""

	
	# INCOMPLETE
	# Creates the initial set of inputs before exploration begins
	def InitialExploration(self, screening):
		# Make sure this is really initial
		points = len(self.slibs) if screening else len(self.flibs)
		if points > 0:
			return
		
		
	
	# Finds the coordinates of next point to sample
	def Exploration(self, screening = False):
		# p_count: number of points to find
		
		print('exploration begins')
		# Get the coordinates of the current database
		if screening:
			coords = [i.Coordinates(self.varied_ips) for i in self.slibs]
		else:
			coords = [i.Coordinates(self.varied_ips) for i in self.flibs]
		if len(coords) == 0:
			#TODO: add basecase lib to database appropriately, this is probably done before this is called
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
		
		print('exploration ends')
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
						
			
	def EstVoronoi(self, s_mult = 500):
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
		samples = len(self.slibs) * s_mult
		p_coords = [i.Coordinates(self.varied_ips) for i in self.slibs]
		p_vol = [0 for i in range(len(self.slibs))]
		dimensions = len(p_coords[0])
		
		# Create samples number of random coordinates
		#random.seed(1) #TODO: remove this line eventually
		s_coords = 	[[random.random() for i in range(dimensions)] for i in range(samples)]
		
		for s in s_coords:
			min_dist = 9999
			p_closest = p_coords[0]
			for p in p_coords:
			# Save the index and its distance if its closest
				x1 = [value for value in p.values()]
				if min_dist > distance.euclidean(x1,s):
					min_dist = distance.euclidean(x1,s)
					p_closest = p_coords.index(p)
			p_vol[p_closest] += 1.0
		
		p_vol = [i/samples for i in p_vol]
		for i in p_vol: print(round(i,6), end=' ')
		return p_vol
		
	def UpdateNeigbors(self, slib=None):
		
		if len(self.slibs) < self.dimensions * 2 + 2:
			# Not enough libs to construct neighborhoods
			return
		
		if slib != None:
			# Save current neighborhood information
			current_score = self.slib_neighbors[slib].neighbor_score
			current_n_neighbors = self.slib_neighbors[slib].lib_numbers
			current_coords = self.slibs[slib].Coordinates(self.varied_ips)
			current_neighborhood = self.slib_neighbors[slib]
			
			# Create a list of all candidate libraries
			list1 = list(range(len(self.slibs)))
			list1.remove(slib)
			# Iterate through each combination
			total_iters = len(itertools.combinations(list1, self.dimensions*2))
			current_iter = 0
			for subset in itertools.combinations(list1, self.dimensions*2):
				current_iter += 1
				n_coordinates = []
				for i in subset:
					n_coordinates.append(self.slibs[i].Coordinates(self.varied_ips))
				new_neighborhood = Neighborhood(current_coords, subset, n_coordinates)
				# Update the current data if new has better score
				if current_score < new_neighborhood.neighbor_score:
					current_score = new_neighborhood.neighbor_score
					print(current_score, subset, 100*current_iter/total_iters,'%')
					current_n_neighbors = subset
					current_coords = new_neighborhood.p_coords
					current_neighborhood = new_neighborhood
			self.slib_neighbors[slib] = current_neighborhood
			return
		
		# Go through each slib and update its neighborhood
		for lib in self.slibs:
			#print('Updating the neighborhood of lib', lib.number)
			
			# Lib doesn't have a defined neighborhood
			if len(self.slib_neighbors) <= lib.number:
				#print(' Constructing an initial neighborhood for lib', lib.number)
				neighbors = list(range(len(self.slibs)))
				neighbors = neighbors[:lib.number] + neighbors[lib.number+1:]
				# selection could be improved, but I predict it never will be worth it (esp if there's saved state data)
				if len(neighbors) > self.dimensions*2:
					neighbors = neighbors[:self.dimensions*2]
				n_coordinates = []
				for i in neighbors:
					n_coordinates.append(self.slibs[i].Coordinates(self.varied_ips))
				#print('coordinates:', n_coordinates)
				self.slib_neighbors.append(Neighborhood(lib.Coordinates(self.varied_ips), neighbors, n_coordinates))
			
			# Compare the current library score to all other possible options
			
			
			
			#TODO: pass the new libs and only check them
		
		
	def Print(self):
		print('Database screening ips:', len(self.slibs), ' full ips:', len(self.flibs))
		print("Database screening run ranges: ")
		print("  Fuel radius range:  ", self.range_fuel_radius)
		print("  Fuel density range: ", self.range_fuel_density)
		print("  Clad denisty range: ", self.range_clad_density)
		print("  Cool density range: ", self.range_cool_density)
		print("  Enrichment range:   ", self.range_enrichment)
		#print("  Max prods: ", self.max_prods)
		#print("  Max desds: ", self.max_dests)
		#print("  Max BUs  : ", self.max_BUs)
		
		





















