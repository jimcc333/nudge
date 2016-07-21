import os
import copy
import math
import itertools

import numpy as np
from matplotlib.mlab import PCA as mlabPCA
from operator import attrgetter
from scipy.spatial import distance


class PathNaming:
	# A class that holds all info about the naming of system file paths
	def __init__(self, database_path = "/home/cem/nudge/db_dbtest1/"):
		self.database_path = database_path
		
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
		
		self.geometry = {		# Geometry inputs
			'fuel_cell_radius': None,	# [cm]
			'void_cell_radius': None,	# [cm]
			'clad_cell_radius': None,	# [cm]
			'void_thickness': None,		# [cm]
			'clad_thickness': None,		# [cm]
			'unit_cell_pitch': None,	# [cm]
			'unit_cell_height': None,	# [cm]
			}
		
		self.density = {			# Density inputs
			'fuel_density': None,	# Fuel density [g/cc]
			'clad_density': None,	# Cladding Density [g/cc]
			'cool_density': None,	# Coolant Density [g/cc]
			}

		self.other = {		# Others
			'enrichment': None,		# Fuel enrichment (Uranium fuel only) as atom fraction
			'flux:': None,	  		# Average reactor flux [n/cm2/s]
			'k_particles': None,	# Number of particles to run per kcode cycle
			'k_cycles': None,		# Number of kcode cycles to run
			'k_cycles_skip': None,	# Number of kcode cycles to run but skip
			'group_structure': None,
			'openmc_group_struct': None,
		}
	
	# Returns the number of defined inputs
	def DefinedCount(self):
		defined_inputs = 0
		defined_inputs += sum(1 for i in self.initial_heavy_metal.values() if i != None)
		defined_inputs += sum(1 for i in self.geometry.values() if i != None)
		defined_inputs += sum(1 for i in self.density.values() if i != None)
		defined_inputs += sum(1 for i in self.other.values() if i != None)
		return defined_inputs
	
	def PrintDefined(self):
		for key, value in self.initial_heavy_metal.items():
			if value != None:
				print(key, value)
		for key, value in self.geometry.items():
			if value != None:
				print(key, value)
		for key, value in self.density.items():
			if value != None:
				print(key, value)
		for key, value in self.other.items():
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
			'fuel_density': None,
			'clad_density': None,
			'cool_density': None,
			'enrichment': None,
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
			if items[0] == 'fuel_cell_radius':
				self.inputs.geometry[items[0]] = float(items[2])
			if items[0] == 'void_cell_radius':
				self.inputs.geometry[items[0]] = float(items[2])
			if items[0] == 'clad_cell_radius':
				self.inputs.geometry[items[0]] = float(items[2])
			if items[0] == 'unit_cell_pitch':
				self.inputs.geometry[items[0]] = float(items[2])
			if items[0] == 'unit_cell_height':
				self.inputs.geometry[items[0]] = float(items[2])
			if items[0] == 'fuel_density':
				self.inputs.density[items[0]] = float(items[2])
			if items[0] == 'clad_density':
				self.inputs.density[items[0]] = float(items[2])
			if items[0] == 'cool_density':
				self.inputs.density[items[0]] = float(items[2])
			if items[0] == 'enrichment':
				self.inputs.other[items[0]] = float(items[2])
			if items[0] == 'flux':
				self.inputs.other[items[0]] = float(items[2])
			if items[0] == 'k_particles':
				self.inputs.other[items[0]] = float(items[2])
			if items[0] == 'k_cycles':
				self.inputs.other[items[0]] = float(items[2])
			if items[0] == 'k_cycles_skip':
				self.inputs.other[items[0]] = float(items[2])
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

		if self.inputs.other['enrichment'] != None:
			if self.inputs.other['enrichment'] - (self.inputs.initial_heavy_metal[922350] \
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
				print('  fuel radius    : ', self.inputs.fuel_cell_radius)
				print('  void radius    : ', self.inputs.void_cell_radius)
				print('  clad radius    : ', self.inputs.clad_cell_radius)
				print('  fuel density   : ', self.inputs.fuel_density)
				print('  clad density   : ', self.inputs.clad_density)
				print('  coolant density: ', self.inputs.cool_density)
			print(' Path:', self.ip_path)
			print(' Normalized values')
			print(self.normalized)
			if(detail):
				print('  clad density: ', self.norm_clad_density)
			print('  cool density: ', self.norm_cool_density)
			print('  enrichment  : ', self.norm_enrichment)		
			print(' output information:')
			print('  max prod: ', self.max_prod)
			print('  max dest: ', self.max_dest)
			print('  max BU  : ', self.max_BU)

	def Coordinates(self, varied_ips):
		coordinates = {}
		
		for key, value in self.normalized.items():
			if key in varied_ips:
				coordinates[key] = value
		
		return coordinates
	
class DBase:
	""" A class that handles all generated libraries """
	slibs = []		# Scout libs
	flibs = []		# Full libs
	
	complete_slibs = 0
	complete_flibs = 0
	
	# Library neighborhoods
	slib_neighbors = []
	flib_neighbors = []
	
	# Scout lib output parameters
	max_prods = []
	max_dests = []
	max_BUs = []
	
	# Base case parameters
	fuel_cell_radius = []
	clad_cell_radius = []
	clad_cell_thickness = []
	void_cell_radius = []
	
	enrichment = []
	
	# Min and max values ("range") in database
	range_fuel_radius = [1,0]	# [min,max]
	
	range_fuel_density = [1,0]
	range_clad_density = [1,0]
	range_cool_density = [1,0]
	
	range_enrichment = [1,0]
	
	
	def __init__(self, paths, dimensions):		
		# Read database, assuming software may have been interrupted and
		#	the folder may have some inputs and outputs 
		
		if not os.path.isdir(paths.database_path):
			raise RuntimeError('The database path does not exist')
			
		if not os.path.exists(paths.database_path + paths.base_input):
			error_message = 'The database base-case input file does not exist. Looking for: ' \
							+ paths.database_path + paths.base_input
			raise RuntimeError(error_message)
		
		self.name = paths.database_name
		self.dimensions = dimensions
		
		# Read database inputs
		self.ReadInput(paths.database_path + paths.dbase_input)
		
		if self.dimensions != self.ip_ranges.DefinedCount():
			raise RuntimeError('Dimensions and database input file varied inputs number mismatch')
		print('Database dimensions:', self.dimensions)
		
		# Check to see if there's an scout library input folder
		if os.path.exists(paths.database_path + paths.SR_Input_folder):
			tot_files = len(os.listdir(paths.database_path + paths.SR_Input_folder))
		# If the input scout library folder doesn't exist create it
		else:
			os.mkdir(paths.database_path + paths.SR_Input_folder)
			tot_files = 0
		
		
		tot_sr_libs = 0
		# If there are files SR_Inputs, read them
		if tot_files > 0:
			print('Attempting to read', tot_files, 'scout libraries')
			for ip_number in range(tot_files):
					ip_path = paths.database_path + paths.SR_Input_folder + '/' + str(ip_number) +'.py'
					op_path = paths.database_path + paths.SR_Output_folder + '/' + paths.xsgen_prefix \
								+ paths.sr_prefix + str(ip_number) + '/' + paths.xsgen_op_folder
								
					if os.path.exists(ip_path):
						inputlib = Library(database_path=paths.database_path, op_path=op_path, ip_path=ip_path, number=ip_number, scout=True)
						self.AddLib(inputlib)
						tot_sr_libs += 1
					else:
						# Could continue here instead of breaking, but at this point this is better
						break
							
		if tot_sr_libs < tot_files:
			print(tot_sr_libs, 'inputs read, however', tot_files - tot_sr_libs, \
						'files in', paths.database_path + paths.SR_Input_folder, 'were ignored')
		
		
	# Database specific immutable parameters
	#TODO: combining fractions
	
	def ReadInput(self, ip_path):
		# Make sure it's there
		if not os.path.exists(ip_path):
			error_message = 'The database input file does not exist. Looking for: ' \
							+ ip_path
			raise RuntimeError(error_message)
		
		self.ip_ranges = xsgenParams()	# Ranges of varied inputs
		self.varied_ips = []	# the varied inputs
		
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
			except (ValueError):
				pass
			
			# Check the rest of defined inputs
			if len(items) == 2 and float(items[1]) > 0:
				if items[0] in self.ip_ranges.geometry:
					self.ip_ranges.geometry[items[0]] = float(items[1])
					self.varied_ips.append(items[0])
				if items[0] in self.ip_ranges.density:
					self.ip_ranges.density[items[0]] = float(items[1])
					self.varied_ips.append(items[0])
				if items[0] in self.ip_ranges.other:
					self.ip_ranges.other[items[0]] = float(items[1])
					self.varied_ips.append(items[0])
		print('Varied inputs:', self.varied_ips)
		
	
	# stored libraries
	#TODO: should check if exists first
	#TODO: should update relevant neighborhoods
	def AddLib(self, added_lib):
		if added_lib.scout:
			if added_lib.number != len(self.slibs):
				raise RuntimeError('Scout library number mismatch when adding to slibs')
			self.slibs.append(copy.deepcopy(added_lib))
		else:
			if added_lib.number != len(self.flibs):
				raise RuntimeError('Full library number mismatch when adding to slibs')
			self.flibs.append(copy.deepcopy(added_lib))
	
	#TODO: this probably isnt used and aint even right
	def Exists(self, number):
		for lib in self.slibs:
			if lib.number == number:
				return True
		return False
	
	#TODO: divide to completed and queued 
	#TODO: currently only recreates, make it update
	# information about current database
	def UpdateData(self):
		# Delete old data
		self.max_prods.clear()
		self.max_dests.clear()
		self.max_BUs.clear()
		
		self.fuel_cell_radius.clear()
		self.clad_cell_radius.clear()
		self.clad_cell_thickness.clear()
		self.void_cell_radius.clear()
		self.enrichment.clear()
			
		# Rebuild lib values
		for i in self.slibs:
			if i.completed:
				complete_slibs += 1
				self.fuel_cell_radius.append(i.inputs.fuel_cell_radius)
				self.clad_cell_radius.append(i.inputs.clad_cell_radius)
				self.clad_cell_thickness.append(i.inputs.clad_cell_radius - i.inputs.fuel_cell_radius)
				self.void_cell_radius.append(i.inputs.void_cell_radius)
				self.enrichment.append(i.inputs.enrichment)
				
				self.max_prods.append(i.max_prod)
				self.max_dests.append(i.max_dest)
				self.max_BUs.append(i.max_BU)
			
	def PCA(self):
		if len(self.max_prods) > 0:
			self.np_prods = np.asarray(self.max_prods)
			self.np_dests = np.asarray(self.max_dests)
			self.np_BUs = np.asarray(self.max_BUs)
			
			self.data_mat = np.column_stack((self.np_prods, self.np_dests, self.np_BUs))
			self.pca_mat = mlabPCA(self.data_mat)	# PCA matrix 
		
	def UpdateMetrics(self):
		# work in progress
		#TODO: have a class to handle metrics
		# Update the range of metrics
		self.range_fuel_radius[0] = min([i.inputs.geometry['fuel_cell_radius'] for i in self.slibs])
		self.range_fuel_radius[1] = max([i.inputs.geometry['fuel_cell_radius'] for i in self.slibs])
		
		self.range_fuel_density[0] = min([i.inputs.density['fuel_density'] for i in self.slibs])
		self.range_fuel_density[1] = max([i.inputs.density['fuel_density'] for i in self.slibs])
		
		self.range_clad_density[0] = min([i.inputs.density['clad_density'] for i in self.slibs])
		self.range_clad_density[1] = max([i.inputs.density['clad_density'] for i in self.slibs])
		
		self.range_cool_density[0] = min([i.inputs.density['cool_density'] for i in self.slibs])
		self.range_cool_density[1] = max([i.inputs.density['cool_density'] for i in self.slibs])
		
		self.range_enrichment[0] = min([i.inputs.other['enrichment'] for i in self.slibs])
		self.range_enrichment[1] = max([i.inputs.other['enrichment'] for i in self.slibs])
		
		# Update the metrics in libraries
		for lib in self.slibs:
			lib.normalized['fuel_cell_radius'] = (lib.inputs.geometry['fuel_cell_radius'] - self.range_fuel_radius[0]) / \
									(self.range_fuel_radius[1] - self.range_fuel_radius[0])
									
			lib.normalized['fuel_density'] = (lib.inputs.density['fuel_density'] - self.range_fuel_density[0]) / \
									(self.range_fuel_density[1] - self.range_fuel_density[0])
									
			lib.normalized['clad_density'] = (lib.inputs.density['clad_density'] - self.range_clad_density[0]) / \
									(self.range_clad_density[1] - self.range_clad_density[0])
									
			lib.normalized['cool_density'] = (lib.inputs.density['cool_density'] - self.range_cool_density[0]) / \
									(self.range_cool_density[1] - self.range_cool_density[0])
									
			lib.normalized['enrichment'] = (lib.inputs.other['enrichment'] - self.range_enrichment[0]) / \
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
			for subset in itertools.combinations(list1, self.dimensions*2):
				n_coordinates = []
				for i in subset:
					n_coordinates.append(self.slibs[i].Coordinates(self.varied_ips))
				new_neighborhood = Neighborhood(current_coords, subset, n_coordinates)
				# Update the current data if new has better score
				if current_score < new_neighborhood.neighbor_score:
					current_score = new_neighborhood.neighbor_score
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
		print("Database scout run ranges: ")
		print("  Fuel radius range:  ", self.range_fuel_radius)
		print("  Fuel density range: ", self.range_fuel_density)
		print("  Clad denisty range: ", self.range_clad_density)
		print("  Cool density range: ", self.range_cool_density)
		print("  Enrichment range:   ", self.range_enrichment)
		#print("  Max prods: ", self.max_prods)
		#print("  Max desds: ", self.max_dests)
		#print("  Max BUs  : ", self.max_BUs)
		
		





















