import os
import copy

import numpy as np
from matplotlib.mlab import PCA as mlabPCA
from operator import attrgetter

class PathNaming:
	def __init__(self, database_path = "/home/cem/nudge/db_dbtest1/"):
		self.database_path = database_path
		
		self.base_input = 'basecase.py'		# This file is used as the base case for the database
		
		self.SR_Input_folder = 'SR_Inputs'	# Where scout run inputs are stored
		self.SR_Output_folder = 'SR_Outputs'	# Where scout run outputs are stored
		self.FR_Input_folder = 'FR_Inputs'	# Where full run inputs are stored
		self.FR_Output_folder = 'FR_Outputs'	# Where full run outputs are stored
		
		self.xsgen_prefix = 'build-'			# The prefix xsgen assigns to output folders
		self.xsgen_op_folder = 'brightlite0'	# The folder xsgen places bright-lite formatted outputs
		
		self.sr_prefix = 'sr'				# The prefix given to scout libs
		self.fr_prefix = 'fr'				# The prefix given to full libs
		
		self.database_name = 'database1'
		
class xsgenParams:
	def __init__(self):
		# Initial heavy metal mass fraction distribution
		initial_heavy_metal = {
			922350: 0.033,
			922380: 0.967,
		}
		
		# Geometry inputs
		fuel_cell_radius = 0.410			# [cm]
		void_cell_radius = 0.4185			# [cm]
		clad_cell_radius = 0.475			# [cm]
		unit_cell_pitch  = 0.65635 * 2.0	# [cm]
		unit_cell_height = 10.0				# [cm]
		
		# Density inputs
		fuel_density = 10.7  # Fuel density [g/cc]
		clad_density = 5.87  # Cladding Density [g/cc]
		cool_density = 0.73  # Coolant Density [g/cc]

		# Others
		flux = 3e14  			# Average reactor flux [n/cm2/s]
		k_particles   = 5000	# Number of particles to run per kcode cycle
		k_cycles      = 130		# Number of kcode cycles to run
		k_cycles_skip = 30		# Number of kcode cycles to run but skip
		group_structure = [1.0e-9, 10]
		#openmc_group_struct = np.logspace(1, -9, 101)


class Library:
	""" A class that holds library information """ 
	#TODO "run" routine to talk to xsgen, needs to be parallizable
	
	def __init__(self, database_path, op_path, ip_path, number, scout):
		
		self.ip_path = ip_path	# Path of the input file
		self.op_path = op_path	# path to the library folder w/ Bright-lite formatted .txt files in it
		self.number = number	# Unique number of the library
		self.scout = scout
		self.inputs = PathNaming(database_path)
		
		self.max_prod = 0
		self.max_dest = 0
		self.max_BU = 0
		
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
		
		for line in doc.readlines():
			items = line.split()
			
			if len(items) < 3:
				continue
			if items[0] == "fuel_cell_radius":
				self.inputs.fuel_cell_radius = float(items[2])
			if items[0] == "void_cell_radius":
				self.inputs.void_cell_radius = float(items[2])
			if items[0] == "clad_cell_radius":
				self.inputs.clad_cell_radius = float(items[2])
			if items[0] == "unit_cell_pitch":
				self.inputs.unit_cell_pitch = float(items[2])
			if items[0] == "unit_cell_height":
				self.inputs.unit_cell_height = float(items[2])
			if items[0] == "fuel_density":
				self.inputs.fuel_density = float(items[2])
			if items[0] == "clad_density":
				self.inputs.clad_density = float(items[2])
			if items[0] == "cool_density":
				self.inputs.cool_density = float(items[2])
			if items[0] == "enrichment":
				self.inputs.enrichment = float(items[2])
			if items[0] == "flux":
				self.inputs.flux = float(items[2])
			if items[0] == "k_particles":
				self.inputs.k_particles = float(items[2])
			if items[0] == "k_cycles":
				self.inputs.k_cycles = float(items[2])
			if items[0] == "k_cycles_skip":
				self.inputs.k_cycles_skip = float(items[2])
		
	def Print(self, detail=0):
		if self.ip_path == "x":
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
			print('  fuel radius : ', self.norm_fuel_radius)
			print('  fuel density: ', self.norm_fuel_density)
			if(detail):
				print('  clad density: ', self.norm_clad_density)
			print('  cool density: ', self.norm_cool_density)
			print('  enrichment  : ', self.norm_enrichment)		
			print(' output information:')
			print('  max prod: ', self.max_prod)
			print('  max dest: ', self.max_dest)
			print('  max BU  : ', self.max_BU)
	
	# --- Inputs ---
	
	# --- Normalized Values ---
	norm_fuel_radius = 0.5
	
	norm_fuel_density = 0.5
	norm_clad_density = 0.5
	norm_cool_density = 0.5
	
	norm_enrichment = 0.5
	
	
	# --- Outputs ---
	
class DBase:
	""" A class that handles all generated libraries """
	slibs = []		# Scout libs
	flibs = []		# Full libs
	
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
	
	
	def __init__(self, paths):
		
		# Read database, assuming software may have been interrupted and
		#	the folder may have some inputs and outputs 
		
		if not os.path.isdir(paths.database_path):
			raise RuntimeError('The database path does not exist')
			
		if not os.path.exists(paths.database_path + paths.base_input):
			error_message = 'The database base-case input file does not exist. Looking for: ' \
							+ paths.database_path + paths.base_input
			raise RuntimeError(error_message)
		
		self.name = paths.database_name
		
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
	
	# stored libraries
	#TODO: should check if exists first
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
		self.range_fuel_radius[0] = min([i.inputs.fuel_cell_radius for i in self.slibs])
		self.range_fuel_radius[1] = max([i.inputs.fuel_cell_radius for i in self.slibs])
		
		self.range_fuel_density[0] = min([i.inputs.fuel_density for i in self.slibs])
		self.range_fuel_density[1] = max([i.inputs.fuel_density for i in self.slibs])
		
		self.range_clad_density[0] = min([i.inputs.clad_density for i in self.slibs])
		self.range_clad_density[1] = max([i.inputs.clad_density for i in self.slibs])
		
		self.range_cool_density[0] = min([i.inputs.cool_density for i in self.slibs])
		self.range_cool_density[1] = max([i.inputs.cool_density for i in self.slibs])
		
		self.range_enrichment[0] = min([i.inputs.enrichment for i in self.slibs])
		self.range_enrichment[1] = max([i.inputs.enrichment for i in self.slibs])
		
		# Update the metrics in libraries
		for lib in self.slibs:
			lib.norm_fuel_radius = (lib.inputs.fuel_cell_radius - self.range_fuel_radius[0]) / \
									(self.range_fuel_radius[1] - self.range_fuel_radius[0])
									
			lib.norm_fuel_density = (lib.inputs.fuel_density - self.range_fuel_density[0]) / \
									(self.range_fuel_density[1] - self.range_fuel_density[0])
									
			lib.norm_clad_density = (lib.inputs.clad_density - self.range_clad_density[0]) / \
									(self.range_clad_density[1] - self.range_clad_density[0])
									
			lib.norm_cool_density = (lib.inputs.cool_density - self.range_cool_density[0]) / \
									(self.range_cool_density[1] - self.range_cool_density[0])
									
			lib.norm_enrichment = (lib.inputs.enrichment - self.range_enrichment[0]) / \
									(self.range_enrichment[1] - self.range_enrichment[0])
		
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
		
		





















