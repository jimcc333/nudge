import os

import numpy as np
from matplotlib.mlab import PCA as mlabPCA

class Library:
	""" A class that holds library information """ 
	#TODO "run" routine to talk to xsgen, needs to be parallizable
	
	def __init__(self, op_path, ip_path, number):
		self.ip_path = ip_path		# Path of the input file
		self.op_path = op_path	# path to the library folder w/ brightlite0 in it
		self.number = number	# Unique number of the library
		
		self.max_prod = 0
		self.max_dest = 0
		self.max_BU = 0
		
		# Read input
		self.ReadInput(ip_path)
		
		#TODO: pass combining fractions (frac) better
		
		if os.path.isdir(op_path + "build-sr" + str(number) + "/"):
			self.scout = True
			self.completed = True
			
			u235_file = op_path + "build-sr" + str(number) + "/brightlite0/922350.txt"
			self.ReadOutput("U235", u235_file, 0.04)
			
			u238_file = op_path + "build-sr" + str(number) + "/brightlite0/922380.txt"
			self.ReadOutput("U235", u238_file, 0.96)
			
			print("Completed reading scout output #" + str(number))
		elif os.path.isdir(op_path + "build-fr" + str(number) + "/"):
			self.scout = False
			self.completed = True
			
			#TODO: read full run output
			
		else:	# library not run yet
			self.completed = False
			
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
				self.fuel_cell_radius = float(items[2])
			if items[0] == "void_cell_radius":
				self.void_cell_radius = float(items[2])
			if items[0] == "clad_cell_radius":
				self.clad_cell_radius = float(items[2])
			if items[0] == "unit_cell_pitch":
				self.unit_cell_pitch = float(items[2])
			if items[0] == "unit_cell_height":
				self.unit_cell_height = float(items[2])
			if items[0] == "fuel_density":
				self.fuel_density = float(items[2])
			if items[0] == "clad_density":
				self.clad_density = float(items[2])
			if items[0] == "cool_density":
				self.cool_density = float(items[2])
			if items[0] == "enrichment":
				self.enrichment = float(items[2])
			if items[0] == "flux":
				self.flux = float(items[2])
			if items[0] == "k_particles":
				self.k_particles = float(items[2])
			if items[0] == "k_cycles":
				self.k_cycles = float(items[2])
			if items[0] == "k_cycles_skip":
				self.k_cycles_skip = float(items[2])
		
	def Print(self):
		print('Lib #', self.number, ' input information:')
		print('  fuel radius: ', self.fuel_cell_radius)
		print('  void radius: ', self.void_cell_radius)
		print('  clad radius: ', self.clad_cell_radius)
		print('  fuel density: ', self.fuel_density)
		print('  clad density: ', self.clad_density)
		print('  coolant density: ', self.cool_density)
		print(' output information:')
		print('  max prod: ', self.max_prod)
		print('  max dest: ', self.max_dest)
		print('  max BU: ', self.max_BU)
	
	# --- Inputs ---
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
	
	# --- Outputs ---
	
class DBase:
	""" A class that handles all generated libraries """
	slibs = []		# Completed scout libs
	max_prods = []
	max_dests = []
	max_BUs = []
	fuel_cell_radius = []
	clad_cell_radius = []
	clad_cell_thickness = []
	void_cell_radius = []
	enrichment = []
	
	def __init__(self, name):
		self.name = name
		
	# Database specific immutable parameters
	#TODO: combining fractions
	
	# stored libraries
	def AddSLib(self, added_lib):
		self.slibs.append(added_lib)
	
	def Exists(self, number):
		for lib in self.slibs:
			if lib.number == number:
				return True
		return False
	
	#TODO: make this actually update and then divide to completed and queued 
	# information about current database
	def UpdateData(self):
		self.max_prods.clear()
		self.max_dests.clear()
		self.max_BUs.clear()
		
		self.fuel_cell_radius.clear()
		self.clad_cell_radius.clear()
		self.clad_cell_thickness.clear()
		self.void_cell_radius.clear()
		self.enrichment.clear()
		
		for i in self.slibs:
			self.fuel_cell_radius.append(i.fuel_cell_radius)
			self.clad_cell_radius.append(i.clad_cell_radius)
			self.clad_cell_thickness.append(i.clad_cell_radius - i.fuel_cell_radius)
			self.void_cell_radius.append(i.void_cell_radius)
			self.enrichment.append(i.enrichment)
			
			self.max_prods.append(i.max_prod)
			self.max_dests.append(i.max_dest)
			self.max_BUs.append(i.max_BU)
			
		
	def PCA(self):
		self.np_prods = np.asarray(self.max_prods)
		self.np_dests = np.asarray(self.max_dests)
		self.np_BUs = np.asarray(self.max_BUs)
		
		self.data_mat = np.column_stack((self.np_prods, self.np_dests, self.np_BUs))
		self.pca_mat = mlabPCA(self.data_mat)	# PCA matrix 
		
	def Print(self):
		print("Database scout run max values: ")
		print("  Max prods: ", self.max_prods)
		print("  Max desds: ", self.max_dests)
		print("  Max BUs  : ", self.max_BUs)
		
		
	



















