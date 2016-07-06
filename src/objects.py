import numpy as np

class Library:
	""" A class that holds library information """ 
	
	def __init__(self, path, number):
		self.path = path		# Path of the library folder w/ brightlite0
		self.number = number	# Unique number of the library
		
		self.max_prod = 0
		self.max_dest = 0
		self.max_BU = 0
		
		#TODO: pass combining fractions (frac) better
		
		u235_file = path + "/brightlite0/922350.txt"
		self.Read("U235", u235_file, 0.04)
		
		u238_file = path + "/brightlite0/922380.txt"
		self.Read("U235", u238_file, 0.96)
	
		
	def Read(self, nuclide, file_path, frac):	
		try:
			doc = open(file_path, "r")
		except IOError:
			print "Could not open ", file_path
					
		for line in doc.readlines():
			items = line.split()
			
			if items[0] == "NEUT_PROD":
				self.max_prod += float(items[len(items)-1]) * frac
			if items[0] == "NEUT_DEST":
				self.max_dest += float(items[len(items)-1]) * frac
			if items[0] == "BUd":
				self.max_BU += sum( [float(i) for i in items[1:]] ) * frac
		doc.close()
	
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
	libs = []
	max_prods = []
	max_dests = []
	max_BUs = []
	
	def __init__(self, name):
		self.name = name
		
	# Database specific immutable parameters
	
	
	# stored libraries
	def Add(self, added_lib):
		self.libs.append(added_lib)
	
	# information about current database
	def UpdateData(self):
		for i in self.libs:
			self.max_prods.append(i.max_prod)
			self.max_dests.append(i.max_dest)
			self.max_BUs.append(i.max_BU)
		
	def PCA(self):
		self.np_prods = np.asarray(self.max_prods)
		self.np_dests = np.asarray(self.max_dests)
		self.np_BUs = np.asarray(self.max_BUs)
		
		self.data_mat = np.column_stack((self.np_prods, self.np_dests, self.np_BUs))
		
		
		
	



















