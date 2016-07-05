class Library:
	""" A class that holds library information """ 
	
	def __init__(self, path, number):
		self.path = path		# Path of the library folder
		self.number = number	# Unique number of the library
	
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
	# Database specific immutable parameters
	
	# stored libraries
	
	# stored information about current database
	

