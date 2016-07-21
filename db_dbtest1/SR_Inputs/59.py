import numpy as np
from numpy import logspace
from xsgen.nuc_track import transmute

#############################
### General specifcations ###
#############################
reactor = "sr59" 
plugins = ['xsgen.pre', 'xsgen.buk']
solver = 'openmc+origen'
formats = ('brightlite',)
burn_regions = 1     # Number of burnup annular regions.
burn_time = 365*10  # Number of days to burn the material [days]
time_step = 100      # Time step by which to increment the burn [days]
burn_times = [0, 3] #+ list(range(10, 101, 10))
burn_times.extend(range(100, 301, 100))  # we now have [0, 3, 100, 200 .. 4000]

batches = 3


################################
### Unit Cell Sepcifications ###
################################
# all in [cm]
fuel_cell_radius = 0.376552
void_cell_radius = 0.385052
clad_cell_radius = 0.443328
unit_cell_pitch  = 0.65635 * 2.0
unit_cell_height = 10.0

fuel_density = 9.43103
clad_density = 5.05172
cool_density = 0.748276

# choose one of the following
flux = 3e14  # the average reactor flux in [n/cm2/s]
#fuel_specific_power = 40.0   # Power garnered from fuel [W / g]

# LEU
initial_heavy_metal = {     # Initial heavy metal mass fraction distribution
    922350: 0.0262069,
    922380: 0.973793,
    }

enrichment = 0.0262069

pnl = 0.96

# UOX
fuel_chemical_form = {                 # Dictionary of initial fuel loading.
    80160: 2.0,
    "IHM": 1.0,
    }

k_particles   = 500       # Number of particles to run per kcode cycle
k_cycles      = 100       # Number of kcode cycles to run
k_cycles_skip = 35        # Number of kcode cycles to run but not tally at the begining.

group_structure = [1.0e-9, 10]
openmc_group_struct = np.logspace(1, -9, 101)

# Temperature
# Should be a positive multiple of 300 K (ie 300, 600, 900, etc)
temperature = 293.6

track_nucs = ["Ac227",
              "Am241",
              "AM242",
              "BA140",
              "C14",
              "CM251",
              "CS141",
              "CS142",
              "CS147",
              "H1",
              "H3",
              "PU236",
              "PU237",
              "Pu238",
              "Pu239",
              "Pu240",
              "Pu241",
              "Th228",
              "Th229",
              "Th230",
              "Th232",
              "U230",
              "U231",
              "U232",
              "U233",
              "U234",
              "U235",
              "U236",
              "U237",
              "U238",
              "U239",
              "Zr93",
              "ZR95",
]