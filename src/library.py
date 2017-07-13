import os
from objects import xsgenParams


class Library:
    """ A class that holds library information """

    def __init__(self, database_path, op_path, ip_path, number, scout):
        # Library constant and inputs
        self.database_path = database_path  # Path to the database
        self.ip_path = ip_path  # Path of the input file
        self.op_path = op_path  # path to the library folder w/ Bright-lite formatted .txt files in it
        self.number = number    # Unique number of the library
        self.scout = scout      # True/False whether it is a scout library
        self.inputs = xsgenParams()

        # Library output information
        self.max_prod = 0       # Maximum neutron production value in library
        self.max_dest = 0       # Maximum neutron destruction value in library
        self.max_BU = 0         # Maximum burnup in library

        # Library analysis values from database
        self.coordinate = []            # The coordinates (based on normalized varied inputs) of the library
        self.voronoi_size = 0           # The Voronoi cell size of the library
        self.excluded_error = 0         # The error of the database if this point is excluded
        self.furthest_point = []        # The furthest found point within the Voronoi cell
        self.furthest_point_dist = 0    # The distance of the furthest point
        self.rank = 0                   # The rank of the library, used during exploitation
        self.proximity_order = []       # Indexes of libraries in database ordered from closest to furthest

        # --- Normalized Input Values ---
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
        self.read_input(ip_path)

        # Read output if exists
        # TODO: pass combining fractions (frac) better
        if os.path.isdir(op_path):
            self.completed = True

            u235_file = op_path + "/922350.txt"
            self.read_output(u235_file, 0.04)

            u238_file = op_path + "/922380.txt"
            self.read_output(u238_file, 0.96)
        elif os.path.exists(op_path):
            self.completed = True
            self.read_output(op_path, 1)
        else:
            self.completed = False
            # TODO: read full run output
            # TODO: add library to queue

    # Returns the normalized [0,1] array of varied inputs (order REALLY matters here)
    def coordinates(self, varied_ips):
        coordinates = []
        for key, value in sorted(self.normalized.items()):
            if key in varied_ips:
                coordinates.append(value)
        self.coordinate = coordinates
        return coordinates

    def print(self, detail=0):
        if self.ip_path == 'x':
            print('Interpolated library output information: ')
            print('  Max prod: ', self.max_prod)
            print('  Max dest: ', self.max_dest)
            print('  Max BU  : ', self.max_BU)
        else:
            print('Lib #', self.number, ' input information:')
            if detail:
                print(self.inputs.xsgen)
            print(' Path:', self.ip_path)
            print(' Normalized values')
            print(self.normalized)

    def read_input(self, ip_path):
        if ip_path == "x":
            return

        try:
            doc = open(ip_path, "r")
        except IOError:
            print("Could not open ", ip_path)
            return

        max_lines = 500
        ip = doc.readlines()

        for line_i in range(len(ip)):
            items = ip[line_i].split()

            if len(items) < 3:
                continue
            if items[0] in self.inputs.xsgen:
                try:
                    self.inputs.xsgen[items[0]] = float(items[2])
                except ValueError:
                    self.inputs.xsgen[items[0]] = str(items[2])
            if items[0] == 'initial_heavy_metal':
                while line_i < max_lines:
                    line_i += 1
                    items = ip[line_i].replace(':',' ').replace(',',' ').split()
                    if len(items) > 2:
                        error_message = 'Input file in ' + self.ip_path + ' has formatting error at initial_heavy_' \
                                                                          'metal. Make sure each NUCID is in a new ' \
                                                                          'line and close bracket (}) at a new line.'
                        raise RuntimeError(error_message)
                    if '}' not in items[0]:
                        continue    # self.inputs.initial_heavy_metal[int(items[0])] = float(items[1])
                    else:
                        # End of initial heavy metal
                        break
        '''
        if self.inputs.xsgen['enrichment'] is not None:
            if self.inputs.xsgen['enrichment'] - (self.inputs.initial_heavy_metal[922350] /
                                                  (self.inputs.initial_heavy_metal[922350] +
                                                   self.inputs.initial_heavy_metal[922380])) > 0.001:
                    error_message = 'Input file in ' + self.ip_path + \
                                    ' has inconsistency between enrichment and given mass compositions'
                    raise RuntimeError(error_message)
        '''

    # Reads the xsgen output of library
    def read_output(self, file_path, frac):
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
                self.max_BU += sum([float(i) for i in items[2:]]) * frac
        doc.close()

    # Resets variables that are not used for calculations and outputs
    def reset(self):
        self.max_BU = 0
        self.max_prod = 0
        self.max_dest = 0
        self.ip_path = 'x'
        self.op_path = 'x'
        self.number = 'x'
