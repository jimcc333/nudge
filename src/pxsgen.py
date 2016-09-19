#!/usr/bin/python3
#
#
#  Copyright 2016 Cem Bagdatli
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#
#  Inputs:
#      1. enrichment
#
#
#
#  Notes:
#    - Does not use fuel composition: just enrichment
#    - All inputs are [0,1]
#
#
#
#

from objects import xsgenParams, Library
import os
import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main(args):
    os.system('cls' if os.name == 'nt' else 'clear')
    print('Placeholder xsgen program')

    # Inputs and constants
    #lib = Library('blank', 'blank', args[1], 1, False)

    inputs = {
        'fuel_radius': None,        # 1
        'void_thickness': None,     # 2
        'clad_thickness': None,     # 3
        'unit_cell_pitch': None,    # 4
        'unit_cell_height': None,   # 5
        'fuel_density': None,       # 6
        'clad_density': None,       # 7
        'cool_density': None,       # 8
        'enrichment': None,            # 9
        'flux': None               # 10
    }

    outputs = {
        'burnup': None,         # 9, 1, 6, 3, 4, 5
        'neutron_prod': None,   # 10, 7, 2, 8, 1, 9
        'neutron_dest': None    # 10, 8, 3, 7, 4, 9
    }

    #print('inputs: ', inputs)

    # assign(inputs, 0.5)
    randomize(inputs)

    print(inputs)

    # Burnup
    outputs['burnup'] = burnup_maker(inputs)



    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    x = []
    y = []
    z = []

    for xx in np.arange(0,1,0.05):
        for yy in np.arange(0,1,0.05):
            inputs['enrichment'] = xx
            inputs['unit_cell_pitch'] = yy
            x.append(xx)
            y.append(yy)
            z.append(burnup_maker(inputs))

    ax1.scatter(x,y,z) # s=sizes, c=colors)

    ax1.set_xlim([0,1])
    ax1.set_ylim([0,1])
    #ax1.set_zlim([0,1])
    plt.show()


    print('\n-TheEnd-')
    return 0

def randomize(ip_data):
    for key, value in ip_data.items():
        if value == None:
            ip_data[key] = round(random.random(), 4)

def assign(ip_data, x):
    for key, value in ip_data.items():
        if value == None:
            ip_data[key] = x

def burnup_maker(inputs):
    x1 = (2 + inputs['enrichment'])**3
    x2 = (2.5635 - inputs['cool_density'])**2.12
    x3 = (3 - inputs['clad_density'])**2.1234
    x4 = (0.5+inputs['fuel_density'])**3.01
    x5 = np.sin(inputs['fuel_radius']*2)*8 + 2

    return ((x1 + x2 + x3 + x4 + x5)**2.3)/100 + 3*inputs['void_thickness'] + 2*inputs['clad_thickness']


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))



























