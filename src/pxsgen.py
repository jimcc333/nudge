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
    print('---Placeholder XSgen program---')

    if len(args) != 3:
        print('Please use "pxsgen input_location output_destination" format')
        return

    # Read inputs
    lib = Library('blank', args[2], args[1], 0, False)

    # Write outputs
    outputs = ''
    outputs = 'BUd = ' + str(burnup_maker(lib.inputs.xsgen)) + '\n'\
              + 'NEUT_PROD = ' + str(prod_maker(lib.inputs.xsgen)) + '\n'\
              + 'NEUT_DEST = ' + str(dest_maker(lib.inputs.xsgen)) + '\n'

    with open(args[2], 'w') as openfile:
        openfile.write(outputs)

    print('Generated output')

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
    random.seed(inputs['enrichment']*inputs['flux']*1245)

    x1 = (2 + inputs['enrichment'])**3
    x2 = (2.5635 - inputs['cool_density'])**2.12
    x3 = (3 - inputs['clad_density'])**2.1234
    x4 = (0.5+inputs['fuel_density'])**3.01
    x5 = np.sin(inputs['fuel_cell_radius']*2)*8 + 2
    x6 = 8 - ((inputs['flux'] - 0.71)*2)**4
    x7 = 3*inputs['void_cell_radius'] + 2*inputs['clad_cell_radius'] + random.random()

    random.seed()

    return ((x1 + x2 + x3 + x4 + x5 + x6)**2.3)/100 + x7


def prod_maker(inputs):
    random.seed(inputs['fuel_cell_radius']*inputs['clad_density']*1245)
    x1 = (1.5-inputs['fuel_cell_radius'])**1.23
    x2 = 9.45**inputs['void_cell_radius']
    x3 = (np.sin(inputs['unit_cell_pitch']*12) + 1)/10
    x4 = (2 - inputs['clad_density'])**2.1234
    x5 = (inputs['enrichment'] + 2)**3
    x6 = (inputs['flux'] + 1)**3.245
    x7 = 12 - ((inputs['enrichment'] - 0.3)*2)**4 + random.random()*1.2

    random.seed()

    return (x2 + x4 + x5 + x6)**0.5 + x1 + x3 + x7


def dest_maker(inputs):
    random.seed(inputs['clad_cell_radius']*inputs['cool_density']*1245)
    x1 = 9.45**inputs['void_cell_radius']
    x2 = (inputs['unit_cell_height'] + 2.12)**4.1 / 100
    x4 = inputs['cool_density']*3.14159
    x5 = (inputs['enrichment'] + 2)**3
    x6 = (inputs['flux'] + 2)**2.254
    x7 = 11 - ((inputs['enrichment'] - 0.4)*2)**4 + random.random()*1.1

    random.seed()

    return (x6 * (x5 - x4 + x1))**0.5 - x2 + x7


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))



























