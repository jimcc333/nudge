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

from library import Library
import random
import numpy as np


def main(args, inputs=None):
    max_dimensions = 9

    # Check if inputs are passed to main
    write_out = False
    if inputs is None:
        # Read inputs from file
        write_out = True
        lib = Library('blank', args[2], args[1], 0, False)
        inputs = lib.inputs.xsgen

    # print('---Placeholder XSgen program---')

    # Write outputs
    all_inputs = [0.5 for i in range(max_dimensions)]
    all_inputs[0] = inputs['fuel_density']
    all_inputs[1] = inputs['clad_density']
    all_inputs[2] = inputs['cool_density']
    all_inputs[3] = inputs['enrichment']
    all_inputs[4] = inputs['flux']
    all_inputs[5] = inputs['fuel_cell_radius']
    all_inputs[6] = inputs['clad_cell_radius']
    all_inputs[7] = inputs['void_cell_radius']
    all_inputs[8] = inputs['unit_cell_pitch']
    output = output_method(all_inputs)

    outputs = 'BUd = ' + str(output) + '\n' + 'NEUT_PROD = ' + str(0) + '\n' + 'NEUT_DEST = ' + str(0) + '\n'

    if write_out:
        with open(args[2], 'w') as openfile:
            openfile.write(outputs)

    return output


def output_method(inputs):
    output = f6(inputs[0], inputs[1], inputs[2]) * 2
    # output += f6(inputs[3], inputs[4], inputs[5])
    # output += f5(inputs[6], inputs[7], inputs[8])
    return output


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


""""
The functions below have been taken from:
Testing Methods for 3D Scattered Data Interpolation
by Mira Bozzini and Milvia Rossini. Monograf´ıas de la Academia de Ciencias de Zaragoza. 20: 111–135, (2002).

"""

def f2(x, y, z):
    return (np.tanh(9 * x - 9 * z - 9 * y) + 1) / 9


# Good
def f3(x, y, z):
    return np.cos(6 * z) * (1.25 + np.cos(5.4 * y)) / (6 + 6 * (3 * x - 1) ** 2)


# Hump at center
def f4(x, y, z):
    return np.exp(-81 / 16 * ((x - 0.5) ** 2 + (y - 0.52) ** 2 + (z - 0.47) ** 2))

# Another hump with steeper drop
def f5(x, y, z):
    return np.sqrt(64 - 81 * ((x - 0.5) ** 2 + (y - 0.52) ** 2 + (z - 0.47) ** 2)) / 9


# perfect at z=0.5
def f6(x, y, z):
    x1 = 0.75 * np.exp(-0.25 * ((9 * x - 2) ** 2 + (9 * y - 2) ** 2 + (9 * z - 2) ** 2))
    x2 = 0.75 * np.exp(- (9 * x + 1) ** 2 / 49 - (9 * y + 1) ** 2 / 10 - (9 * z + 1) ** 2 / 10)
    x3 = 0.50 * np.exp(-0.25 * ((9 * x - 7) ** 2 + (9 * y - 3) ** 2 + (9 * z - 5) ** 2))
    x4 = -0.2 * np.exp(-(9 * x - 4) ** 2 - (9 * y - 7) ** 2 - (9 * z - 5) ** 2)
    return x1 + x2 + x3 + x4 + 1


def f7(x, y, z):
    return 1 / np.sqrt(1 + 2 * np.exp(-3 * (np.sqrt(x ** 2 + y ** 2 + z ** 2) - 6.7)))


# Ignores z
def f8(x, y, z):
    return x * (1 - x) * np.cos(4 * np.pi * x) * np.sin(4 * np.pi * y ** 2) ** 2 + 1


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))


















