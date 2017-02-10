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
    write_out = False

    # Check if inputs are passed to main
    if inputs is None:
        # Read inputs from file
        write_out = True
        lib = Library('blank', args[2], args[1], 0, False)  # The inputs are flipped here and this is fixed elsewhere...
        inputs = lib.inputs.xsgen

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
    output = output_method(all_inputs, inputs['placeholder_function'])

    if write_out:
        outputs = 'BUd = ' + str(output) + '\n' + 'NEUT_PROD = ' + str(0) + '\n' + 'NEUT_DEST = ' + str(0) + '\n'
        with open(args[2], 'w') as openfile:
            openfile.write(outputs)

    return output


def output_method(inputs, function):
    output = None
    if function == 'f2':
        output = f2(inputs[0], inputs[1], inputs[2])
    if function == 'f3':
        output = f3(inputs[0], inputs[1], inputs[2])
    if function == 'f4':
        output = f4(inputs[0], inputs[1], inputs[2])
    if function == 'f5':
        output = f5(inputs[0], inputs[1], inputs[2])
    if function == 'f6':
        output = f6(inputs[0], inputs[1], inputs[2])
    if function == 'f7':
        output = f7(inputs[0], inputs[1], inputs[2])
    if function == 'f8':
        output = f8(inputs[0], inputs[1], inputs[2])
    if function == 'f9':
        output = f9(inputs[0], inputs[1], inputs[2])

    if output is None:
        error_message = 'Placeholder xsgen does not recognize function type ' + str(function)
        raise RuntimeError(error_message)

    return output


def assign(ip_data, x):
    for key, value in ip_data.items():
        if value == None:
            ip_data[key] = x


""""
Some of the functions below have been taken from:
Testing Methods for 3D Scattered Data Interpolation
by Mira Bozzini and Milvia Rossini. Monograf´ıas de la Academia de Ciencias de Zaragoza. 20: 111–135, (2002).

"""


def f2(x, y, z):
    return (np.tanh(0.2 * x - 0.2 * z - 0.2 * y) + 1) / 2


# Good luck figuring this out
def f3(x, y, z):
    return np.cos(6 * z) * (1.25 + np.cos(5.4 * y)) / (6 + 6 * (3 * x - 1) ** 2) + f8(x, y, z) / 4 + f6(x, y, z) / 4


# Hump at center
def f4(x, y, z):
    return np.exp(-81 / 16 * ((x - 0.5) ** 2 + (y - 0.52) ** 2 + (z - 0.47) ** 2))


# Triple hump
def f5(x, y, z):
    A = np.sqrt(64 - 1 * ((x - 0.5) ** 2 + (y - 0.52) ** 2 + (z - 0.47) ** 2)) / 4
    return A + f2(x, y, z)/10 + f9(x, y, z) / 75


# Small hump at z=0.5
def f6(x, y, z):
    x1 = 0.75 * np.exp(-0.25 * ((9 * x - 2) ** 2 + (9 * y - 2) ** 2 + (9 * z - 2) ** 2))
    x2 = 0.75 * np.exp(- (9 * x + 1) ** 2 / 49 - (9 * y + 1) ** 2 / 10 - (9 * z + 1) ** 2 / 10)
    x3 = 0.50 * np.exp(-0.25 * ((9 * x - 7) ** 2 + (9 * y - 3) ** 2 + (9 * z - 5) ** 2))
    x4 = -0.2 * np.exp(-(9 * x - 4) ** 2 - (9 * y - 7) ** 2 - (9 * z - 5) ** 2)
    return x1 + x2 + x3 + x4 + 1


# 3D
def f7(x, y, z):
    return 1 / np.sqrt(1 + 2 * np.exp(-3 * (np.sqrt(x ** 2 + y ** 2 + z ** 2) - 6.7)))


# Ignores z
def f8(x, y, z):
    return x * (1 - x) * np.cos(4 * np.pi * x) * np.sin(4 * np.pi * y ** 2) ** 2 + 1


# 1 full 3 half peaks (taken from https://github.com/wgurecky/pCRTree/blob/master/tests/test_gradient_boost.py)
def f9(x, y, z):
    x *= 5
    y *= 10
    alpha, phi_ext = 0.7, 2 * np.pi * 0.5
    return (2 + alpha - 2 * np.cos(x) * np.cos(y) - alpha * np.cos(phi_ext - 2 * x)) / 12

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))


















