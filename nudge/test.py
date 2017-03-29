import unittest
import os
import shutil
import warnings

from nudge.dbase import DBase
from nudge.objects import PathNaming
from nudge.objects import xsgenParams
from nudge.library import Library


class Tester(unittest.TestCase):
    os.system('cls' if os.name == 'nt' else 'clear')
    test_folder = 'NUDGE_test_folder'

    path = PathNaming(os.name, 'notARealPathNoWay')
    xsgen = xsgenParams()

    # Create test database
    cwd = os.getcwd() + path.slash
    dbd = cwd + test_folder  # database directory
    try:
        os.mkdir(dbd)
    except FileExistsError:
        shutil.rmtree(cwd + path.slash + test_folder)
        os.mkdir(dbd)

    #   Add inputs and basecase
    ip_file = 'fuel_density 0.5' + os.linesep + 'clad_density 0.5' + os.linesep
    base_file = 'reactor = ' + os.linesep + 'fuel_density = 0.5' + os.linesep + 'clad_density = 0.5' + os.linesep +  \
                'cool_density = 0.5' + os.linesep + 'enrichment = 0.5' + os.linesep + 'unit_cell_pitch = 0.5' + \
                os.linesep + 'flux = 0.5' + os.linesep + \
                'k_particles = 0.5' + os.linesep + 'unit_cell_height = 0.5' + \
                os.linesep + 'void_cell_radius = 0.5' + os.linesep + 'fuel_cell_radius = 0.5' + os.linesep + \
                'clad_cell_radius = 0.5' + os.linesep
    with open(dbd + path.slash + path.dbase_input, 'w') as openfile:
        openfile.write(ip_file)
    openfile.close()
    with open(dbd + path.slash + path.base_input, 'w') as openfile:
        openfile.write(base_file)
    openfile.close()

    # Test database
    db_path = PathNaming(os.name, dbd + path.slash)
    database = DBase(db_path)

    # Test functions
    def test_2_PathNaming(self):
        self.assertEqual(self.path.base_input, 'basecase.py')

    def test_2_xsgen_inputs(self):
        self.assertEqual(self.xsgen.xsgen['fuel_cell_radius'], None)

    def test_2_xsgen_counter(self):
        self.assertEqual(self.xsgen.defined_count(), 0)
        self.xsgen.xsgen['clad_density'] = 1
        self.assertEqual(self.xsgen.defined_count(), 1)

    def test_3_database_add(self):
            self.assertIsNone(self.database.initial_exploration(False))

    def test_4_database_calcs(self):
            self.assertIsNone(self.database.voronoi())
            self.assertIsNone(self.database.estimate_error())

    def test_9_cleanup(self):
        shutil.rmtree(self.cwd + self.path.slash + self.test_folder)
        pass


unittest.main()
