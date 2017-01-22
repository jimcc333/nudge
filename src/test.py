import unittest
import os

from dbase import DBase
from objects import PathNaming
from library import Library


class Tester(unittest.TestCase):
    os.system('cls' if os.name == 'nt' else 'clear')

    paths = PathNaming(os.name, 'notARealPathNoWay')

    def test_database(self):

        self.assertEqual(self.paths.base_input, 'basecase.py')


    def test_second(self):
        with self.assertRaises(RuntimeError):
            database = DBase(self.paths)




unittest.main()
