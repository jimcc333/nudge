#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(name='nudge',
      version='0.0.1',
      description='Nuclear data generation software.',
      author='Cem Bagdatlioglu',
      packages=find_packages(),
      test_suite="tests",
      install_requires=['numpy>=1.8.0', 'scipy', 'matplotlib'],
      package_data={'': ['*.txt']},
      license='GPL-3',
      author_email='cem@cem-VirtualB',
      entry_points={
          'console_scripts': [
              'nudge=nudge.nudge:run_nudge'
          ]
      }
)
