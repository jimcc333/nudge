from nudge.pca import NudgePCA
from matplotlib import pyplot as plt
from matplotlib.mlab import PCA as mlabPCA
import unittest
import os
import numpy as np
pwd_ = os.path.dirname(os.path.abspath(__file__))
dataDir = pwd_ + "/data/"


class TestPCA(unittest.TestCase):
    def test_pca_iris(self):
        # load the iris dataset
        iris = np.loadtxt(dataDir + "iris.csv", skiprows=1, usecols=(1,2,3,4,6), delimiter=',')
        # class response
        y = iris[:, -1].astype(int)
        # explanatory vars
        x = iris[:, 0:-1]

        # apply pca
        npca = NudgePCA(x)
        npca.plot_ranks('test_pca.png')

        print("Iris dataset eigen vales:")
        print(npca.eig_vals)
        print("Iris dataset eigen vectors:")
        print(npca.eig_vecs)
        print("Nudge fractioinal varience:")
        print(npca.frac_explained_var)

        # Check against mlabPCA for "validation"
        mpca = mlabPCA(x)
        print("Matplotlib fractioinal varience:")
        print(mpca.fracs)

        # check against expected result
        self.assertAlmostEqual(npca.frac_explained_var[0], mpca.fracs[0], delta=1e-4)
        self.assertAlmostEqual(npca.frac_explained_var[1], mpca.fracs[1], delta=1e-4)
        self.assertAlmostEqual(npca.frac_explained_var[2], mpca.fracs[2], delta=1e-4)
        self.assertAlmostEqual(npca.frac_explained_var[3], mpca.fracs[3], delta=1e-4)

        # check projecting matrix
        w = npca.pcw()
        print("Projection matrix:")
        print(w)

        # project original data onto pca principal axes
        x_transformed = npca.project(retain_frac_var=0.95)

        # plot transformed data
        plt.figure()
        plt.scatter(x_transformed[:, 0], x_transformed[:, 1])
        plt.savefig('project_test.png')
        plt.xlabel('pc 1')
        plt.ylabel('pc 2')
        plt.close()
