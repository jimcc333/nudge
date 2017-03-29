import numpy as np
from matplotlib import pyplot as plt


class NudgePCA(object):
    """!
    @brief Nudge PCA class.
    Extends functionality of matplotlib's PCA submodule.
        - Adds ability to supply weights alongside data
        - Adds ability to quickly compute a projection matrix
        which retains a given fraction of the total original data varience.
    """
    def __init__(self, data, weights=None):
        self._weights = weights
        self._data = data
        self._eig_pairs = None
        self._frac_explained_var = None
        self._cum_frac_explained_var = None
        self._pca()

    @property
    def eig_vals(self):
        """!
        @brief Convinience accessor to eigen values
        """
        pass
        return np.array([self._eig_pairs[i][0] \
               for i in range(len(self._eig_pairs))])

    @property
    def eig_vecs(self):
        """!
        @brief Convinience accessor to eigen vectors.
        The i_th column is an eigenvector corrosponding to
        the i_th eigenvalue returned by self.eig_vals
        """
        return np.array([self._eig_pairs[i][1] \
               for i in range(len(self._eig_pairs))]).T

    @property
    def eig_pairs(self):
        """!
        @brief Read only eig_pairs attribute
        Determined by self._pca() method
        """
        return self._eig_pairs

    @property
    def frac_explained_var(self):
        """!
        @brief Read only frac_explained_var attribute
        Determined by self._pca() method
        """
        return self._frac_explained_var

    @property
    def cum_frac_explained_var(self):
        """!
        @brief Read only cum_frac_explained_var attribute
        Determined by self._pca() method
        """
        return self._cum_frac_explained_var

    def pcw(self, retain_frac_var=0.95):
        """!
        @brief Computes a pca projection matrix.
        Maps from the original data space to a reduced dim space.
        @param retainFracVar <b>double</b> Desired fraction of explained varience to retain
        @param reducedDim <b>int</b> Target reduced data dimension
        @return <b>np_ndarray</b> PC projection matrix, W
        """
        retained_eig_vecs = []
        for i in range(len(self._eig_pairs)):
            retained_eig_vecs.append(self._eig_pairs[i][1].reshape(len(self._eig_pairs), 1))
            if self._cum_frac_explained_var[i] >= retain_frac_var:
                break
        pcW = np.hstack(retained_eig_vecs)
        return pcW

    def project(self, data_in=None, retain_frac_var=0.95):
        """!
        @brief Apply pca projection matrix to data.
        \f[
           y_{reduced} = W.T \cdot x
        \f]
        @param data_in  Optional custom data to project onto principle axes.
            Default is to project the original data onto principle axes.
        @return <b>np_ndarray</b> transformed output data
        """
        if data_in is None:
            data_in = self._data
        else:
            assert(data_in.shape[1] == self._data.shape[1])
        pcW = self.pcw(retain_frac_var)
        return data_in.dot(pcW)

    def plot_ranks(self, fig_name='pca_default.png'):
        """!
        @brief Plots relative explained varience of each
        principle component.
        """
        plt.figure()
        plt.bar(range(len(self._eig_pairs)),
                self._frac_explained_var,
                align='center',
                label='relative explained varience')
        plt.xlabel('Principal Components')
        plt.ylabel('Explained Varience Ratio')
        plt.savefig(fig_name)
        plt.close()

    def _pca(self):
        """!
        @brief Computes principal components of multivariate data set.
        Provides a measure of explained varience per principal component,
        the principal compenent directions and magnitudes.
        @param nd_data  <b>np_ndarray</b>  n-dim data array
        @param weights  <b>np_1darray</b>  optional data weights

        note: eigen_pairs are formated as (eig_val, eig_vector) tuples
        """
        # scale and shift data so that mean ==0 and var==1 in all columns
        std_mvdData = (self._data - self._data.mean(axis=0)) / self._data.std(axis=0)
        # compute cov matrix and eig vals
        std_cov = np.cov(std_mvdData.T, aweights=self._weights)
        eig_vals, eig_vecs = np.linalg.eig(std_cov)
        # Create sorted eigen-pairs
        eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:, i]) for i in range(len(eig_vals))]
        eig_pairs.sort(key=lambda x: x[0], reverse=True)
        eig_sum = np.sum(eig_vals)
        frac_explained_var = [(s / eig_sum) for s in sorted(eig_vals, reverse=True)]
        cum_frac_explained_var = np.cumsum(frac_explained_var)
        # Assign to read only attributes
        self._eig_pairs = eig_pairs
        self._frac_explained_var = frac_explained_var
        self._cum_frac_explained_var = cum_frac_explained_var

