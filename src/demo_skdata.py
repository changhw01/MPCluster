"""
=========================================================
Comparing different clustering algorithms on toy datasets
=========================================================

This example aims at showing characteristics of different
clustering algorithms on datasets that are "interesting"
but still in 2D. The last dataset is an example of a 'null'
situation for clustering: the data is homogeneous, and
there is no good clustering.

While these examples give some intuition about the algorithms,
this intuition might not apply to very high dimensional data.

The results could be improved by tweaking the parameters for
each clustering strategy, for instance setting the number of
clusters for the methods that needs this parameter
specified. Note that affinity propagation has a tendency to
create many clusters. Thus in this example its two parameters
(damping and per-point preference) were set to to mitigate this
behavior.
"""
#print __doc__

import time

import numpy as np
import pylab as pl

from sklearn import cluster, datasets
from sklearn.metrics import euclidean_distances
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler

from scipy import stats
import scipy.spatial.distance as scipy_dist
import matplotlib.pyplot as plt
from matplotlib import rc
import MPCA
reload(MPCA)


rc('text',usetex=True)
np.random.seed(0)

# Generate datasets. We choose the size big enough to see the scalability
# of the algorithms, but not too big to avoid too long running times
n_samples = 1500
noisy_circles = datasets.make_circles(n_samples=n_samples, factor=.5,
                                      noise=.05)
noisy_moons = datasets.make_moons(n_samples=n_samples, noise=.05)
blobs = datasets.make_blobs(n_samples=n_samples, random_state=8)
no_structure = np.random.rand(n_samples, 2), None

'''
colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
colors = np.hstack([colors] * 20)

pl.figure(figsize=(14, 9.5))
pl.subplots_adjust(left=.001, right=.999, bottom=.001, top=.96, wspace=.05,
                   hspace=.01)
'''

plot_num = 1
#for i_dataset, dataset in enumerate([noisy_circles, noisy_moons, blobs,
#                                     no_structure]):

for i_dataset, dataset in enumerate([blobs]):

    X, y = dataset
    # normalize dataset for easier parameter selection
    X = StandardScaler().fit_transform(X)
    
    # density estimate
    kernel = stats.gaussian_kde(X.T)
    dsty_est = kernel(X.T)
    I_sorting = np.argsort(dsty_est)
    I_sorting = I_sorting[::-1]

    dsty_est_sorted = dsty_est[I_sorting]
    X = X[I_sorting]

    #ax_dsty = fig.add_subplot(111)
    #ax_dsty.plot(dsty_est_sorted)

    dsty_thres = np.linspace(0.3,0.5,10)

    # distance matrix
    pD = scipy_dist.squareform(scipy_dist.pdist(X))

    dist_thres = np.linspace(0.05,0.7,10)


    # MPCA
    mpca = MPCA.MPCA(pD,dsty_est_sorted)
    mpca.fit(dsty_thres,dist_thres)
    
    pdiag = mpca.get_pdiagram(0)
    fig, ax = pdiag.show(max_size=mpca.size_levelset[0])
    ax.axis([0, pdiag.n_r+1, 0, pdiag.n_h+1])
    ax.set_xticks([2,4,6,8,10])
    ax.set_xticklabels(['$r_2$','$r_4$','$r_6$','$r_8$','$r_{10}$'])
    ax.set_yticks([2,4,6,8,10])
    ax.set_yticklabels(['$h_2$','$h_4$','$h_6$','$h_8$','$h_{10}$'])
    ax.grid(True, linestyle='-.', linewidth=1.2)
    fig.show()

