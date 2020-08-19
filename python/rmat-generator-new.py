"""
R-MAT generator for hypergraphs
-------------------------------
These algorithms are adapted from  Chakrabarti's work,
*R-MAT: A Recursive Model for Graph Mining*, https://www.cs.cmu.edu/~christos/PUBLICATIONS/siam04.pdf
    ```
    Chakrabarti, Deepayan & Zhan, Yiping & Faloutsos, Christos. (2004).
    R-MAT: A recursive model for graph mining. SIAM Proceedings Series.
    6. 10.1137/1.9781611972740.43.
    ```
The principle is to construct an incidence matrix based on two community distributions,
one for nodes and one for edges.
"""

import random as rand
import numpy as np
import itertools as it
from scipy.sparse import coo_matrix
import scipy
import scipy.io
from connected_components import  *
import time

def weight_distribution_independent(arr):
    """
    Create a conditional weight distribution for product of
    independent random variables and return pairing of coordinates
    and weights for recursive construction of incidence matrix

    Parameters
    ----------
    arr : list
        list of probability distributions, one for each dimension

    Returns
    -------
    coords, weights : list, list
        coords are a list of coordinates of dimension len(arr)
        weights are the probabilities of each set of
        coords
    """
    dims = [len(a) for a in arr]
    coords = list(it.product(*[range(a) for a in dims]))

    def f(c):
        return np.prod([arr[i][c[i]] for i in range(len(c))])

    weights = [f(c) for c in coords]
    return coords, weights


def _parr(arr, cp):
    def f(x):
        for idx, v in enumerate(cp):
            if x < v:
                return arr[idx]
    return f


def rmat_hypergraph(nrows, ncols, p, res, nv, return_matrix=True):
    """
    Generate an incidence matrix for a hypergraph using rmat
    algorithm with nrows x ncols seed matrix

    Parameters
    ----------
    nrows : int
        number of rows = number of nodes
    ncols : int
        number of columns = number of edges
    p : arr
        arr of nrows*ncols float values, probability distribution
    res : int
        Number of iterations before a cell index is resolved
    nv : int
        average number of nodes per edges

    Returns
    -------
    ndx : arr
        integer index of row of filled cell
    edx : arr
        integer index of column of filled cell
    """
    coords = np.array(list(it.product(range(nrows), range(ncols)))).transpose()
    assert len(p) == len(coords[0]), "size mismatch of seed matrix to probability"
    p = p/np.sum(p)

    args = np.argsort(p)[::-1]
    ii = coords[0][args]
    jj = coords[1][args]
    p = np.array(p)[args]
    cp = np.cumsum(p)
    nc = ncols**res * nv

    pi = np.frompyfunc(_parr(ii, cp), 1, 1)
    pj = np.frompyfunc(_parr(jj, cp), 1, 1)
    ndx = np.zeros(nc, dtype=int)
    edx = np.zeros(nc, dtype=int)

    for idx in range(res):
        rdx = np.random.rand(nc)
        iidx = pi(rdx)
        jjdx = pj(rdx)
        ndx = ndx + np.power(nrows, idx) * iidx
        edx = edx + np.power(ncols, idx) * jjdx
    if return_matrix:
        d = np.ones(nc, dtype=int)
        m = coo_matrix((d, (edx, ndx)))
        return edx, ndx, m
    else:
        return edx, ndx


nrows = 5
ncols = 3
nres  = 5
nv    = 7

p = np.array(
    [
        [0.5, 0.1, 0.1],
        [0.4, 0.0, 0.3],
        [0.1, 0.0, 0.0],
        [0.0, 0.2, 0.0],
        [0.0, 0.1, 0.3]
    ]
)

for nv in range(4,44,4):
  for nres in range(3,10):
    edx, ndx, c = rmat_hypergraph(nrows, ncols, p.flatten(), nres, nv, True)

    d = c.tocsr();
    c = d.tocoo();

    t = time.perf_counter()
    b = comps(c.col, c.row) 
    t = time.perf_counter()-t

    ncomps = len(b)

    fname = 'hrmat_%d_%d_%d_%d' % (nrows, ncols, nres, nv)

    comment_string = '''

Hypergraph generated by rmat-generator.py

  File name: %s
  Stored in COO format: 
      first column is hyperedge number
      second column is hypernode number
      the data has been sorted and duplicates removed

  nrows = %d
  ncols = %d
  p = 
  %s
  res = %d
  nv = %d

  #components = %d
  python time to compute = %f

''' % (fname, nrows, ncols, np.array2string(p), nres, nv, ncomps, t)

    print(fname)
    scipy.io.mmwrite(fname, c, comment=comment_string, field='pattern', precision=32, symmetry='general')




