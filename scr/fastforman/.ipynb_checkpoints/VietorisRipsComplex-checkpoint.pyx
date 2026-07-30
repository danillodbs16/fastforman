import numpy as np
from itertools import combinations
from scipy.special import binom
from math import dist as distance
from libc.math cimport pow
from sklearn.neighbors import NearestNeighbors
from time import time
import networkx as nx
from typing import Tuple


def knn_distance_matrix(X, k, metric='euclidean'):
    nbrs = NearestNeighbors(n_neighbors=k+1, metric=metric).fit(X)
    
    # Get the sparse distance matrix (k neighbors per row)
    knn_graph = nbrs.kneighbors_graph(X, mode='distance').toarray()
    
    return knn_graph

def cliques(G, int max_dim):
    cdef list boundary
    cdef tuple edge
    cdef dict N = {}
    cdef list E = sorted(
        [((u, v), data['weight']) for u, v, data in G.edges(data=True)],
        key=lambda x: x[1]
    )

    cdef tuple cell, new_cell
    cdef set node
    cdef frozenset b
    cdef set W
    cdef double w0  # Initial edge weight (floating point)

    def func(tuple cell, double w0)-> tuple[tuple,double]:
        cdef int d = len(cell)
        
        
        yield (cell, w0)

        if d < max_dim:
            boundary = [frozenset(b0) for b0 in combinations(cell, d - 1)]
            for b in boundary:
                node = set(cell) - set(b)
                try:
                    N[b].update(node)
                except:
                    N[b] = set()
                    N[b].update(node)

            W = set.intersection(*[N[b] for b in boundary])
            for w in W:
                new_cell = (*cell, w)
                yield from func(new_cell, w0)

    for e in E:
        edge,w0=e
        #cell = tuple(edge)  # Sort to ensure consistent order
        yield from func(edge,w0)


cpdef int delta_cython(bint val, int d):
    return (d + 1) if val else -1


def compute_local_FRC(M, int max_dim, float max_dist, int k_nearest=-1, int precision=2, metric="euclidean",verbose=False) -> tuple[dict,dict]: 
    """Computes the d-th Forman-Ricci curvature (FRC) per node for d in {1,...,max_dim} from the Vietoris-Rips (VR) complex generated from the filtration distance max_dist.
    Inputs:
 
    - M: np.array matrix of points;
    - max_dim: int, the maximum simplex dimension;
    - max_dist: float, the maximum cutoff distance for the VR filtration;
    - k_nearest: int, the number of nearest neighbours for neighbourhood filtration. If not provided, it will be understood as unlimited;
    - precision: int, the precision of decimal places for the distance cutoff;
    - metric: string, the metric space for computing pairwise distance, as in scipy.spatial.distance.cdist parameters;
    - verbose: str, if the benchmar/debugging mode is activated.

    Outputs:

    - Output1: dict, the total d-th FRC per node, for d in {1, ..., max_dim};
    - Output2: dict, the average d-th FRC per node, for d in {1, ..., max_dim}.
    """
    from itertools import combinations
    from math import pow
    import numpy as np
    import networkx as nx
    from sklearn.neighbors import NearestNeighbors
    from scipy.special import binom
    from time import time
    cdef int d,u,v
    cdef double f,w0
    cdef list stack
    cdef set H

    cdef int nnodes = M.shape[0]
    cdef float delta_dist = pow(10.0, -precision)
    cutoffs = np.round(np.arange(0, max_dist + delta_dist, delta_dist), decimals=precision)

    # Build distance matrix and graph
    if k_nearest == -1:
        D = knn_distance_matrix(M, nnodes - 1, metric=metric)
    else:
        D = knn_distance_matrix(M, k_nearest, metric=metric)

    D = np.round(D, precision)
    D = np.where(D <= max_dist, D, 0)
    G = nx.from_numpy_array(D)
    cdef list nodes = list(range(nnodes))

    # Preallocate structures
    cdef dict combs = {d: binom(nnodes, d + 1) for d in range(1, max_dim + 1)}
    cdef dict F = {d: 0.0 for d in range(1, max_dim + 1)}
    cdef dict C = {d: 0 for d in range(1, max_dim + 1)}
    cdef dict N_faces = {d: {} for d in range(1, max_dim + 1)}
    cdef dict nF = {d: {n: 0.0 for n in nodes} for d in range(1, max_dim + 1)}
    cdef dict L1 = {d: [] for d in range(1, max_dim + 1)}
    cdef dict  L2 = {d: [] for d in range(1, max_dim + 1)}
    cdef dict Aux_Dist = {d: 0.0 for d in range(1, max_dim + 1)} 
    cdef dict Aux_Line1= {d: 0.0 for d in range(1, max_dim + 1)} 
    cdef dict Aux_Line2 = {d: 0.0 for d in range(1, max_dim + 1)} 
    cdef list B,out1,out2
    cdef int neigh
    cdef dict headers = {
        d: (["cutoff", "density", "avg_forman"] + nodes,
            ["cutoff", "density", "total_forman"] + nodes)
        for d in range(1, max_dim + 1)
    }

    # Sort edges by distance
    cdef list edges = sorted([((u, v), data['weight']) for u, v, data in G.edges(data=True)], key=lambda x: x[1])
    N_clique = {}

    t0 = time()

    # Explicit stack for manual recursion
    for edge, w0 in edges:
        stack = [tuple(sorted(edge))]

        while stack:
            cell = stack.pop()
            d = len(cell) - 1
            if d > max_dim:
                continue

            # === Begin Ricci computation ===
            dist = round(w0, precision)
            if d >= 1:
                C[d] += 1
                B = [frozenset(b) for b in combinations(cell, d)]

                for b in B:
                    N_faces[d].setdefault(b, set()).update(set(cell) - set(b))

                H = set.intersection(*[N_faces[d][b] for b in B])
                neigh = sum(len(N_faces[d][b]) for b in B)
                f = (d + 2) * len(H) + 2 * (d + 1) - neigh

                for node in cell:
                    nF[d][node] += f / (d + 1)

                for b in B:
                    for a in N_faces[d][b] - set(cell):
                        delta_val = ((d + 1) if a in H else -1) / (d + 1)
                        for node in b:
                            nF[d][node] += delta_val
                        nF[d][a] += delta_val
                        f += delta_val * (d + 1)

                F[d] += f

                out1 = [dist, C[d] / combs[d], F[d] / C[d]] + [nF[d][i] / C[d] for i in nodes]
                out2 = [dist, C[d] / combs[d], F[d]] + [nF[d][i] for i in nodes]

                if Aux_Dist[d] != dist:
                    if Aux_Line1[d]:
                        L1[d].append(Aux_Line1[d])
                        L2[d].append(Aux_Line2[d])

                Aux_Dist[d] = dist
                Aux_Line1[d] = out1
                Aux_Line2[d] = out2
            # === End Ricci computation ===

            # Build higher-dimensional cliques
            if d < max_dim:
                boundary = [frozenset(b) for b in combinations(cell, d)]
                for b in boundary:
                    new_node = set(cell) - set(b)
                    N_clique.setdefault(b, set()).update(new_node)

                W = set.intersection(*[N_clique[b] for b in boundary])
                for w in W:
                    new_cell = tuple(sorted(set(cell) | {w}))
                    stack.append(new_cell)

    # Final push
    for d in range(1, max_dim + 1):
        if Aux_Line1[d]:
            L1[d].append(Aux_Line1[d])
            L2[d].append(Aux_Line2[d])

    # Post-processing
    def post_process(L, header):
        Aux = {d: {row[0]: dict(zip(header[d][0], row)) for row in L[d]} for d in range(1, max_dim + 1)}
        null_row = dict(zip(header[1][0], [np.nan] * len(header[1][0])))
        Output = {d: {cut: null_row.copy() for cut in cutoffs} for d in range(1, max_dim + 1)}

        for d in range(1, max_dim + 1):
            keys = sorted(Aux[d].keys())
            if not keys:
                continue
            last_cut = keys[0]
            for cut in cutoffs:
                if cut in Aux[d]:
                    Output[d][cut] = Aux[d][cut].copy()
                    last_cut = cut
                elif cut >= last_cut:
                    Output[d][cut] = Aux[d][last_cut].copy()
                    Output[d][cut]["cutoff"] = cut
                    Output[d][cut]["density"] = np.nan
        return Output

    t1 = time()
    Output1 = post_process(L1, headers)
    Output2 = post_process(L2, headers)
    t2 = time()

    if verbose:
        print(f"Total ricci computation time: {t1 - t0:.4f} seconds.")
        print(f"Total post-processing time: {t2 - t1:.4f} seconds.")
        print(f"Total computation time: {t2 - t0:.4f} seconds.")

    return Output1, Output2

