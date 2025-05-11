import networkx as nx
from scipy.special import binom
from itertools import combinations
from math import dist
import numpy as np
import re
from scipy.spatial.distance import cdist
from sklearn.neighbors import kneighbors_graph

import numpy as np
import networkx as nx
from sklearn.neighbors import kneighbors_graph

import numpy as np
import networkx as nx

def knn_graph(W, k):
    """Filter symmetric adjacency matrix to keep k strongest neighbors per node, return NetworkX graph."""
    W = W.copy()
    n = W.shape[0]
    
    # Set self-connections to -inf to avoid selecting them
    np.fill_diagonal(W, -np.inf)
    
    # For each row, find indices of top k entries
    rows = np.arange(n)[:, None]
    top_k_idx = np.argpartition(-W, kth=k, axis=1)[:, :k]
    
    # Create a new zero matrix
    W_k = np.zeros_like(W)

    # Populate top-k entries per row
    W_k[rows, top_k_idx] = W[rows, top_k_idx]

    # Make symmetric: max(W_k, W_k.T) to preserve strongest edges in either direction
    W_k = np.maximum(W_k, W_k.T)
    
    # Create a NetworkX graph from the filtered adjacency matrix
    G = nx.from_numpy_array(W_k)
    
    return G



# Global neighborhood dictionary
cdef dict Neigh

def parse_graphml_to_networkx(filename):
    file = open(filename)
    Graph = nx.Graph()
    Dict = {}
    counter = -1
    weight_key = None
    key_info = "None"
    quoted_strings = ""

    for line in file:
        if "<key id" in line:
            quoted_strings = re.findall(r'"(.*?)"', line)
        if "weight" in quoted_strings:
            weight_key = quoted_strings[0]
            key_info = "<data key=" + '"' + weight_key + '"'
        if "<node id" in line:
            quoted_strings = re.findall(r'"(.*?)"', line)
            node = quoted_strings[0]
            counter += 1
            Dict[node] = counter
        if "<edge id" in line:
            quoted_strings = re.findall(r'"(.*?)"', line)
            Id, x, y = quoted_strings
            if x != y:
                Graph.add_edge(*[Dict[x], Dict[y]])
        if key_info in line:
            quoted_strings = re.findall(r'>(.*?)<', line)
            if x != y:
                Graph[Dict[x]][Dict[y]]["weight"] = float(quoted_strings[0])

    Graph.remove_edges_from(nx.selfloop_edges(Graph))
    Dict = {Dict[i]: i for i in Dict.keys()}
    return Graph, Dict    

def filter_graph(Obj, cutoff=np.inf,n_neighbors=None):
    if isinstance(Obj, nx.Graph):
        L = list(Obj.nodes())
        mapping = {i: L[i] for i in range(len(L))}
        M = nx.to_numpy_array(Obj, weight="weight")
        M = np.where(M <= cutoff, M, 0)
        H = nx.from_numpy_array(M)
        if n_neighbors!=None:
            H=knn_graph(nx.to_numpy_array(H),k=n_neighbors)
        return H, mapping
    elif isinstance(Obj, dict):
        L = list(Obj.keys())
        mapping = {i: L[i] for i in range(len(L))}
        points = np.array(list(Obj.values()))
        M = cdist(points, points)
        M = np.where(M <= cutoff, M, 0)
        H = nx.from_numpy_array(M)
        if n_neighbors!=None:
            H=knn_graph(nx.to_numpy_array(H),k=n_neighbors)
        return H, mapping
    elif isinstance(Obj, np.ndarray):
        mapping = {i: i for i in range(len(Obj))}
        M = np.where(Obj <= cutoff, Obj, 0)
        H = nx.from_numpy_array(M)
        if n_neighbors!=None:
            H=knn_graph(nx.to_numpy_array(H),k=n_neighbors)
        return H, mapping
    elif isinstance(Obj, str):
        H, mapping = parse_graphml_to_networkx(Obj)
        H = nx.to_numpy_array(H)
        H = np.where(H <= cutoff, H, 0)
        H = nx.from_numpy_array(H)
        if n_neighbors!=None:
            H=knn_graph(nx.to_numpy_array(H),k=n_neighbors)
        return H, mapping
    else:
        raise ValueError(f"Object type ({type(Obj)}) is not allowed to provide a network input: Object types allowed: nx.Graph, dict, np.array and graphml (str)")

# Manual set intersection counter for speed
cdef int intersect_size(set a, set b):
    cdef int count = 0
    for x in a:
        if x in b:
            count += 1
    return count

cpdef int ricci_cell(tuple c):
    cdef int d = len(c)
    cdef int H = -1
    cdef int N = 0
    cdef int i
    cdef set temp

    # Compute intersection size of all neighbors in c
    temp = Neigh[c[0]]
    for i in range(1, d):
        temp = temp & Neigh[c[i]]
    H = len(temp)

    # For each (d-1)-face, compute intersection size
    for face in combinations(c, d - 1):
        temp = Neigh[face[0]]
        for i in range(1, d - 1):
            temp = temp & Neigh[face[i]]
        N += len(temp)

    return (d + 1) * H + 2 * d - N

def compute_FRC(object D, float cutoff=np.inf,n_neighbors=None, int max_dim=1, mapped_nodes=False):
    """Computes the Forman-Ricci Curvature (FRC) up to dimension dim from the object provided.
   
     Input: D, an object that can be:

    - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
    - A nx.Graph object. The edge's attibutes `weight` will be considered.
    - A np.ndarray object that represents a symmetric matrix of float numbers
    cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
    (if D is a nx.Graph) of float values (if D is a np.matrix);
    -n_neighbors: Interger, the number of the nearest neighbours to consider. If None, it is not applied.
    - A string, the file name of the graphml/xml file.
    - mapped_nodes: True or False - if the mapping of the original node labels are kept

    dim: integer, the maximum simplex dimension allowed to the computation.

    Output:
    a dictionary whose keys are the dimensions and the values are dictionaries with the FRC for each cell."""
    global Neigh
    G, mapping = filter_graph(D, cutoff,n_neighbors)

    cdef int i, n, d, x, y, old_node,w
    cdef tuple old_cell
    cdef dict N = {d: {n: set() for n in G.nodes()} for d in range(1, max_dim + 1)}
    Neigh = {i: set(G.neighbors(i)) for i in G.nodes()}
    cdef dict F = {i: {} for i in range(1, max_dim + 1)}
    cdef tuple e
    cdef list E = list(G.edges())

    def func(tuple cell):
        d = len(cell)
        f = ricci_cell(cell)
        F[d - 1][cell] = f
       

        if d < max_dim + 1:
            if d == 2:
                x, y = cell
                N[d][x].update({y})
                N[d][y].update({x})
            else:
                old_cell = cell[:-1]
                old_node = cell[-1]
                for x in old_cell:
                    N[d][x].update({old_node})

            W = set.intersection(*[N[d][node] for node in cell])
            for w in W:
                new_cell = (*cell, w)
                func(new_cell)

    for e in E:
        func(e)

    if mapped_nodes:
        F = {d: {tuple(mapping[j] for j in cell): F[d][cell] for cell in F[d]} for d in F.keys()}

    return F

def compute_FRC_node(object D, float cutoff=np.inf,n_neighbors=None, int max_dim=1, mapped_nodes=False):
    """Computes the Forman-Ricci Curvature (FRC) for the nodes (up to dimension dim from the object provided).
   
    Input: D, an object that can be:

    - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
    - A nx.Graph object. The edge's attibutes `weight` will be considered.
    - A np.ndarray object that represents a symmetric matrix of float numbers
    cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
    (if D is a nx.Graph) of float values (if D is a np.matrix);
    - A string, the file name of the graphml/xml file.
     -n_neighbors: Interger, the number of the nearest neighbours to consider. If None, it is not applied.
    dim: integer, the maximum simplex dimension allowed to the computation.
     - mapped_nodes: True or False - if the mapping of the original node labels are kept

    Output:
    a dictionary whose keys are the dimensions and the values are dictionaries with the FRC for each node.
    If the node does not have curvature, the value information for the node will not be provided."""
    global Neigh
    G, mapping = filter_graph(D, cutoff,n_neighbors)

    cdef int i, n, d, x, y, old_node,node,w
    cdef tuple old_cell
    cdef dict N = {d: {n: set() for n in G.nodes()} for d in range(1, max_dim + 1)}
    Neigh = {i: set(G.neighbors(i)) for i in G.nodes()}
    cdef dict F = {i: {n:0 for n in G.nodes()} for i in range(1, max_dim + 1)}
    cdef dict C = {i: 0 for i in range(1, max_dim + 1)}
    cdef tuple e
    cdef list E = list(G.edges())
    cdef long counter

    def func(tuple cell):
        d = len(cell)
        f = ricci_cell(cell)
        C[d-1]+=1
        
        for node in cell:
            #print(d-1,cell,node)
            #print(d,node)
            F[d-1][node] += f
        
        if d < max_dim +1:
            if d == 2:
                x, y = cell
                N[d][x].update({y})
                N[d][y].update({x})
            else:
                old_cell = cell[:-1]
                old_node = cell[-1]
                for x in old_cell:
                    N[d][x].update({old_node})

            W = set.intersection(*[N[d][node] for node in cell])
            for w in W:
                new_cell = (*cell, w)
                func(new_cell)

    for e in E:
        func(e)
    for d in F.keys():
        for node in F[d].keys():
            try:
                F[d][node]=F[d][node]/((d+1)*C[d])
            except:
                next
        
    if mapped_nodes:
        F = {d: {mapping[node]:F[d][node] for node in F[d].keys()} for d in F.keys()}

    return F

def compute_average_FRC(object D, float cutoff=np.inf, n_neighbors=None, int max_dim=1):
    """Computes the average Forman-Ricci Curvature (FRC) up to dimension dim from the object provided.
    The average FRC is known as the sum of all local FRC divided by the number of d-cells.
    
    Input: D, an object that can be:

    - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
    - A nx.Graph object. The edge's attibutes `weight` will be considered.
    - A np.ndarray object that represents a symmetric matrix of float numbers
    cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
    (if D is a nx.Graph) of float values (if D is a np.matrix);
     -n_neighbors: Interger, the number of the nearest neighbours to consider. If None, it is not applied.
    - A string, the file name of the graphml/xml file.
     - mapped_nodes: True or False - if the mapping of the original node labels are kept

    dim: integer, the maximum simplex dimension allowed to the computation.
    
    Output:
    a dictionary whose keys are the dimensions and the values are the average FRC for respective dimension."""
    
    global Neigh
    G, mapping = filter_graph(D, cutoff,n_neighbors)

    cdef int i, n, d, x, y, old_node,node,w
    cdef tuple old_cell
    cdef dict N = {d: {n: set() for n in G.nodes()} for d in range(1, max_dim + 1)}
    Neigh = {i: set(G.neighbors(i)) for i in G.nodes()}
    cdef dict F = {i: 0.0 for i in range(1, max_dim + 1)}
    cdef dict C = {i: 0 for i in range(1, max_dim + 1)}
    cdef tuple e
    cdef list E = list(G.edges())
    cdef long counter

    def func(tuple cell):
        d = len(cell)
        f = ricci_cell(cell)
        C[d-1]+=1
        
      
        F[d-1] += f
        
        if d < max_dim +1:
            if d == 2:
                x, y = cell
                N[d][x].update({y})
                N[d][y].update({x})
            else:
                old_cell = cell[:-1]
                old_node = cell[-1]
                for x in old_cell:
                    N[d][x].update({old_node})

            W = set.intersection(*[N[d][node] for node in cell])
            for w in W:
                new_cell = (*cell, w)
                func(new_cell)

    for e in E:
        func(e)
    for d in F.keys():
      
        try:
            F[d]=F[d]/C[d]
        except:
            next
        
   # if mapped_nodes:
    #    F = {d: {mapping[node]:F[d][node] for node in F[d].keys()} for d in F.keys()}

    return F

def compute_FRC_node_frequency(object D, float cutoff=np.inf,n_neighbors=None, int max_dim=1, mapped_nodes=False):
    """Computes the frequency of Forman-Ricci Curvature (FRC) values up to dimension dim from the object provided.
  
   Input: D, an object that can be:

    - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
    - A nx.Graph object. The edge's attibutes `weight` will be considered.
    - A np.ndarray object that represents a symmetric matrix of float numbers
    cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
    (if D is a nx.Graph) of float values (if D is a np.matrix);
     -n_neighbors: Interger, the number of the nearest neighbours to consider. If None, it is not applied.
    - A string, the file name of the graphml/xml file.
     - mapped_nodes: True or False - if the mapping of the original node labels are kept
    dim: integer, the maximum simplex dimension allowed to the computation.

    Output:
    a dictionary whose keys are the dimensions and the values are the average FRC for respective dimension."""
    
   
    global Neigh
    G, mapping = filter_graph(D, cutoff,n_neighbors)

    cdef int i,j, n, d, x, y, old_node,w
    cdef int num_nodes=G.number_of_nodes()
    cdef tuple old_cell
    cdef dict N = {d: {n: set() for n in G.nodes()} for d in range(1, max_dim + 1)}
    Neigh = {i: set(G.neighbors(i)) for i in G.nodes()}
    cdef dict F = {i: {node:{j:0 for j in range(-num_nodes+2*(i+1),num_nodes+1)} for node in G.nodes()} for i in range(1, max_dim + 1)}
    cdef tuple e
    cdef list E = list(G.edges())

    def func(tuple cell):
        d = len(cell)
        f = ricci_cell(cell)
        for node in cell:
            F[d - 1][node][f] +=1       

        if d < max_dim + 1:
            if d == 2:
                x, y = cell
                N[d][x].update({y})
                N[d][y].update({x})
            else:
                old_cell = cell[:-1]
                old_node = cell[-1]
                for x in old_cell:
                    N[d][x].update({old_node})

            W = set.intersection(*[N[d][node] for node in cell])
            for w in W:
                new_cell = (*cell, w)
                func(new_cell)

    for e in E:
        func(e)

    if mapped_nodes:
        F = {d: {mapping[node]:{i:F[d][node][i] for i in F[d][node].keys()} for node in G.nodes()} for d in F.keys()}

    return F


def compute_FRC_frequency(object D, float cutoff=np.inf,n_neighbors=None, int max_dim=1):
    """Computes the frequency of Forman-Ricci Curvature (FRC) values up to dimension dim from the object provided.
  
   Input: D, an object that can be:

    - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
    - A nx.Graph object. The edge's attibutes `weight` will be considered.
    - A np.ndarray object that represents a symmetric matrix of float numbers
    cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
    (if D is a nx.Graph) of float values (if D is a np.matrix);
     -n_neighbors: Interger, the number of the nearest neighbours to consider. If None, it is not applied.
    - A string, the file name of the graphml/xml file.
    
     - mapped_nodes: True or False - if the mapping of the original node labels are kept

    dim: integer, the maximum simplex dimension allowed to the computation.

    Output:
    a dictionary whose keys are the dimensions and the values are the average FRC for respective dimension."""
    
  
    global Neigh
    G, mapping = filter_graph(D, cutoff,n_neighbors)

    cdef int i, n, d, x, y, old_node,w
    cdef int num_nodes=G.number_of_nodes()
    cdef tuple old_cell
    cdef dict N = {d: {n: set() for n in G.nodes()} for d in range(1, max_dim + 1)}
    Neigh = {i: set(G.neighbors(i)) for i in G.nodes()}
    cdef dict F = {i: {j:0 for j in range(-num_nodes+2*(i+1),num_nodes+1)} for i in range(1, max_dim + 1)}
    cdef tuple e
    cdef list E = list(G.edges())

    def func(tuple cell):
        d = len(cell)
        f = ricci_cell(cell)
       
        F[d - 1][f] +=1       

        if d < max_dim + 1:
            if d == 2:
                x, y = cell
                N[d][x].update({y})
                N[d][y].update({x})
            else:
                old_cell = cell[:-1]
                old_node = cell[-1]
                for x in old_cell:
                    N[d][x].update({old_node})

            W = set.intersection(*[N[d][node] for node in cell])
            for w in W:
                new_cell = (*cell, w)
                func(new_cell)

    for e in E:
        func(e)



    return F