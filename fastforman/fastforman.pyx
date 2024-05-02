import numpy as np
from itertools import combinations
from math import dist
import gudhi as gd
import networkx as nx



def cliques_gudhi(data,double f,int d):
    """An adaptation of Gudhi's algorithm for computing the cliques from a poitnt cloud data.
    Input: data - A list of n-dimensional points;
    f - the filtration value for the maximum distance between points;
    d - the maximum clique dimension allowed.
    Output: All cliques um to dimension d and diameter f."""

    cdef tuple i
    skeleton = gd.RipsComplex(points = data, max_edge_length = f)
    Rips_simplex_tree_sample = skeleton.create_simplex_tree(max_dimension = d)
    rips_generator = Rips_simplex_tree_sample.get_filtration()
      
    return (i[0] for i in rips_generator if len(i[0])>1)
   


cpdef int ricci_cell(list c,dict Neigh):

    cdef int d=len(c)
    cdef int j
    
    
    cdef int H=len(set.intersection(*[Neigh[j] for j in c]))
    cdef tuple i
      
    cdef list l=[i for i in combinations(c,d-1)]

    cdef tuple k
   

    cdef int N=sum([len(set.intersection(*[Neigh[j] for j in k])) for k in l])

   
    cdef int ricci=(d+1)*H+2*d-N    
    return ricci


cpdef dict FRC(L, dict Neigh):
    cdef list c
    cdef int ricci
    cdef tuple c0
    cdef int d

    forman=dict()
    for c in L:
        d=len(c)-1
        ricci=ricci_cell(c,Neigh)
        c0=tuple(c)
        try:
            forman[d][c0]=ricci
        except:
            forman[d]={}
            forman[d][c0]=ricci
          
    return forman


cpdef dict frequency(L, dict Neigh, int max_dim,int n):
    cdef list c
    cdef int ricci,i,k
    cdef tuple c0
  

    forman={i:{} for i in range(1,max_dim+1)}
    for c in L:
        d=len(c)-1
        ricci=ricci_cell(c,Neigh)
        #c0=tuple(c)
        try:
            forman[d][ricci]+=1
        except:
              
            forman[d][ricci]=0
            forman[d][ricci]+=1
    
    max_f=n
    for i in range(1,max_dim+1):
        min_f=-n+2*(i+1)
        for k in range(min_f,max_f+1):
            try:
                forman[i][k]
            except:
                forman[i][k]=0
        forman[i]= dict(sorted(forman[i].items()))       
    
                
           
    return forman


cpdef  dict compute_avg(L, dict Neigh,int max_dim):
    cdef list c
    cdef int ricci
    cdef int d
    cdef dict output={d:0 for d in range(1,max_dim+1)}
    cdef dict counter={d:0 for d in range(1,max_dim+1)}
    

    for c in L:
        d=len(c)-1
        ricci=ricci_cell(c,Neigh)
        counter[d]+=1
        output[d]+=ricci
    for d in range(1,max_dim+1):
        try:
            output[d]=output[d]/counter[d]
        except:
            output[d]=0.0
    return output

cpdef dict dict2neigh(dict D,double cutoff):
    cdef int i,x,y
    cdef double val
    cdef dict Neigh={i:set() for i in D.keys()}
    cdef tuple pair
    for pair in combinations(D.keys(),2):
        x,y=pair
        val=dist(D[x],D[y])
        if val <=cutoff:
            Neigh[x].add(y)
            Neigh[y].add(x)
    return Neigh

def compute_FRC(Obj,dim):
    """Computes the Forman-Ricci Curvature (FRC) up to dimension dim from the object provided.
    Input: 
    Obj: nx.Graph or tuple. If a nx.graph, returns the average FRC from the graph provided. Weights of nodes and edges are complitely ignored.
    if a tuple, the first entry has to be dictionary of point cloud data and the second entry is the cutoff distance (float)
    dim: the maximum dimension for computing FRC.

    Output:
    a dictionary whose keys are the dimensions and the values are dictionaries with the FRC for each cell."""
    if isinstance(Obj,nx.Graph):
        Neigh={i:set(Obj.neighbors(i)) for i in Obj.nodes()}
        D={i:[0] for i in Obj.nodes()}
        cutoff=0
        C=cliques_gudhi(D.values(),cutoff,dim)
        return FRC(C,Neigh)
    else:
        
        D,cutoff=Obj
        Neigh=dict2neigh(D,cutoff)
        C=cliques_gudhi(D.values(),cutoff,dim)
        return FRC(C,Neigh)


def compute_average_FRC(Obj,dim):
    """Computes the average Forman-Ricci Curvature (FRC) up to dimension dim from the object provided.
    The average FRC is known as the sum of all local FRC divided by the number of d-cells.
    Input: 
    Obj: nx.Graph or tuple. If a nx.graph, returns the average FRC from the graph provided. Weights of nodes and edges are complitely ignored.
    if a tuple, the first entry has to be dictionary of point cloud data and the second entry is the cutoff distance (float)
    dim: the maximum dimension for computing FRC.

    Output:
    a dictionary whose keys are the dimensions and the values are the average FRC for respective dimension."""
    if isinstance(Obj,nx.Graph):
        Neigh={i:set(Obj.neighbors(i)) for i in Obj.nodes()}
        D={i:[0] for i in Obj.nodes()}
        cutoff=0
        C=cliques_gudhi(D.values(),cutoff,dim)
        return compute_avg(C,Neigh,dim)
    else:
        
        D,cutoff=Obj
        Neigh=dict2neigh(D,cutoff)
        C=cliques_gudhi(D.values(),cutoff,dim)
        return compute_avg(C,Neigh,dim)

def compute_FRC_frequency(Obj,dim):
    """Computes the frequency of Forman-Ricci Curvature (FRC) values up to dimension dim from the object provided.
    Input: 
    Obj: nx.Graph or tuple. If a nx.graph, returns the average FRC from the graph provided. Weights of nodes and edges are complitely ignored.
    if a tuple, the first entry has to be dictionary of point cloud data and the second entry is the cutoff distance (float)
    dim: the maximum dimension for computing FRC.

    Output:
    a dictionary whose keys are the dimensions and the values are the average FRC for respective dimension."""
    if isinstance(Obj,nx.Graph):
        Neigh={i:set(Obj.neighbors(i)) for i in Obj.nodes()}
        n=Obj.number_of_nodes()
        D={i:[0] for i in Obj.nodes()}
        cutoff=0
        C=cliques_gudhi(D.values(),cutoff,dim)
        return frequency(C,Neigh,dim,n)
    else:
        
        D,cutoff=Obj
        Neigh=dict2neigh(D,cutoff)
        n=len(D)
        C=cliques_gudhi(D.values(),cutoff,dim)
        return frequency(C,Neigh,dim,n)
