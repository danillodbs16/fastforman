import numpy as np
from itertools import combinations,tee
from math import dist
import gudhi as gd
import networkx as nx
import re



def filter_graph(G,cutoff):
    for e in G.edges():
        x,y=e
        try:
            w=G[x][y]["weight"]
        except:
            w=1.
        if w>cutoff:
            G.remove_edge(*e)

def cliques_networkx(G,k):
    """An adaptation of  netowkrx algorithm for finding cliques up to fixed dimensinon.
    Inputs: G, a networkx simple undirected network.
        k, the maximum simplex dimension.
    Output: a generator with all the clique with size at most k+1."""

    cdef int i
    cdef tuple e,c,j0
    if k==1:
        C=(e for e in G.edges())
        for e in C:
            yield list(e)
    else:        
        Cl = (c for c in nx.find_cliques(G))
        C = (tuple(sorted(c)) for c in Cl)
        C=tee(C,k+1)
        for i in range(k+1):
            K=( j0 for j0 in set(c for mc in C[i] for c in combinations(mc, i+1)) )
            for c in K:
                if len(c)>1:
                    yield list(c)


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

def mapping_graph(G):
    H=nx.Graph()
    Nodes=list(G.nodes())
    mapping={Nodes[i]:i for i in range(len(Nodes))}
    H.add_nodes_from([i for i in range(len(Nodes))])
    for e in G.edges():
        x,y=e
        mx,my=mapping[x],mapping[y]
        H.add_edge(mx,my)
        try:
            H[mx][my]["weight"]=G[x][y]["weight"]
        except:
            next
    mapping={mapping[i]:i for i in mapping.keys()}
    return H,mapping        


def parse_graphml_to_networkx(filename):
    """Parses a netowrkx undirected graph from the graphml/xml file provided.
    Input: filename: a string of the file name.
    Output: A netowrkx undirected graph with numerical enumeration of nodes and a dictionary for mapping the original nodes labels"""
   
    file=open(filename)
    Graph=nx.Graph()
  
    Dict={}
    counter=-1
    weight_key=None
    key_info="None"
    quoted_strings=""
    
   
    for line in file:
        
        if "<key id" in line:
            quoted_strings = re.findall(r'"(.*?)"', line)
        if "weight" in quoted_strings:
            weight_key=quoted_strings[0]
            key_info="<data key="+'"'+weight_key+'"'
        #print(line)
        if "<node id" in line:
            quoted_strings = re.findall(r'"(.*?)"', line)
            node=quoted_strings[0]
          
            counter+=1
            Dict[node]=counter
        
    
        if "<edge id" in line:
            quoted_strings = re.findall(r'"(.*?)"', line)
            Id,x,y=quoted_strings
            if x!=y:
                Graph.add_edge(*[Dict[x],Dict[y]])
           # print(quoted_strings)
        if key_info in line:
            quoted_strings = re.findall(r'>(.*?)<', line)
            if x!=y:
                Graph[Dict[x]][Dict[y]]["weight"]=float(quoted_strings[0])
                
    Graph.remove_edges_from(nx.selfloop_edges(Graph))
    Dict={Dict[i]:i for i in Dict.keys()}
    return Graph,Dict    
    
        


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


cpdef dict FRC_node(L, dict Neigh):
    cdef list c
    cdef int ricci
    cdef tuple c0
    cdef int d,n
    cdef  dict counter=dict()
    cdef dict forman=dict()
    for c in L:
        d=len(c)-1
        ricci=ricci_cell(c,Neigh)

        try:
            counter[d]+=1
        except:
            counter[d]=0
            counter[d]+=1

        try:
            forman[d]
        except:
            forman[d]={}    

                
        for n in c:
            try:
                forman[d][n]+=ricci
            except:
                forman[d][n]=0
                forman[d][n]+=ricci

    for d in forman.keys():
        for n in forman[d].keys():
            try:
                forman[d][n]=forman[d][n]/((d+1)*counter[d])
            except:
                next    

          
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

cpdef dict frequency_node(L, dict Neigh, int max_dim,int n):
    cdef list c
    cdef int ricci,i,k
    cdef tuple c0
    cdef int node
  

    forman={i:{} for i in range(1,max_dim+1)}
    for c in L:
        d=len(c)-1
        ricci=ricci_cell(c,Neigh)
       
        for node in c:
            try:
                forman[d][node][ricci]+=1
            except:
                try:
                    forman[d][node][ricci]=0
                    forman[d][node][ricci]=+1
                except:
                    forman[d][node]={}
                    forman[d][node][ricci]=0
                    forman[d][node][ricci]+=1

              
    max_f=n
    for i in range(1,max_dim+1):
        min_f=-n+2*(i+1)
        for node in range(n):
            for ricci in range(min_f,max_f+1):
                try:
                    forman[i][node][ricci]

                except:
                    try:
                        forman[i][node][ricci]=0

                    except:
                        try:
                            forman[i][node]={}
                            forman[i][node][ricci]=0
                        except:
                            forman[i]={}
                            forman[i][node]={}
                            forman[i][node][ricci]=0 

        
            forman[i][node]= dict(sorted(forman[i][node].items()))                   

         
    
                
           
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

def mapping_dict(D):
    keys=list(D.keys())

    mapping={i:keys[i] for i in range(len(keys))}
    Dicto=dict(tuple(zip([i for i in range(len(D))],list(D.values()))))
    return Dicto,mapping

def compute_FRC(D,cutoff,dim):
    """Computes the Forman-Ricci Curvature (FRC) up to dimension dim from the object provided.
   
   Input: D, an object that can be:

    - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
    - A nx.Graph object. The edge's attibutes `weight` will be considered.
    - A np.ndarray object that represents a symmetric matrix of float numbers
    cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
    (if D is a nx.Graph) of float values (if D is a np.matrix);
    - A string, the file name of the graphml/xml file.

    dim: integer, the maximum simplex dimension allowed to the computation.

    Output:
    a dictionary whose keys are the dimensions and the values are dictionaries with the FRC for each cell."""

    def map_tuple(t, mapping):
        return tuple(mapping.get(i, i) for i in t)
    
    if isinstance(D, np.ndarray):
       # n=len(D)
        M=np.where(D<=cutoff,D,0)
        M=nx.from_numpy_array(M)
        Neigh={i:set(M.neighbors(i)) for i in M.nodes()}
        C=cliques_networkx(M,dim)
        F=FRC(C,Neigh)

        for d in range(1,dim+1):
            if d not in F.keys():
                F[d]={}

    elif isinstance(D,nx.Graph):
        G,Dicto=mapping_graph(D)
        G.remove_edges_from(nx.selfloop_edges(G))
        filter_graph(G,cutoff)
        #n=D.number_of_nodes()
        #M=nx.to_numpy_array(D)
        #M=np.where(M<=cutoff,M,0)
        #M=nx.from_numpy_array(M)
        Neigh={i:set(G.neighbors(i)) for i in G.nodes()}
        C=cliques_networkx(G,dim)
        F=FRC(C,Neigh)

        for d in range(1,dim+1):
            if d not in F.keys():
                F[d]={}

        for d in F.keys():
            try:
                F[d]={map_tuple(k, Dicto): v for k, v in F[d].items()}
            except:
                next 

    

    elif isinstance(D,str):
        G,Dicto= parse_graphml_to_networkx(D)
        #n=G.number_of_nodes()
        filter_graph(G,cutoff)
        Neigh={i:set(G.neighbors(i)) for i in G.nodes()}
        C=cliques_networkx(G,dim)
        F=FRC(C,Neigh)

        for d in range(1,dim+1):
            if d not in F.keys():
                F[d]={}

        
        for d in F.keys():
            try:
                F[d]={map_tuple(k, Dicto): v for k, v in F[d].items()}
            except:
                next    
                    

        

    else:
       # n=len(D)
        Dicto,mapping=mapping_dict(D)
        Neigh=dict2neigh(Dicto,cutoff)
        C=cliques_gudhi(Dicto.values(),cutoff,dim)
        F=FRC(C,Neigh)
        for d in range(1,dim+1):
            if d not in F.keys():
                F[d]={}

        for d in F.keys():
            try:
                F[d]={map_tuple(k, mapping): v for k, v in F[d].items()}
            except:
                next  

    
   
    return F

cpdef dict fill_info(dict D,int n):
    Out=D.copy()
    cdef int k
    for k in range(n):
        try:
            Out[k]
        except:
            Out[k]=0
    return Out        

def compute_FRC_node(D,cutoff,dim):
    """Computes the Forman-Ricci Curvature (FRC) for the nodes (up to dimension dim from the object provided).
   
   Input: D, an object that can be:

    - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
    - A nx.Graph object. The edge's attibutes `weight` will be considered.
    - A np.ndarray object that represents a symmetric matrix of float numbers
    cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
    (if D is a nx.Graph) of float values (if D is a np.matrix);
    - A string, the file name of the graphml/xml file.
    dim: integer, the maximum simplex dimension allowed to the computation.

    Output:
    a dictionary whose keys are the dimensions and the values are dictionaries with the FRC for each node.
    If the node does not have curvature, the value information for the node will not be provided."""
    
    if isinstance(D, np.ndarray):
        n=len(D)
        M=np.where(D<=cutoff,D,0)
        M=nx.from_numpy_array(M)
        Neigh={i:set(M.neighbors(i)) for i in M.nodes()}
        C=cliques_networkx(M,dim)
        F=FRC_node(C,Neigh)

        for d in range(1,dim+1):
            try:
                F[d]
            except:
                F[d]={}
               
            for i in range(n):
                try:
                    F[d][i]
                except:
                    F[d][i]=np.nan
            F[d]=dict(sorted(F[d].items())) 


    elif isinstance(D,nx.Graph):
        G,Dicto=mapping_graph(D)
        G.remove_edges_from(nx.selfloop_edges(G))
        filter_graph(G,cutoff)
        n=G.number_of_nodes()
        #M=nx.to_numpy_array(D)
        #M=np.where(M<=cutoff,M,0)
        #M=nx.from_numpy_array(M)
        Neigh={i:set(G.neighbors(i)) for i in G.nodes()}
        C=cliques_networkx(G,dim)
        F=FRC_node(C,Neigh)

        for d in range(1,dim+1):
            try:
                F[d]
            except:
                F[d]={}
               
            for i in range(n):
                try:
                    F[d][i]
                except:
                    F[d][i]=np.nan
            F[d]=dict(sorted(F[d].items()))

        for d in F.keys():
            F[d]={Dicto.get(k, k): v for k, v in F[d].items()}  


    elif isinstance(D,str):
        G,Dicto= parse_graphml_to_networkx(D)
        n=G.number_of_nodes()
        filter_graph(G,cutoff)
        Neigh={i:set(G.neighbors(i)) for i in G.nodes()}
        C=cliques_networkx(G,dim)
        F=FRC_node(C,Neigh)


        for d in range(1,dim+1):
            try:
                F[d]
            except:
                F[d]={}
               
            for i in range(n):
                try:
                    F[d][i]
                except:
                    F[d][i]=np.nan
            F[d]=dict(sorted(F[d].items()))

        for d in F.keys():
            F[d]={Dicto.get(k, k): v for k, v in F[d].items()}      



    else:
        n=len(D)
        Dicto,mapping=mapping_dict(D)
        Neigh=dict2neigh(Dicto,cutoff)
        C=cliques_gudhi(Dicto.values(),cutoff,dim)
        F=FRC_node(C,Neigh)
        for d in range(1,dim+1):
            try:
                F[d]
            except:
                F[d]={}
               
            for i in range(n):
                try:
                    F[d][i]
                except:
                    F[d][i]=np.nan
            F[d]=dict(sorted(F[d].items()))

        for d in F.keys():
            F[d]={mapping.get(k, k): v for k, v in F[d].items()}      


    

   
    return F

def compute_average_FRC(D,cutoff,dim):
    """Computes the average Forman-Ricci Curvature (FRC) up to dimension dim from the object provided.
    The average FRC is known as the sum of all local FRC divided by the number of d-cells.
    
    Input: D, an object that can be:

    - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
    - A nx.Graph object. The edge's attibutes `weight` will be considered.
    - A np.ndarray object that represents a symmetric matrix of float numbers
    cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
    (if D is a nx.Graph) of float values (if D is a np.matrix);
    - A string, the file name of the graphml/xml file.

    dim: integer, the maximum simplex dimension allowed to the computation.
    
    Output:
    a dictionary whose keys are the dimensions and the values are the average FRC for respective dimension."""
    
    if isinstance(D, np.ndarray):
        M=np.where(D<=cutoff,D,0)
        M=nx.from_numpy_array(M)
        Neigh={i:set(M.neighbors(i)) for i in M.nodes()}
        C=cliques_networkx(M,dim)
        F=compute_avg(C,Neigh,dim)

    elif isinstance(D,nx.Graph):
        #M=nx.to_numpy_array(D)
        #M=np.where(M<=cutoff,M,0)
        M,Dicto=mapping_graph(D)
        M.remove_edges_from(nx.selfloop_edges(M))
        filter_graph(M,cutoff)
        Neigh={i:set(M.neighbors(i)) for i in M.nodes()}
        C=cliques_networkx(M,dim)
        F=compute_avg(C,Neigh,dim)

    elif isinstance(D,str):
        G,Dicto= parse_graphml_to_networkx(D)
        G.remove_edges_from(nx.selfloop_edges(G))
        n=G.number_of_nodes()
        filter_graph(G,cutoff)
        Neigh={i:set(G.neighbors(i)) for i in G.nodes()}
        C=cliques_networkx(G,dim)
        F=compute_avg(C,Neigh,dim)
        

    else:
        Dicto,mapping=mapping_dict(D)
        Neigh=dict2neigh(Dicto,cutoff)
        C=cliques_gudhi(Dicto.values(),cutoff,dim)
        F=compute_avg(C,Neigh,dim)
    
    
    
    return F

def compute_FRC_frequency(D,cutoff,dim):
    """Computes the frequency of Forman-Ricci Curvature (FRC) values up to dimension dim from the object provided.
  
   Input: D, an object that can be:

    - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
    - A nx.Graph object. The edge's attibutes `weight` will be considered.
    - A np.ndarray object that represents a symmetric matrix of float numbers
    cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
    (if D is a nx.Graph) of float values (if D is a np.matrix);
    - A string, the file name of the graphml/xml file.

    dim: integer, the maximum simplex dimension allowed to the computation.

    Output:
    a dictionary whose keys are the dimensions and the values are the average FRC for respective dimension."""
    
       
    if isinstance(D, np.ndarray):
        n=len(D)
        M=np.where(D<=cutoff,D,0)
        M=nx.from_numpy_array(M)
        Neigh={i:set(M.neighbors(i)) for i in M.nodes()}
        C=cliques_networkx(M,dim)
        F=frequency(C,Neigh,dim,n)
        for d in range(1,dim+1):
            if d not in F.keys():
                F[d]={}

    elif isinstance(D,nx.Graph):
        n=len(D)
        M,Dicto=mapping_graph(D)
        M.remove_edges_from(nx.selfloop_edges(M))
        filter_graph(M,cutoff)
        #M=nx.to_numpy_array(D)
        #M=np.where(M<=cutoff,M,0)
        #M=nx.from_numpy_array(M)
        Neigh={i:set(M.neighbors(i)) for i in M.nodes()}
        C=cliques_networkx(M,dim)
        F=frequency(C,Neigh,dim,n)
        for d in range(1,dim+1):
            if d not in F.keys():
                F[d]={}

    elif isinstance(D,str):
        G,Dicto= parse_graphml_to_networkx(D)
        G.remove_edges_from(nx.selfloop_edges(G))
        n=G.number_of_nodes()
        filter_graph(G,cutoff)
        Neigh={i:set(G.neighbors(i)) for i in G.nodes()}
        C=cliques_networkx(G,dim)
        F=frequency(C,Neigh,dim,n)
        for d in range(1,dim+1):
            if d not in F.keys():
                F[d]={}


    else:
        n=len(D)
        Dicto,mapping=mapping_dict(D)
        Neigh=dict2neigh(Dicto,cutoff)
        C=cliques_gudhi(Dicto.values(),cutoff,dim)
        F=frequency(C,Neigh,dim,n)


        for d in range(1,dim+1):
            if d not in F.keys():
                F[d]={}
    
    
    
    return F


def compute_FRC_node_frequency(D,cutoff,dim):
    """Computes the frequency of Forman-Ricci Curvature (FRC) values up to dimension dim from the object provided.
  
   Input: D, an object that can be:

    - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
    - A nx.Graph object. The edge's attibutes `weight` will be considered.
    - A np.ndarray object that represents a symmetric matrix of float numbers
    cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
    (if D is a nx.Graph) of float values (if D is a np.matrix);
    - A string, the file name of the graphml/xml file.

    dim: integer, the maximum simplex dimension allowed to the computation.

    Output:
    a dictionary whose keys are the dimensions and the values are also dictionaries with the FRC frequency per node"""
    
       
    if isinstance(D, np.ndarray):
        n=len(D)
        M=np.where(D<=cutoff,D,0)
        M=nx.from_numpy_array(M)
        Neigh={i:set(M.neighbors(i)) for i in M.nodes()}
        C=cliques_networkx(M,dim)
        F=frequency_node(C,Neigh,dim,n)

    elif isinstance(D,nx.Graph):
        n=len(D)
        G,Dicto=mapping_graph(D)
        G.remove_edges_from(nx.selfloop_edges(G))
        n=G.number_of_nodes()
        filter_graph(G,cutoff)
        #M=nx.to_numpy_array(D)
        #M=np.where(M<=cutoff,M,0)
        #M=nx.from_numpy_array(M)
        Neigh={i:set(G.neighbors(i)) for i in G.nodes()}
        C=cliques_networkx(G,dim)
        F=frequency_node(C,Neigh,dim,n)
        for d in F.keys():
            F[d]={Dicto.get(k, k): v for k, v in F[d].items()} 


    elif isinstance(D,str):
        G,Dicto= parse_graphml_to_networkx(D)
        G.remove_edges_from(nx.selfloop_edges(G))
        n=G.number_of_nodes()
        filter_graph(G,cutoff)
        Neigh={i:set(G.neighbors(i)) for i in G.nodes()}
        C=cliques_networkx(G,dim)
        F=frequency_node(C,Neigh,dim,n)
        for d in F.keys():
            F[d]={Dicto.get(k, k): v for k, v in F[d].items()} 


    else:
        n=len(D)
        Dicto,mapping=mapping_dict(D)
        Neigh=dict2neigh(Dicto,cutoff)
        C=cliques_gudhi(Dicto.values(),cutoff,dim)
        F=frequency_node(C,Neigh,dim,n)
        for d in F.keys():
            F[d]={mapping.get(k, k): v for k, v in F[d].items()} 


    
    
    
    return F
