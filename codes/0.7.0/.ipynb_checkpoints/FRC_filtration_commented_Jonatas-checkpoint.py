#!/usr/bin/env python
# coding: utf-8

# In[1]:


import gudhi as gd
import pandas as pd
import numpy as np
import networkx as nx
from time import time
from scipy.special import binom
import matplotlib.pyplot as plt
from itertools import combinations
from math import dist


# In[13]:

#Compute cliques from distace matrix. d is the maximum clique dimension. f is the maximum distance allowed.
def cliques_gudhi_matrix(distance_matrix,f,d=2):
  
    
    rips_complex = gd.RipsComplex(distance_matrix=distance_matrix, max_edge_length=f)
    Rips_simplex_tree_sample = rips_complex.create_simplex_tree(max_dimension = d)
    rips_generator = Rips_simplex_tree_sample.get_filtration()
 
    return (i for i in rips_generator if len(i[0])>1 and i[1]!=0)


#Same as before, but with the points in a list as input instead.
def cliques_gudhi(data,f,d=2):
#    
    skeleton = gd.RipsComplex(points = data, max_edge_length = f)
    Rips_simplex_tree_sample = skeleton.create_simplex_tree(max_dimension = d)
    rips_generator = Rips_simplex_tree_sample.get_filtration()

    return (i for i in rips_generator if len(i[0])>1)
    


#A function to recover point cloud data from a specific dataframe. YOu will not need to use it.
def recover_data(dataframe):
    
    header=list(dataframe)
    rows=len(dataframe)
    T={h:[] for h in header}
    for column in header:
        for row in range(rows):
            t=dataframe[column][row][1:-1].split()
           # print(t)
            t=np.array([float(i) for i in t])
            #print(t)
            T[column].append(t)
    return T


#THe delta function as in iur manuscrit. 
def delta(val, d):
    if val==True:
        return (d+1)
    else:
        return -1


#Computes FRC from a distance matrix (or measure matrix, for instance Pearson)
def compute_local_ricci_from_matrix(M,max_dim,max_dist,nodes=None,precision=2):
   #nodes: a list with the labels of the nodes in the data provided. This will be useful to relabel the nodes that are called in the cliques found. If not provided, it will be understood as intergers (from 0 to n)
   #precision: The number of decimal places for setting the cutoff values
   
    if nodes==None:
        nodes=list(range(len(M)))
    
    
    #L_dict=list(D.keys())
    node_dict={i:nodes[i] for i in range(len(nodes))}
    
    assert len(nodes)==len(M), "Node labels size does not match with matrix size."
    
    #Setting the density of point from the precision provided
    len_den=len(precision*" ")
    delta_dist=10**(-precision)
    
    #Computing generator for the cliques
    WC=cliques_gudhi_matrix(np.round(M,precision),max_dist,max_dim)
   
    
    
    t0=time()

    #computing number of nodes
    nnodes=len(nodes)
    
    #computing headers for average and total FRC:
    heads1=["cutoff","density","avg_forman"]
    heads2=["cutoff","density","total_forman"]
    header1=heads1+nodes
    header2=heads2+nodes
   
    #Defining the set of neighbors per boundarty for each dimension dim. 
    N={dim:{} for dim in range(1,max_dim+1)}
   
    #dictionary for the number of cliques of dimension dim
    C={dim:0 for dim in range(1,max_dim+1)}
   
    #binomial of nnodes choose 2
    binominal=int(binom(nnodes,2))
    
    #auxiliar dictionary for computing global FRC for dimension dim:
    F={dim:0 for dim in range(1,max_dim+1)}
   
    #list of all the cutoffs for computing FRC:
    cutoffs=[round(i,len_den) for i in np.arange(0,max_dist+delta_dist,delta_dist)]
   
    #Local FRC (per node) for dimesions dim:
    nF={dim:{n:0 for n in nodes} for dim in range(1,max_dim+1)}

   
    t0=time()
    
    #All possible binom of nnodes choose dim (this will be used for computing the density):
    combs={d:binom(nnodes,d+1) for d in range(1,max_dim+1)}
   
    #Dictionary for computing Auxiliary distances (Aux_dist), auxiliary instant avg FRC (Aux_line1) and total FRC (Aux_line2) for each dim. The intial values for distances are set as zerom the auxiliary lines are empty strings to start the process:
    Aux_Dist={dim:0 for dim in range(1,max_dim+1)}
    Aux_Line1={dim:"" for dim in range(1,max_dim+1)}
    Aux_Line2={dim:"" for dim in range(1,max_dim+1)}
  
    #Auxiliary lines for Avg (L1) and total (L2) FRC for dimension d:
    L1={d:[] for d in range(1,max_dim+1)}
    L2={d:[] for d in range(1,max_dim+1)}
  

    #Startng the process of computing FRC's in function of the cutoff and the cell provided:
    for line in WC:
        #From gudhi output, each line is a tuple, where the first entry is the clique and the second is the weigh/distance/value associated to the clique provided:
    
        dist=line[1]
        
        #Relaibeling the clique from the list of node labels provided as input:
        clique=[node_dict[i] for i in line[0]]
        
        #computing the clique dimension
        d=len(clique)-1
        
        
        #Counting the clique for this dimension:
        C[d]+=1
            
        #Computing the boundary of the clique:   
        B=[frozenset(b) for b in combinations(clique,d)]

            
        #for each element in the boundary, do:
        for b0 in B:
            
            #update the neighborhood of the boundary
            try:
                N[d][b0]=N[d][b0].union(set(clique)-set(b0))
            except:
                N[d][b0]=set()
                N[d][b0]=N[d][b0].union(set(clique)-set(b0))
    
        #auxiliary set for the boundary of dimension d:
        N_aux=N[d]

        #computing the intersection of neighbors of the boundaries of the cell:
        H=set.intersection(*[N_aux[b0] for b0 in B])
        #computing the number of neighbors of the cell
        neigh=sum([len(N_aux[b0]) for b0 in B])

        #computing FRC of the cell:
        f=(d+2)*len(H)+2*(d+1)-neigh

        #update the FRC for each node of the cell:
        for node in clique:
              
            nF[d][node]+=f/(d+1)

        #update the contribution of the new cell in the simplicial complex in the neighborhood of the new clique:
        #FOr each element in the boundaty of the cell, do:
        for b0 in B:
            #define an auxiliary set for the neighbouring nodes of the boundary (except the ones in the cell itself)
            aux=N_aux[b0]-set(clique)
           
            #for each node in this auxiliary neighborhood, do:
            for a in aux:

                #compute delta as in the manuscript: if a is in H, the new neighbor is transverse (delta=d+1), otherwise is parallel (delta=-1):
                aux_delta=delta(a in H,d)
                val=(aux_delta)/(d+1)

                #update this information for each node in the boundary:
                for node in b0:
                    nF[d][node]+=val


                #update the information for the node itself:
                nF[d][a]+=val

                #update the contribution to the FRC:
                f+=aux_delta    

                   
    
        #Update the d-th total FRC
        F[d]+=f

           
        #create output list for avg (out1) and total (out2) FRC, globally and for each node in the graph:
        out1=[dist,C[d]/combs[d],F[d]/C[d]]+[nF[d][i]/C[d] for i in nF[d].keys()]
        out2=[dist,C[d]/combs[d],F[d]]+[nF[d][i] for i in nF[d].keys()]
            
            

            
        #Update the result in the dictionary only if the previous distance is different from the one in question: 
        if Aux_Dist[d]!=dist:
            aux1=Aux_Line1[d]
            aux2=Aux_Line2[d]
                
            if aux1!="":
                L1[d].append(aux1)
                L2[d].append(aux2)
                  

        Aux_Dist[d]=dist
        Aux_Line1[d]=out1
        Aux_Line2[d]=out2

    t1=time()
    #Finally, update the result for the last point of cutoff, for each dim:
    for dim in range(1,max_dim+1):
        L1[dim].append(Aux_Line1[dim])
        L2[dim].append(Aux_Line2[dim])
       
    #The FRC computation ends here. The following command are for pos processing, in order to fufill the missing informatrion between cutoffs and make the output a fency dictionary:
    Aux={d:{dict(zip(header1,i))["cutoff"]:dict(zip(header1,i)) for i in L1[d]} for d in range(1,max_dim+1)}
    null_dict=dict(zip(header1,[np.nan for i in range(len(header1))]))
    
    Output1={d:{cut:null_dict for cut in cutoffs} for d in range(1,max_dim+1)}
    
    for d in range(1,max_dim+1):
        min_cut=next((i for i in Aux[d].keys()))
        last_cutoff=min_cut
        for cut in cutoffs:
            if cut>=min_cut:
                if cut in Aux[d]:
            
                    Output1[d][cut]=Aux[d][cut].copy()
                    last_cutoff=cut
        
                else:
            
                    Output1[d][cut]=Aux[d][last_cutoff].copy()
                    Output1[d][cut]["cutoff"]=cut
                    Output1[d][cut]["density"]=np.nan
     
    Aux={d:{dict(zip(header2,i))["cutoff"]:dict(zip(header2,i)) for i in L2[d]} for d in range(1,max_dim+1)}
    null_dict=dict(zip(header2,[np.nan for i in range(len(header2))]))
    
    Output2={d:{cut:null_dict for cut in cutoffs} for d in range(1,max_dim+1)}
    
    for d in range(1,max_dim+1):
        min_cut=next((i for i in Aux[d].keys()))
        last_cutoff=min_cut
        for cut in cutoffs:
            if cut>=min_cut:
                if cut in Aux[d]:
            
                    Output2[d][cut]=Aux[d][cut].copy()
                    last_cutoff=cut
        
                else:
            
                    Output2[d][cut]=Aux[d][last_cutoff].copy()
                    Output2[d][cut]["cutoff"]=cut
                    Output2[d][cut]["density"]=np.nan
           
          


    t2=time()

    print("Total ricci computation time: "+str(t1-t0)+" seconds."+"\n")
    print("Total pos processing time: "+str(t2-t1)+" seconds."+"\n")
    print("Total computation time: "+str(t2-t0)+" seconds."+"\n")
  

    return Output1,Output2




def compute_local_ricci(D,max_dim,max_dist,nodes=None,precision=2):
  
    len_den=len(precision*" ")
    delta_dist=10**(-precision)
    
    
    
    WC=cliques_gudhi(D.values(),max_dist,max_dim)
   
    if nodes==None:
        nodes=list(D.keys())
        L_dict=list(D.values())
        node_dict={i:L_dict[i] for i in range(len(L_dict))}
        WC=cliques_gudhi(D.values(),max_dist,max_dim)
    else:
        node_dict={i:nodes[i] for i in range(len(nodes))}
        WC=cliques_gudhi(D.values(),max_dist,max_dim)

            
    assert len(D)==len(nodes), "Node labels size does not match with dictionary size."
  
    
    t0=time()

    
    nnodes=len(nodes)
   
    heads1=["cutoff","density","avg_forman"]
    heads2=["cutoff","density","total_forman"]
    header1=heads1+nodes
    header2=heads2+nodes
   
    N={dim:{} for dim in range(1,max_dim+1)}
   
    C={dim:0 for dim in range(1,max_dim+1)}
    #Setting FRC for all dims
    binominal=int(binom(nnodes,2))
    #R={dim:0 for dim in range(1,max_dim+1)}
    F={dim:0 for dim in range(1,max_dim+1)}
   
    cutoffs=[round(i,len_den) for i in np.arange(0,max_dist+delta_dist,delta_dist)]

    
    nF={dim:{n:0 for n in nodes} for dim in range(1,max_dim+1)}

   
    t0=time()
   
    combs={d:binom(nnodes,d+1) for d in range(1,max_dim+1)}
    
    Aux_Dist={dim:0 for dim in range(1,max_dim+1)}
    Aux_Line1={dim:"" for dim in range(1,max_dim+1)}
    Aux_Line2={dim:"" for dim in range(1,max_dim+1)}
   
    L1={d:[] for d in range(1,max_dim+1)}
    L2={d:[] for d in range(1,max_dim+1)}
   
    for line in WC:
    
     
        dist=round(line[1],len_den)
      
        clique=[node_dict[k] for k in line[0]]
       
    
        #clique dimension:
        d=len(clique)-1
        
        
        #update number of cliques:
        C[d]+=1
            

        B=[frozenset(b) for b in combinations(clique,d)]

           

        for b0 in B:
            

            try:
                N[d][b0]=N[d][b0].union(set(clique)-set(b0))
            except:
                N[d][b0]=set()
                N[d][b0]=N[d][b0].union(set(clique)-set(b0))
    
        
        N_aux=N[d]

        H=set.intersection(*[N_aux[b0] for b0 in B])
        neigh=sum([len(N_aux[b0]) for b0 in B])

        f=(d+2)*len(H)+2*(d+1)-neigh

        for node in clique:
              
            nF[d][node]+=f/(d+1)

        for b0 in B:
            aux=N_aux[b0]-set(clique)
           
        
            for a in aux:

                aux_delta=delta(a in H,d)
                val=(aux_delta)/(d+1)

                for node in b0:
                    nF[d][node]+=val


                nF[d][a]+=val

                f+=aux_delta    

                   
    

        F[d]+=f

           
        metric=[dist,C[d]/combs[d]]
        res=[F[d]]+[nF[d][i] for i in nF[d].keys()]
        out2=metric+res
        out1=metric+[i/C[d] for i in res]
       
      
        if Aux_Dist[d]!=dist:
            aux1=Aux_Line1[d]
            aux2=Aux_Line2[d]
                
            if aux1!="":
                L1[d].append(aux1)
                L2[d].append(aux2)
                    

        Aux_Dist[d]=dist
        Aux_Line1[d]=out1
        Aux_Line2[d]=out2


    for dim in range(1,max_dim+1):
        L1[dim].append(Aux_Line1[dim])
        L2[dim].append(Aux_Line2[dim])
        
    t1=time()
    
    Aux={d:{dict(zip(header1,i))["cutoff"]:dict(zip(header1,i)) for i in L1[d]} for d in range(1,max_dim+1)}
    null_dict=dict(zip(header1,[np.nan for i in range(len(header1))]))
    
    Output1={d:{cut:null_dict for cut in cutoffs} for d in range(1,max_dim+1)}
    
    for d in range(1,max_dim+1):
        min_cut=next((i for i in Aux[d].keys()))
        last_cutoff=min_cut
        for cut in cutoffs:
            if cut>=min_cut:
                if cut in Aux[d]:
            
                    Output1[d][cut]=Aux[d][cut].copy()
                    last_cutoff=cut
        
                else:
            
                    Output1[d][cut]=Aux[d][last_cutoff].copy()
                    Output1[d][cut]["cutoff"]=cut
                    Output1[d][cut]["density"]=np.nan
     
#
    Aux={d:{dict(zip(header2,i))["cutoff"]:dict(zip(header2,i)) for i in L2[d]} for d in range(1,max_dim+1)}
    null_dict=dict(zip(header2,[np.nan for i in range(len(header2))]))
    
    Output2={d:{cut:null_dict for cut in cutoffs} for d in range(1,max_dim+1)}
    
    for d in range(1,max_dim+1):
        min_cut=next((i for i in Aux[d].keys()))
        last_cutoff=min_cut
        for cut in cutoffs:
            if cut>=min_cut:
                if cut in Aux[d]:
            
                    Output2[d][cut]=Aux[d][cut].copy()
                    last_cutoff=cut
        
                else:
            
                    Output2[d][cut]=Aux[d][last_cutoff].copy()
                    Output2[d][cut]["cutoff"]=cut
                    Output2[d][cut]["density"]=np.nan
           
          


    t2=time()

    print("Total ricci computation time: "+str(t1-t0)+" seconds."+"\n")
    print("Total pos processing time: "+str(t2-t1)+" seconds."+"\n")
    print("Total computation time: "+str(t2-t0)+" seconds."+"\n")
   
    return Output1,Output2


def compute(obj,max_dim,max_dist,nodes=None,precision=2):
    
    if isinstance(obj,dict):
        if nodes==None:
            nodes=list(obj.keys())
        return compute_local_ricci(obj,max_dim,max_dist,nodes,precision)
    elif isinstance(obj,np.ndarray):
        if nodes==None:
            nodes=[i for i in range(len(obj))]
        return compute_local_ricci_from_matrix(obj,max_dim,max_dist,nodes,precision)



def matrix_from_dict(D):
    n=len(D)
    T={}
    L=list(D.values())
    for i in range(len(L)):
        T[i]=L[i]
    M=np.zeros((n,n))
    for pair in combinations([i for i in range(n)],2):
        i,j=pair
        M[i][j]=M[j][i]=dist(T[i],T[j])
    return M    
        
