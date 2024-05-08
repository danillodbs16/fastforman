#importing package needed:
from itertools import combinations

#defining Forman-Ricci curvature for each cell:
cpdef int ricci_cell(list c,dict Neigh):

    cdef int d=len(c)
    cdef set c0=set(c)
    
    cdef int H=len(set.intersection(*[Neigh[j] for j in c]))
    cdef tuple i
      
    cdef list l=[i for i in combinations(c,d-1)]

    cdef tuple k
    cdef int j

    N=sum([len(set.intersection(*[Neigh[j] for j in k])) for k in l])

    cdef int ricci=(d+1)*H+2*d-N
           
    return ricci
#defining Forman-Ricci curvature for the entire simplicial complex:
cpdef dict fastforman(L, dict Neigh):
    cdef list c
    cdef int ricci
    cdef tuple c0

    forman=dict()
    for c in L:
        try:
            ricci=ricci_cell(c,Neigh)
            c0=tuple(c)
            forman[c0]=ricci
        except:
            pass
    return forman        