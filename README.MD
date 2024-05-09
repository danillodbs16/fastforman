# FastForman - An efficient Forman-Ricci Curvature computation for higher-order faces in Simplicial Complexes

***
   **Author:** Danillo Barros de Souza
   -- 
   **ORCID: [0000-0002-7792-8862](https://orcid.org/0000-0002-7762-8862)**

***
This work is based on our set theoretical approach for computing Forman-Ricci curvature [[1]](https://arxiv.org/abs/2308.11763). We implemented our code in Python language and coinded **fastforman**. We profiled the time processing and memory usage and benchmarked with other alogorithms in the literature, namelly, **GeneralisedFormanRicci** and **HodgeLaplacians**. The benchmark codes and reports can be found at [[2]](https://www.kaggle.com/datasets/danillosouza2020/forman-ricci-curvature-benchmark-report). 

 Here, we provide a set-theoretical Python coding for efficiently computing Forman-Ricci Curvature (FRC) for simplicial complexes.

## Content:

- `fastforman.pyx`
- `compiler.py`
- `setup.py`

## Python version: 
    3.8.5
## Package requirement:

- `numpy`
- `gudhi`
- `cython`



## Installation:
   ### Option 1: local installation:

   The will need to compile the file `compiler.py` in cython by running the following command:
   
   ```python3 compiler.py build_ext --inplace```
    
    
  
   ### Option 2: global installation:
   Execute `setup.py` file by running in the command line:
   
   ```pip install . ```
   

After successful installation, the `fastforman` can be imported as a package, for instance:

```python
import fastforman as ff
```
## Functions:
- ``compute_average_FRC``:

    ***Input:*** The input parameters can be in two ways:
        
    - `Tuple, max_dim`, where  `Tuple=(D,Thres)` is a tuple containing a dictionary of points (the keys are the points enumeration and the values are the n-dimensional points coordinates) and the distance threshold,`Thres`,  and `(D,max_dim)` is the maximum face dimension, `max_dim`;
    
    - `G, max_dim`, where `G` is a simple undirected network graph, and is `max_dim` the maximum dimension.
    
  ***Output:*** A dictionary in which the keys are the dimension and the values are the average FRC for fixed        dimension.
    
    
- ``compute_FRC_frequency``:

    ***Input:*** Same as in `compute_average_FRC`.
    
     ***Output:*** A dictionary in which the keys are the dimension and the values is a dictionary of the frequency of FRC values.
    
Below, we provide some examples for clarification.

## Contact and Support:

danillo.dbs16@gmail.com, dbarros@bcamath.org

## References: 

If you appreciate our work, please, cite us:
[1] Efficient set-theoretic algorithms for computing high-order Forman-Ricci curvature on abstract simplicial complexes; Barros de Souza, Danillo; da Cunha, Jontatas, A.N. Santos, Fernando; Jost, Juergen & Rodigues, Serafim; Link: https://arxiv.org/html/2308.11763v4

[2] Forman-Ricci curvature Benchmark report, Kaggle repository, Link: https://www.kaggle.com/datasets/danillosouza2020/forman-ricci-curvature-benchmark-report