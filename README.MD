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
- `Example.ipynb`

## Last version:
    0.4.0

## Python version: 
    3.8.5
## Package requirement:

- `numpy`
- `networkx`
- `gudhi`
- `cython`


## Installation:

The files `fastforman.pyx`, `compiler.py` and `setup.py` must be placed at the same directory.

 ### Local installation:
To this purpose, the user will need to compile the file `compiler.py` in cython by running the command:

```
python3 compiler.py build_ext --inplace
```

### Global installation:

Run `setup.py` by executing the command:

```
pip install .
```

After successful compilation, the generated `fastforman` file can be imported as a package, for instance:

```python
import fastforman as ff
```
***
## Functions:

- ``compute_FRC``
- ``compute_FRC_node``
- ``compute_average_FRC``
- ``compute_FRC_frequency``
- ``compute_FRC_node_frequency``
***

- ### ``compute_FRC``:

  ***Input:*** 
        
    - `D`, can be:
     - A `dictionary` in which the keys are integer numbers from 0 to len(D)-1 and the values are the N-dimensional points;

     - A symmetric `numpy.matrix` of float numbers;

     - A simple undirected `nx.Graph` that may include the feature `weights`on edges attribute,
     
     - A string of a graphml/xml file.

    - `cutoff`: Float number for the threshold distance between points.
    - `dim`: Integer, the maximum simplex dimension to compute FRC.
   
    
  ***Output:*** 
      A dictionary in which the keys are the dimension and the values are dictionary whose values keys are the simpleces and the values are the respective FRC.

***
    
- ### ``compute_FRC_node``:

  ***Input:*** 
        
    - `D`, can be:
     - A `dictionary` in which the keys are integer numbers from 0 to len(D)-1 and the values are the N-dimensional points;

     - A symmetric `numpy.matrix` of float numbers;

     - A simple undirected `nx.Graph` that may include the feature `weights`on edges attributes,
     
     - A string of a graphml/xml file.

    - `cutoff`: Float number for the threshold distance between points.
    - `dim`: Integer, the maximum simplex dimension to compute FRC.
   
    
   ***Output:*** 
   A dictionary in which the keys are the dimension and the values are dictionary whose values keys are the nodes and the values are the respective FRC. When the node does not have curvature it will be set as numpy.nan.
    
***

- ### ``compute_average_FRC``:

  ***Input:*** 
        
    - `D`, can be:
     - A `dictionary` in which the keys are integer numbers from 0 to len(D)-1 and the values are the N-dimensional points;

     - A symmetric `numpy.matrix` of float numbers;

     - A simple undirected `nx.Graph` that may include the feature `weights`on edges attributes,
     
     - A string of a graphml/xml file.

    - `cutoff`: Float number for the threshold distance between points.
    - `dim`: Integer, the maximum simplex dimension to compute FRC.
   
    
  ***Output:*** 
     
   A dictionary in which the keys are the dimension and the values are the average FRC for fixed dimension.
    
***

- ### ``compute_FRC_frequency``:

     ***Input:*** 
         
    - `D`, can be:
     - A `dictionary` in which the keys are integer numbers from 0 to len(D)-1 and the values are the N-dimensional points;

     - A symmetric `numpy.matrix` of float numbers;

     - A simple undirected `nx.Graph` that may include the feature `weights`on edges attributes,
     
     - A string of a graphml/xml file.

    - `cutoff`: Float number for the threshold distance between points.
    - `dim`: Integer, the maximum simplex dimension to compute FRC.
    
  ***Output:*** A dictionary in which the keys are the dimension and the values is a dictionary of the frequency of FRC values. 

***

- ### ``compute_FRC_node_frequency``:

     ***Input:*** 
         
    - `D`, can be:
     - A `dictionary` in which the keys are integer numbers from 0 to len(D)-1 and the values are the N-dimensional points;

     - A symmetric `numpy.matrix` of float numbers;

     - A simple undirected `nx.Graph` that may include the feature `weights`on edges attributes,
     
     - A string of a graphml/xml file.

    - `cutoff`: Float number for the threshold distance between points.
    - `dim`: Integer, the maximum simplex dimension to compute FRC.
    
  ***Output:*** A dictionary in which the keys are the dimension and the values are dictionaries of the frequency of FRC values per node. 
 
***
Check out our examples.

## Contact and Support:

danillo.dbs16@gmail.com, dbarros@bcamath.org

## References: 

If you appreciate our work, please, cite us:

[1] Efficient set-theoretic algorithms for computing high-order Forman-Ricci curvature on abstract simplicial complexes; Barros de Souza, Danillo; da Cunha, Jontatas, A.N. Santos, Fernando; Jost, Juergen & Rodigues, Serafim; Link: https://arxiv.org/html/2308.11763v4

[2] Forman-Ricci curvature Benchmark report, Kaggle repository, Link: https://www.kaggle.com/datasets/danillosouza2020/forman-ricci-curvature-benchmark-report