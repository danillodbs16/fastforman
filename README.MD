# FastForman - An efficient Forman-Ricci Curvature computation for higher-order faces in Simplicial Complexes

 Here, we provide a set-theoretical Python coding for effifiently compute Forman-Ricci Curvature (FRC) for simplicial complexes.

This code development was based on our work on [1]. A benchmark comparison with other Python algorithms in the lilterature for computing FRC can be found at [2].

## Content:

- `fastforman.pyx`
-`compyler.py`
## Python version: 
    3.8.5
## Package requirement:

- `numpy==1.23.0`
-`networkx==2.6.3`
- `gudhi==3.5.0`
- `cython==0.29.21`
- `distutils==3.8.5`


## Instalation:
To this purpuse, the user will need to compile the file `compiler.py` in cython by running the command:

```python
python3 compiler.py build_ext --inplace
```

After successfull compilation, the generated `fastforman` file can be inported as a package, for instance:

```python
import fastforman as ff
```
## Functions:
- ``compute_average_FRC``:

    ***Input:*** Either a dictionary of points (the keys are the points enumeration and the values are the n-dimensional points coordinates) or a networkx graph `G`, and the maximum dimension, `(D,max_dim)`.

    ***Output:*** A dictionary in which the keys are the dimension and the values are the average FRC for fixed dimension.
    
    
- ``compute_FRC_frequency``:

    ***Input:*** Either a dictionary of points (the keys are the points enumeration and the values are the n-dimensional points coordinates) or a networkx graph `G`, and the maximum dimension, `(D,max_dim)`.

    ***Output:*** A dictionary in which the keys are the dimension and the values is a dictionary of the frequency of FRC values.
    
Bellow, we provide some examples.

## Contact and Support:

danillo.dbs16@gmail.com, dbarros@bcamath.org

## References: 

[1] Efficient set-theoretic algorithms for computing high-order Forman-Ricci curvature on abstract simplicial complexes; Barros de Souza, Danillo; da Cunha, Jontatas, A.N. Santos, Fernando; Jost, Juergen & Rodigues, Serafim; Link: https://arxiv.org/html/2308.11763v4

[2] Forman-Ricci curvature Benchmark report, Kaggle repository, Link: https://www.kaggle.com/datasets/danillosouza2020/forman-ricci-curvature-benchmark-report
