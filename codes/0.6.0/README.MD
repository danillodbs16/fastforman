# FastForman - An efficient Forman-Ricci Curvature computation for higher-order faces in Simplicial Complexes

---

**Author:** Danillo Barros de Souza

**ORCID:** [0000-0002-7792-8862](https://orcid.org/0000-0002-7792-8862)

---

This work is based on our set-theoretical approach for computing Forman-Ricci curvature [\[1\]](https://royalsocietypublishing.org/doi/10.1098/rspa.2024.0364). We implemented our code in Python and coined **fastforman**. We profiled time and memory usage and benchmarked against other algorithms in the literature, namely, **GeneralisedFormanRicci** and **HodgeLaplacians**. The benchmark codes and reports can be found at [\[2\]](https://www.kaggle.com/datasets/danillosouza2020/forman-ricci-curvature-benchmark-report).

Here, we provide a set-theoretical Python implementation for efficiently computing Forman-Ricci Curvature (FRC) for simplicial complexes. We use our own methods for computing the cliques in the simplicial complexes. For more details, see [\[3\]](https://arxiv.org/abs/2502.14593). For alternative topological computations, see [\[4\]](https://github.com/danillodbs16/fastknill).

## Current version:

```
0.6.0
```

## New Features:

* Added FRC computations of filtrations on Vietoris-Rips Complexes;
* No need for importing external packages for computing cliques (e.g., `networkx` and `gudhi`);
* Faster computations and reduced memory usage;
* Input parameters `n_neighbors` and `mapped_nodes` added;

**Caution: The root package structure has changed. See previous versions for comparison.**

## Content:

* `setup.py`
* `fastforman/`

  * `GeometricSimplicialComplex.pyx`
  * `VietorisRipsComplex.pyx`
* `Example.ipynb`

## Python version:

```
3.8.5
```

## Package requirements:

* `numpy`
* `networkx`
* `cython`
* `scikit-learn`

## Installation:

The folder `fastforman` and the file `setup.py` must be placed in the same directory.

### Local installation:

To compile `setup.py` with Cython, run:

```
python3 setup.py build_ext --inplace
```

### Global installation:

Run:

```
pip install .
```

After successful compilation, the package can be imported as:

```python
import fastforman as ff
```

## Functions:

### 1. Computations on Geometric Simplicial Complexes

Use case: Study of geometric information from weighted and unweighted abstract simplicial complexes (in a static perspective).

---

* `ff.GeometricSimplicialComplex.compute_FRC`
  Computes the Forman-Ricci Curvature (FRC) up to dimension `dim`.

  **Input:**

  * `D`: One of the following:

    * A `dict` mapping integers `0` to `len(D)-1` to N-dimensional points;
    * A symmetric `numpy.matrix` of floats;
    * An undirected `nx.Graph` (optionally weighted);
    * A file path to a `graphml` file.
  * `cutoff`: Float threshold distance.
  * `n_neighbors`: Integer for nearest neighbors (optional).
  * `dim`: Maximum simplex dimension.
  * `mapped_nodes`: Boolean, preserve original node labels.

  **Output:**
  A `dict` where keys are dimensions and values are `dict`s mapping simplices to FRC values.

---

* `ff.GeometricSimplicialComplex.compute_FRC_node`
  Computes the FRC for nodes up to dimension `dim`.

  **Input:** Same as above.

  **Output:**
  A `dict` with dimensions as keys, and values being `dict`s mapping nodes to FRC values. Unavailable values are `np.nan`.

---

* `ff.GeometricSimplicialComplex.compute_average_FRC`
  Computes the average FRC up to a given dimension.

  **Input:** Same as above.

  **Output:**
  A `dict` mapping dimensions to average FRC values.

---

* `ff.GeometricSimplicialComplex.compute_FRC_frequency`
  Computes frequency distribution of FRC values up to a given dimension.

  **Input:** Same as above.

  **Output:**
  A `dict` where keys are dimensions and values are dictionaries of frequency distributions.

---

* `compute_FRC_node_frequency`
  Computes frequency of FRC values for nodes.

  **Input:**

  * `D`: As described previously.
  * `cutoff`, `n_neighbors`, `mapped_nodes`, `dim`: Same definitions.

  **Output:**
  A `dict` mapping dimensions to frequency dictionaries for node-level FRC.

---

### 2. Computations on Vietoris-Rips Complexes

Use case: Study of geometric information across filtrations.

* `ff.VietorisRipsComplex.compute_local_FRC`

  **Input:**

  * `M`: `np.ndarray`, point cloud;
  * `max_dim`: Max simplex dimension;
  * `max_dist`: Cutoff distance;
  * `k_nearest`: Max number of neighbors;
  * `precision`: Decimal rounding for distances;
  * `metric`: Distance metric (e.g., "euclidean");
  * `verbose`: Enable debugging output.

  **Output:**

  * `Output1`: `dict`, total d-FRC per node;
  * `Output2`: `dict`, average d-FRC per node.

---

## Contact and Support:

[danillo.dbs16@gmail.com](mailto:danillo.dbs16@gmail.com), [dbarros@bcamath.org](mailto:dbarros@bcamath.org)

---

## References:

\[1] Barros de Souza, Danillo; da Cunha, Jontatas; Santos, Fernando A.N.; Jost, Juergen; Rodrigues, Serafim (2024).
**Efficient set-theoretic algorithms for computing high-order Forman-Ricci curvature on abstract simplicial complexes**.
*Proc. R. Soc. A.* [https://doi.org/10.1098/rspa.2024.0364](https://royalsocietypublishing.org/doi/10.1098/rspa.2024.0364)

\[2] Barros de Souza, Danillo (2024).
**Forman-Ricci curvature benchmark report**. Kaggle Repository.
[https://www.kaggle.com/datasets/danillosouza2020/forman-ricci-curvature-benchmark-report](https://www.kaggle.com/datasets/danillosouza2020/forman-ricci-curvature-benchmark-report)

\[3] Barros de Souza, Danillo (2024).
**Alternative set-theoretical algorithms for efficient computations of cliques in Vietoris-Rips complexes**.
[https://arxiv.org/abs/2502.14593](https://arxiv.org/abs/2502.14593)

\[4] Barros de Souza, Danillo.
**FastKnill - Fast computation of Euler characteristics and Knill curvature**.
[https://github.com/danillodbs16/fastknill](https://github.com/danillodbs16/fastknill)
