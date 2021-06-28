# Repertoire-Comparison-Algorithm (RCA)

A macro-structure comparison algorithm for TCR repertoires developed by the ASU Biodesign Center for Biocomputing, Security and Society.
The formal description of the RCA algorithm is given in [this paper](TODO).

The main components of RCA are:

1. **Clustering**. Instead of comparing repertoires at the individual amino acid sequence level, RCA first identifies *clusters* of amino acid sequences.
A cluster is a connected component of the graph with amino acid sequences as nodes and edges between sequences that are within [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance) 1.
The relevant code for clustering can be found in `cluster.py`.

2. **Matching**. RCA then compares the clusters of two individuals, where a pair of clusters *match* if at least one sequence from the first cluster is within Hamming distance 1 of some sequence in the second cluster.
Any cluster that does not match is considered *missing*.
The relevant code for matching can be found in `match.py`.

**TODO**: Discuss sampling.

Finally, this repository also contains the analysis and plotting code used to create the paper's tables and figures in `analyze.py`.


## Getting Started

### 1. Prerequisites

You'll need a command line (Unix-based, Windows Command Prompt, or macOS Terminal) and any Python installation version 3.8 or newer.
You will also need (**TODO**: list packages, or make a `requirements.txt` file).

Clone this repository or download the latest [release](https://github.com/fesponda/Repertoire-Comparison-Algorithm/releases).


### 2. TCR Repertoire Data

The directory structure for RCA is:

```
data      # TCR sequencing data.
|--- dataset1
|--- |--- id1.tsv
|--- |--- id2.tsv
|--- |--- ...
|--- dataset2
|--- ...
clusters  # Clusters files from cluster.py.
|--- dataset1
|--- |--- id1.pkl
|--- |--- id2.pkl
|--- |--- ...
|--- dataset2
|--- ...
matches   # Matches files from match.py.
|--- dataset1.pkl
|--- dataset2.pkl
|--- ...
```

The `data/` directory contains the raw TCR high-throughput sequencing input data.
The contents of `clusters/` and `matches/` will be created after running the clustering and matching code, respectively.

The RCA paper uses the [esponda-2020]() and [emerson-2017-natgen](https://clients.adaptivebiotech.com/pub/emerson-2017-natgen) datasets.
You will need to obtain your input dataset of interest (e.g., `mydata`) and create and populate the corresponding `data/mydata` directory.
The code is compatible with any `.tsv` dataset that contains (1) the CDR3 amino acid sequences, e.g., `CASSSTGVEQFF` and (2) the number of times those sequences were read.
If your data is in a different format than the datasets used in the RCA paper (e.g., it has different column names), then you will need to modify the `read_data` function in `cluster.py` so the right columns are chosen.


### 3. Clustering

Clustering is performed by `cluster.py`.
At minimum, this code needs the dataset name corresponding to the folder containing the sequencing data (e.g., `mydata` for `data/mydata`) and the format of that data (see the `read_data` function).
You may optionally specify a minimum cluster size in terms of the fraction of sequences it must contain to be counted as a cluster (the default is `0.001`).

To cluster a single individual's sequencing data &mdash; e.g., those in `data/mydata/janedoe.tsv` &mdash; run:

```
# Option 1: Command line.
python cluster.py --id janedoe --dataset mydata --datafmt myfmt

# Option 2: Within the python shell.
>>> from cluster import *
>>> cluster_one('janedoe', 'mydata', 'myfmt', cluster_frac=0.001)
```

This will output a `ClusterSet` written as a binary to `clusters/mydata/janedoe.pkl`.

Alternatively, to cluster all individuals' data in your `mydata` dataset, run:

```
# Option 1: Command line.
python cluster.py --dataset mydata --datafmt myfmt

# Option 2: Within the python shell.
>>> from cluster import *
>>> cluster_all('mydata', 'myfmt', cluster_frac=0.001)
```


### 4. Matching

Matching is performed by `match.py`.
To find all matching clusters for each pair of individuals in a dataset `mydata`, run:

```
# Option 1: Command line.
python match.py --dataset mydata

# Option 2: Within the python shell.
>>> from match import *
>>> matches = match_all('mydata')
```

Both options write a dictionary detailing the matches to `matches/mydata.pkl`.


## Contributing

If you'd like to leave feedback, feel free to open a [new issue](https://github.com/fesponda/Repertoire-Comparison-Algorithm/issues/new/). If you'd like to contribute, please submit your code via a [pull request](https://github.com/fesponda/Repertoire-Comparison-Algorithm/pulls).
