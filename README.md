# RCA
 Clustering matching and comparison algorithms for TCR repertoires


## Getting Started

TODO: describe. File structure is:

```
data  # Sequencing data.
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
matches  # Matches files from match.py.
|--- dataset1.pkl
|--- dataset2.pkl
|--- ...
```

TODO: Explain how to configure new datasets. Need to add folders as above, choice to argparse, and an elif case in cluster::read_data


## Terminology

The following is terminology that will be used throughout the documentation and code.

| Term | Description | Domain |
| --- | --- | --- |
| amino acid | one of the twenty amino acids | `{'A', 'R', ..., 'V'}` |
| sequence | a string of one of more amino acids | `{'A', 'R', ..., 'V'}^+` |
| repertoire | a pandas DataFrame of indexed amino acid sequences | N/A |

- TODO: Cluster: a list of indices corresponding to sequences that are connected via a BFS search where sequences are neighbors if and only if their (Hamming) distance is 1.

IDEA: Create a new class for clusters that contain as members:
- the data file name that it corresponds to
- the sequence indices that the cluster comprises
- statistics, like the number of reads for this cluster's sequences
