# RCA 
 Clustering matching and comparison algorithms for TCR repertoires

## Terminology

The following is terminology that will be used throughout the documentation and code.

| Term | Description | Domain |
| --- | --- | --- |
| (DNA) nucleobase | the four bases of DNA | `{'G', 'A', 'C', 'T'}` |
| codon | a string of three nucleobases | `{'G', 'A', 'C', 'T'}^3` |
| amino acid | one of the twenty amino acids | `{'A', 'R', ..., 'V'}` |
| sequence | a string of one of more amino acids | `{'A', 'R', ..., 'V'}^+` |
| repertoire | a pandas DataFrame of indexed amino acid sequences | N/A |

- TODO: Cluster: a list of indices corresponding to sequences that are connected via a BFS search where sequences are neighbors if and only if their (Hamming) distance is 1.

IDEA: Create a new class for clusters that contain as members:
- the data file name that it corresponds to
- the sequence indices that the cluster comprises
- statistics, like the number of reads for this cluster's sequences
