# Project:   PublicClusters
# Filename:  cluster.py
# Authors:   Fernando Esponda (fernando.esponda@itam.mx) and Joshua J. Daymude
#            (jdaymude@asu.edu).

"""
cluster: Cluster repertoires of amino acid sequences.
"""

import argparse
from os import listdir
from os.path import join
import pandas as pd
import pickle
from tqdm import tqdm


def amino_acids():
    """
    :returns: a list of the 20 amino acids and the wildcard '*'
    """
    return ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', \
            'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']


def read_data(fname, datafmt):
    """
    Reads data from the given file as a DataFrame, removing all columns except
    the amino acid sequences and counts, grouping and sorting rows by count.

    :param fname: a string file name to read from (should be .tsv)
    :param datafmt: a string data format ('esponda2020' or 'emerson2017')
    :returns: a DataFrame containing amino acid sequences and their counts
    """
    # Read the data as a DataFrame, renaming columns according to the dataset.
    df = pd.read_csv(fname, sep='\t')
    if datafmt == 'esponda2020':
        col_map = {'aminoAcid': 'amino_seq', 'count (templates/reads)': 'count'}
    elif datafmt == 'emerson2017':
        col_map = {'amino_acid': 'amino_seq', 'reads': 'count'}
    else:
        raise ValueError('unrecognized data format \"' + datafmt + '\"')
    df = df.rename(columns=col_map)

    # Group by amino acid sequence and sort from most to least reads.
    df = df.groupby('amino_seq', as_index=False)['count'].sum()
    df = df.sort_values(by=['count'], ascending=False)

    return df


def close_to(amino_seq, repertoire, dist_metric='Hamming'):
    """
    Computes the set of amino sequences in the given repertoire set that are
    within distance 1 of the given amino sequence.

    :param amino_seq: a string representing an amino acid sequence
    :param repertoire: a dict mapping amino acid sequences to their indices
    :param dist_metric: a string distance metric ('Hamming' or 'Levenshtein')
    :returns: a set of indices of amino acid sequences within distance 1
    """
    nbr_seq_idxs = set()
    if dist_metric == 'Hamming':
        for i in range(len(amino_seq)):
            for a in amino_acids():
                # Substitute the i-th amino acid with amino acid a.
                one_sub = amino_seq[0:i] + a + amino_seq[i+1:]
                if one_sub in repertoire.keys():
                    nbr_seq_idxs.add(repertoire[one_sub])
                    repertoire.pop(one_sub)
    elif dist_metric == 'Levenshtein':
        for i in range(len(amino_seq)):
            for a in amino_acids():
                # Substitute the i-th amino acid with amino acid a.
                one_sub = amino_seq[0:i] + a + amino_seq[i+1:]
                if one_sub in repertoire.keys():
                    nbr_seq_idxs.add(repertoire[one_sub])
                    repertoire.pop(one_sub)
                # Insert amino acid a just before the i-th amino acid.
                one_insert = amino_seq[0:i] + a + amino_seq[i:]
                if one_insert in repertoire.keys():
                    nbr_seq_idxs.add(repertoire[one_insert])
                    repertoire.pop(one_insert)
            # Delete the i-th amino acid.
            one_delete = amino_seq[0:i] + amino_seq[i+1:]
            if one_delete in repertoire.keys():
                nbr_seq_idxs.add(repertoire[one_delete])
                repertoire.pop(one_delete)
        for a in amino_acids():
            # Insert amino acid a at the end of the sequence.
            end_insert = amino_seq + a
            if end_insert in repertoire.keys():
                nbr_seq_idxs.add(repertoire[end_insert])
                repertoire.pop(end_insert)

    return nbr_seq_idxs


def dbscan(df_amino, min_cluster_size, verbose):
    """
    Implements the DBSCAN algorithm to find at most the given maximum number of
    clusters each containing at least the given minimum number of amino acid
    sequences from the amino acid data.

    :param df_amino: a DataFrame containing a column of amino acid sequences
    :param min_cluster_size: an integer minimum number of sequences per cluster
    :param verbose: True iff progress info should be printed
    :returns: a list of clusters, each of which are lists of sequence indices
    """
    repertoire = dict(zip(df_amino['amino_seq'], df_amino.index))
    repertoire_idxs = set(df_amino.index.values)
    clusters = []

    # While there are enough unclustered sequences to possibly form another
    # cluster, find one cluster from among the unclustered sequences.
    while len(repertoire_idxs) > min_cluster_size:
        explore_idxs = set([repertoire_idxs.pop()])
        cluster = []
        while len(explore_idxs) > 0:
            seq_idx = explore_idxs.pop()
            cluster.append(seq_idx)
            nbr_seq_idxs = close_to(df_amino['amino_seq'][seq_idx], repertoire)
            repertoire_idxs -= nbr_seq_idxs
            explore_idxs |= (nbr_seq_idxs - set([seq_idx]))

        if len(cluster) >= min_cluster_size:
            clusters.append(cluster.copy())
            if verbose:
                print('Found cluster #{} containing {} sequences'\
                      .format(len(clusters), len(cluster)))

    return clusters


class ClusterSet:
    """
    A collection of clusters from the same individual's sequencing data.
    """

    def __init__(self, id, dataset, datafmt, cluster_frac):
        """
        Initializes a new ClusterSet, reading an individual's sequencing data
        from data/<dataset>/<id>.tsv.

        :param id: a string identifier for this ClusterSet
        :param dataset: a string dataset name
        :param datafmt: a string data format ('esponda2020' or 'emerson2017')
        :param cluster_frac: a float representing the fraction of total
        sequences that must appear in a cluster for it to be recorded
        """
        self.id = id
        self.dataset = dataset
        self.seq_df = read_data(join('data', dataset, id + '.tsv'), datafmt)
        self.min_cluster_size = max(int(cluster_frac * len(self.seq_df)), 5)
        self.clusters = []           # List of lists of sequence indices.
        self.clustered_seqs = set()  # Union of sequence indices in clusters.


    def find_clusters(self, verbose=True):
        """
        Finds all clusters in the ClusterSet's sequence data that contain at
        least the given fraction of total sequences.

        :param verbose: True iff progress info should be printed
        """
        self.clusters = dbscan(self.seq_df, self.min_cluster_size, verbose)
        for cluster in self.clusters:
            self.clustered_seqs |= set(cluster)


    def info(self):
        """
        :returns: a string description of this ClusterSet
        """
        return '\"' + self.dataset + '/' + self.id + '\" has ' + \
               '{} clusters clustering {} of {} total sequences'\
               .format(len(self.clusters), len(self.clustered_seqs), \
                       len(self.seq_df))


    def save(self):
        """
        Saves this ClusterSet object to clusters/<dataset>/<id>.pkl.
        """
        with open(join('clusters', self.dataset, self.id + '.pkl'), 'wb') as f:
            pickle.dump(self, f)


def load_clusterset(fname):
    """
    Loads a ClusterSet from the given file.

    :param fname: a string file name to load the ClusterSet from
    :returns: the ClusterSet loaded from file
    """
    with open(fname, 'rb') as f:
        return pickle.load(f)


def cluster_one(id, dataset, datafmt, cluster_frac, verbose=True):
    """
    Perform clustering on the given individual, writing results to file.

    :param id: a string identifier for the individual
    :param dataset: a string dataset name
    :param datafmt: a string data format ('esponda2020' or 'emerson2017')
    :param cluster_frac: a float threshold for minimum cluster size in terms of
    the total number of sequences
    :param verbose: True iff progress info should be printed
    """
    cs = ClusterSet(id, dataset, datafmt, cluster_frac)
    cs.find_clusters(verbose)
    tqdm.write(cs.info())
    cs.save()


def cluster_all(dataset, datafmt, cluster_frac):
    """
    Iteratively find clusters for all individuals in the given dataset.

    :param dataset: a string dataset name
    :param datafmt: a string data format ('esponda2020' or 'emerson2017')
    :param cluster_frac: a float threshold for minimum cluster size in terms of
    the total number of sequences
    """
    files = listdir(join('data', dataset))
    verbose = len(files) == 1
    for f in tqdm(files, desc='Clustering individual'):
        cluster_one(f.split('.')[0], dataset, datafmt, cluster_frac, verbose)


if __name__ == "__main__":
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-I', '--id', type=str, default='all', \
                        help='ID of individual to cluster or \'all\'')
    parser.add_argument('-D', '--dataset', type=str, required=True, \
                        help='Dataset to perform clustering on')
    parser.add_argument('-F', '--datafmt', type=str, required=True, \
                        help='Format of the sequencing data')
    parser.add_argument('-C', '--cluster_frac', type=float, default=0.001, \
                        help='Minimum cluster size as fraction of #sequences')
    args = parser.parse_args()

    if args.id == 'all':  # Cluster all individuals in the given dataset.
        cluster_all(args.dataset, args.datafmt, args.cluster_frac)
    else:  # Cluster a single individual.
        cluster_one(args.id, args.dataset, args.datafmt, args.cluster_frac)
