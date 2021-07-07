# Project:  PublicClusters
# Filename: cluster.py
# Authors:  Fabio Calo (fcalodiz@itam.mx), Joshua J. Daymude (jdaymude@asu.edu),
#           and Fernando Esponda (fernando.esponda@itam.mx).

"""
cluster: Cluster repertoires of amino acid sequences.
"""

from helper import *

import argparse
import multiprocessing as mp
import numpy as np
from tqdm import tqdm, trange


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


def dbscan(df_amino, min_cluster_size):
    """
    Implements the DBSCAN algorithm to find clusters of amino acid sequences
    containing at least the given minimum number of sequences.

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

    return clusters


class ClusterSet:
    """
    A collection of clusters from the same individual's sequencing data.
    """

    def __init__(self, id, dataset, datafmt, cluster_frac, sample_num=None, \
                 sample_size=None, rng=None):
        """
        Initializes a new ClusterSet, reading an individual's sequencing data
        from data/<dataset>/<id>.tsv.

        :param id: a string identifier for this ClusterSet
        :param dataset: a string dataset name
        :param datafmt: a string data format ('esponda2020' or 'emerson2017')
        :param cluster_frac: a float representing the fraction of total
        sequences that must appear in a cluster for it to be recorded
        :param sample_num: an int sample identifier, or None if not a sample
        :param sample_size: an int number of sequences for the sample
        :param rng: a numpy.random Generator for random numbers
        """
        self.id = id
        self.dataset = dataset
        self.seq_df = read_data(osp.join('data', dataset, id + '.tsv'), datafmt)
        self.min_cluster_size = max(int(cluster_frac * len(self.seq_df)), 5)
        self.sample_num = sample_num
        self.clusters = []           # List of lists of sequence indices.
        self.clustered_seqs = set()  # Union of sequence indices in clusters.

        # Perform sampling, if applicable.
        if self.sample_num is not None:
            idxs = list(self.seq_df.index)
            probs = list(self.seq_df['count'] / self.seq_df['count'].sum())
            idxs = rng.choice(idxs, size=sample_size, replace=False, p=probs)
            self.seq_df = self.seq_df.loc[idxs]


    def find_clusters(self):
        """
        Finds all clusters in the ClusterSet's sequence data that contain at
        least the given fraction of total sequences.
        """
        self.clusters = dbscan(self.seq_df, self.min_cluster_size)
        for cluster in self.clusters:
            self.clustered_seqs |= set(cluster)


    def info(self):
        """
        :returns: a string description of this ClusterSet
        """
        info = '\"' + self.dataset + '/'
        info += 'full/' if self.sample_num is None else \
                'sample' + str(self.sample_num) + '/'
        info += self.id + '\" has {} clusters clustering {} of {} sequences'\
                .format(len(self.clusters), len(self.clustered_seqs), \
                        len(self.seq_df))
        return info


    def save(self):
        """
        Saves this ClusterSet object to clusters/<dataset>/<sample>/<id>.pkl.
        """
        sample_str = 'full' if self.sample_num is None else \
                     'sample' + str(self.sample_num)
        fname = osp.join('clusters', self.dataset, sample_str, self.id + '.pkl')
        dump_obj(fname, self)


def cluster_chunk(ids_chunk, seeds_chunk, dataset, datafmt, cluster_frac, \
                  num_samples, sample_size):
    """
    Perform clustering on a chunk of individuals, writing results to file.

    :param ids_chunk: a list of string identifiers for individuals
    :param id_seeds_chunk: a list of int seeds for random number generation
    :param dataset: a string dataset name
    :param datafmt: a string data format ('esponda2020' or 'emerson2017')
    :param cluster_frac: a float threshold for minimum cluster size in terms of
    the total number of sequences
    :param num_samples: an int number of samples to generate
    :param sample_size: an int number of sequences per sample
    """
    # First, cluster the individual's entire repertoire; then generate and
    # cluster repertoire samples.
    pbar = tqdm(list(zip(ids_chunk, seeds_chunk)))
    for id, seed in pbar:
        pbar.set_description('PID{} clustering {}'.format(os.getpid(), id))
        rng = np.random.default_rng(seed)
        for i in range(num_samples + 1):
            cs = ClusterSet(id, dataset, datafmt, cluster_frac) if i == 0 else \
                 ClusterSet(id, dataset, datafmt, cluster_frac, i, sample_size, rng)
            cs.find_clusters()
            tqdm.write(cs.info())
            cs.save()


def cluster_all(dataset, datafmt, cluster_frac, num_samples, sample_size, seed,\
                num_procs=1):
    """
    Iteratively find clusters for all individuals in the given dataset.

    :param dataset: a string dataset name
    :param datafmt: a string data format ('esponda2020' or 'emerson2017')
    :param cluster_frac: a float threshold for minimum cluster size in terms of
    the total number of sequences
    :param num_samples: an int number of samples to generate
    :param sample_size: an int number of sequences per sample
    :param seed: an int seed for random number generation
    :param num_procs: an int number of processors to parallelize over
    """
    # Get list of all individuals and generate a random seed for each.
    ids = [f.split('.')[0] for f in os.listdir(osp.join('data', dataset))]
    seeds = np.random.default_rng(seed).integers(0, 2**32, size=len(ids))

    # Partition these lists over the number of processors.
    ids_chunks = np.array_split(ids, num_procs)
    seeds_chunks = np.array_split(seeds, num_procs)

    # Start all processes, clustering each individual and their samples.
    procs = []
    for ids_chunk, seeds_chunk in zip(ids_chunks, seeds_chunks):
        proc = mp.Process(target=cluster_chunk, args=(ids_chunk, seeds_chunk, \
                          dataset, datafmt, cluster_frac, num_samples, \
                          sample_size,))
        proc.start()
        procs.append(proc)

    # Wait until all processes have finished.
    for proc in procs:
        proc.join()


if __name__ == "__main__":
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-D', '--dataset', type=str, required=True, \
                        help='Dataset to perform clustering on')
    parser.add_argument('-F', '--datafmt', type=str, required=True, \
                        help='Format of the sequencing data')
    parser.add_argument('-C', '--cluster_frac', type=float, default=0.001, \
                        help='Minimum cluster size as fraction of #sequences')
    parser.add_argument('-N', '--num_samples', type=int, default=0, \
                        help='Number of samples to generate')
    parser.add_argument('-S', '--sample_size', type=int, default=41000, \
                        help='Number of sequences per sample')
    parser.add_argument('-R', '--rand_seed', type=int, default=None, \
                        help='Seed for random number generation')
    parser.add_argument('-P', '--num_procs', type=int, default=1, \
                        help='Number of processors to parallelize over')
    args = parser.parse_args()

    cluster_all(args.dataset, args.datafmt, args.cluster_frac, \
                args.num_samples, args.sample_size, args.rand_seed, \
                args.num_procs)
