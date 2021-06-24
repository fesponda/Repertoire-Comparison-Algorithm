# Project:   PublicClusters
# Filename:  cluster.py
# Authors:   Fernando Esponda (fernando.esponda@itam.mx) and Joshua J. Daymude
#            (jdaymude@asu.edu).

import argparse
import datetime
from os import listdir
from os.path import isfile, join
import pandas as pd
import pickle
from time import time


def amino_acids():
    """
    :returns: a list of the 20 amino acids and the wildcard '*'
    """
    return ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', \
            'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']


def codon2amino():
    """
    :returns: a dict mapping codons to amino acids
    """
    return {'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', \
            'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R', \
            'AAT':'N', 'AAC':'N', \
            'GAT':'D', 'GAC':'D', \
            'TGT':'C', 'TGC':'C', \
            'CAA':'Q', 'CAG':'Q', \
            'GAA':'E', 'GAG':'E', \
            'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', \
            'CAT':'H', 'CAC':'H', \
            'ATT':'I', 'ATC':'I', 'ATA':'I', \
            'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'TTA':'L', 'TTG':'L', \
            'AAA':'K', 'AAG':'K', \
            'ATG':'M', \
            'TTT':'F', 'TTC':'F', \
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', \
            'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S', \
            'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', \
            'TGG':'W', \
            'TAT':'Y', 'TAC':'Y', \
            'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', \
            'TAA':'Stop', 'TAG':'Stop', 'TGA':'Stop'}


def amino_acid_dists():
    """
    For each pair of amino acids, compute the minimum Hamming distance between
    pairs of their codons.

    :returns: a dict mapping amino acid pairs (e.g., 'IL', 'HH', or '*T') to the
    minimum Hamming distance between any pair of their codons.
    """
    # Invert codon2amino to obtain a mapping of amino acids to lists of codons.
    amino2codon = {}
    for amino in codon2amino().values():
        amino2codon[amino] = []
    for codon, amino in codon2amino().items():
        amino2codon[amino].append(codon)

    # For each pair of amino acids, loop over all pairs of codons and compute
    # their minimum Hamming distance.
    min_dist_amino = {'**': 0}
    for amino1 in amino2codon.keys():
        for amino2 in amino2codon.keys():
            min_dist = 3
            for codon1 in amino2codon[amino1]:
                for codon2 in amino2codon[amino2]:
                    hamming = 0
                    for i in range(3):
                        if codon1[i] != codon2[i]:
                            hamming += 1
                    min_dist = min(min_dist, hamming)
            min_dist_amino[amino1 + amino2] = min_dist
        min_dist_amino['*' + amino1] = 0
        min_dist_amino[amino1 + '*'] = 0
    return min_dist_amino


def read_data(fname, dataset='esponda2020'):
    """
    Reads data from the given file as a DataFrame, removing all columns except
    the amino acid sequences and counts, grouping and sorting rows by count.

    :param fname: a string file name to read from (should be a .csv or .tsv)
    :param dataset: a string dataset name ('esponda2020' or 'emerson2017')
    :returns: a DataFrame containing amino acid sequences and their counts
    """
    # Determine the file extension and corresponding separator.
    fext = fname.split('.')[-1]
    if fext == 'csv':
        sep = ','
    elif fext == 'tsv':
        sep = '\t'
    else:
        raise ValueError('fname \"' + fname + '\" must be .csv or .tsv')
    
    # Read the data as a DataFrame, renaming columns according to the dataset.
    df = pd.read_csv(fname, sep=sep)
    if dataset == 'esponda2020':
        col_map = {'aminoAcid': 'amino_seq', 'count (templates/reads)': 'count'}
    elif dataset == 'emerson2017':
        col_map = {'amino_acid': 'amino_seq', 'reads': 'count'}
    else:
        raise ValueError('unrecognized dataset \"' + dataset + '\"')
    df = df.rename(columns=col_map)
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


def dbscan(df_amino, max_clusters, min_cluster_size):
    """
    Implements the DBSCAN algorithm to find at most the given maximum number of
    clusters each containing at least the given minimum number of amino acid
    sequences from the amino acid data.

    :param df_amino: a DataFrame containing a column of amino acid sequences
    :param max_num_clusters: an integer maximum number of clusters to find
    :param min_cluster_size: an integer minimum number of sequences per cluster
    :returns: a list of clusters, each of which are lists of sequence indices
    """
    repertoire = dict(zip(df_amino['amino_seq'], df_amino.index))
    repertoire_idxs = set(df_amino.index.values)
    clusters = []

    while len(clusters) < max_clusters and len(repertoire_idxs) > 0:
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


def extract_equivalent_cluster(cluster, cluster_fname, other_fname):
    """
    TODO: Documentation. This is for extracting equivalent groups.

    :param cluster: a cluster TODO
    :param cluster_fname: a string file name TODO
    :param other_fname: a string file name TODO
    :returns: TODO
    """
    # Create both DataFrames and compute their total number of reads.
    cluster_df, other_df = read_data(cluster_fname), read_data(other_fname)
    cluster_count_tot = cluster_df['count'].sum()
    other_count_tot = other_df['count'].sum()

    # Drop all rows except those in the cluster and mark their indices with 'a'.
    cluster_df = cluster_df.loc[cluster]
    cluster_df.index = [str(i) + 'a' for i in range(len(cluster_df))]

    # Combine these marked cluster rows with the other dataset and use DBSCAN to
    # find a cluster in the combined data. Because this combined data includes
    # the original cluster, the new cluster contains the original cluster and
    # neighboring sequences from the other data. Dropping the marked rows leaves
    # behind an "equivalent" cluster of rows from the other data.
    combined_df = pd.concat([cluster_df, other_df])
    equiv_cluster = dbscan(combined_df[['amino_seq']], max_clusters=1, min_cluster_size=30)
    equiv_cluster = [idx for idx in equiv_cluster[0] if str(idx)[-1] != 'a']

    # Compute the number of reads for the original and equivalent clusters'
    # sequences.
    cluster_count = cluster_df['count'].sum()
    equiv_count = other_df.loc[equiv_cluster]['count'].sum()

    return equiv_cluster, other_df, cluster_count_tot, other_count_tot, cluster_count, equiv_count


class cluster_obj:
    """
    TODO: Documentation.
    """

    def __init__(self, desc, seq_idxs):
        self.desc = desc
        self.cluster = seq_idxs
        self.equivalent_groups = {}


    def set_within_metric(self, number_of_amino_in, number_of_amino_out, number_of_amino_tot, number_of_amino_leaked):
        self.number_of_amino_in = number_of_amino_in
        self.number_of_amino_out = number_of_amino_out
        self.number_of_amino_tot = number_of_amino_tot
        self.number_of_amino_leaked = number_of_amino_leaked


    def set_equivalent_group(self, desc, fname, other_fname):
        """
        TODO: Documentation.

        :param desc: a string description of TODO
        :param fname: a string file name of TODO
        :param other_fname: a string file name of TODO
        """
        equiv_cluster, other_df, cluster_count_tot, other_count_tot, cluster_count, equiv_count = \
            extract_equivalent_cluster(self.cluster, fname, other_fname)
        exp = cluster_obj(desc, equiv_cluster)

        exp.set_within_metric(*evaluate_within(other_df, equiv_cluster)) #comented for other data
        self.equivalent_groups[other_fname] = {'cluster_obj': exp, \
                                                   'sizes': (other_count_tot, cluster_count_tot, equiv_count,
                                                             cluster_count),
                                                   'relative_sizes': (1.0 * equiv_count / other_count_tot,
                                                                      1.0 * cluster_count / cluster_count_tot)}


class experiment:
    """
    TODO: Documentation.
    """

    def __init__(self, description, fname):
        self.description = description
        self.fname = fname
        self.exp_number = 0
        self.experiment_log = {}
        self.explored_set = []


    def add_experiment(self, cluster, df=''):
        self.explored_set += cluster
        self.exp_number += 1
        exp = cluster_obj(self.description, cluster)
        if len(df) == 0:
            df = read_data(self.fname)
        exp.set_within_metric(*evaluate_within(df, cluster))  # method outside class for now
        self.experiment_log[self.exp_number] = exp


    def set_equivalent_groups(self, description, other_file_name, groups=[]):
        """
        TODO: Only used in unused functions (commented out in a used function).
        """
        for g in groups:
            if g in self.experiment_log.keys():
                self.experiment_log[g].set_equivalent_group(description, self.fname, other_file_name)

    
    def save(self, fname):
        with open(fname, 'wb') as f:
            pickle.dump(self, f)


    def find_clusters(self, cluster_frac):
        """
        TODO: Documentation.
        """
        df = read_data(self.fname)
        min_cluster_size = max(int(cluster_frac * len(df)), 5)  # Clusters must have at least 5 sequences.
        num_clusters = 1

        while num_clusters < 5000: ##maximim number of clusters to find
            df = select_data_subset(df, self.explored_set)

            if len(df) < min_size: #should be min size
                print('done')
                break
            grupos = dbscan(df[['aminoAcid']], max_clusters=1, min_cluster_size=min_size)
            if len(grupos) == 0:
                break
            self.add_experiment(grupos[0],df)  # avoid reading I think
            '''
            experiment1.set_equivalent_groups('Using nucleotide sequences, distance=1' + exp_number,
                                            path + 'Spleen_2.tsv', [experiment1.exp_number])
            experiment1.set_equivalent_groups('Using nucleotide sequences, distance=1' + exp_number,
                                            path + 'Spleen_3.tsv', [experiment1.exp_number])
            '''
            num_clusters += 1
            #print('Iteration = ', cluster_num, ' Clusters processed = ', experiment1.exp_number)
        name = self.fname.split('/')[-1].split('.')[0] + 'Experiment' + self.exp_number + '-' + datetime.datetime.fromtimestamp(time()).strftime(
                '%Y-%m-%d-%H')
        self.save(name)


# Calculate the  aminoacid 'leakage'
def evaluate_within(df, cluster):
    if 'vMaxResolved' not in df.columns:
        return 1, 1, 1, 1
    indexes_not_in_group = set(df.index.values).difference(set(cluster))

    all_data = df['aminoAcid'].unique()
    in_group = df.loc[cluster]['aminoAcid'].unique()
    ##filter v segments
    v_segments = df.loc[cluster]['vMaxResolved'].unique()
    ###filter j segments
    j_segments = df.loc[cluster]['jMaxResolved'].unique()
    # print('I havent included j because of testing')
    not_in_group = df.loc[
        df.index.isin(indexes_not_in_group) & df['vMaxResolved'].isin(v_segments) & df['jMaxResolved'].isin(
            j_segments)]['aminoAcid'].unique()

    # not_in_group=df.loc[df.index.isin(indexes_not_in_group) & df['vMaxResolved'].isin(v_segments)]['aminoAcid'].unique()

    amino_leaked = set(in_group).intersection(set(not_in_group))

    number_of_amino_in = len(in_group)
    number_of_amino_out = len(not_in_group)
    number_of_amino_tot = len(all_data)
    number_of_amino_leaked = len(amino_leaked)
    if (len(cluster) > 0):
        print(1.0 * number_of_amino_leaked / number_of_amino_in)
    print(number_of_amino_in, number_of_amino_out, number_of_amino_tot, number_of_amino_leaked)
    return number_of_amino_in, number_of_amino_out, number_of_amino_tot, number_of_amino_leaked


def load_data_object(name,dataPath='',experimentPath=''):
    """
    TODO: Unused, but clearly the complement of experiment.save().
    """
    name=experimentPath+name.split('/')[-1]
    with open(name,'rb') as fp:
        data_object=pickle.load(fp)
    data_object.fname = dataPath + data_object.fname.split('/')[-1]
    return data_object


def select_data_subset(df,explored_indexes):
    indexes_not_explored=set(df.index.values).difference(set(explored_indexes))
    return df.loc[indexes_not_explored]


if __name__ == "__main__":
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-E', '--exps', type=str, default='exps/', \
                        help='path to experiments directory')
    parser.add_argument('-D', '--data', type=str, default='data/', \
                        help='path to data directory')
    parser.add_argument('-F', '--file', type=str, default='all', \
                        help='TODO: what is this?')
    parser.add_argument('-S', '--size', type=float, default=0.001, \
                        help='Minimum cluster size as fraction of repertoire')
    args = parser.parse_args()

    if args.file == 'all':  # Running all files in path.
        files = [f for f in listdir(args.data) if isfile(join(args.data, f))]
        for file in files:
            exp = experiment('Using amino acid sequences, distance=1', join(args.data, file))  # TODO: better description.
            exp.exp_number = '1'
            exp.find_clusters(args.size)
    else:  # Running one file.
        exp = experiment('Using amino acid sequences, distance=1', join(args.data, args.file))
        exp.exp_number = '1'
        exp.find_clusters(args.size)
