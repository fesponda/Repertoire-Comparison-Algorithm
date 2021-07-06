# Project:  PublicClusters
# Filename: analyze.py
# Authors:  Fernando Esponda (fernando.esponda@itam.mx) and Joshua J. Daymude
#           (jdaymude@asu.edu).

"""
analyze: Compute metrics regarding repertoire clusters and matches.
"""

# TODO: From the paper, it seems that this file should contain functions to get:
# Figure 4: Earth Mover's Distance (EMD)
# Figure 5: Anomaly Detection using EMD

from match import *


def traverse_matches(matches1, matches2, cluster, clusters1, clusters2, \
                     visited1, visited2, forward):
    """
    A recursive helper function for match_components that performs depth-first
    search (DFS) on the union of two symmetric, directed, bipartite graphs
    representing the cluster matches between two individuals (ID1 <-> ID2).

    :param matches1: a dict mapping cluster indices of ID1 to the indices of
    clusters in ID2 that they match with
    :param matches2: a dict mapping cluster indices of ID1 to the indices of
    clusters in ID2 that they match with
    :param cluster: a cluster index from which to start the traversal
    :param clusters1: a list of visited clusters indices of ID1
    :param clusters2: a list of visited clusters indices of ID2
    :param visited1: a dict mapping clusters indices of ID1 to True iff visited
    :param visited2: a dict mapping clusters indices of ID2 to True iff visited
    :param forward: True iff the next edge to traverse is from ID1 to ID2
    """
    if forward:  # Traverse a forward edge from cset1 to cset2.
        visited1[cluster] = True
        clusters1 += [cluster]
        for match in matches1[cluster]:
            if not visited2[match]:
                traverse_matches(matches1, matches2, match, clusters1, \
                                 clusters2, visited1, visited2, not forward)
    else:  # Traverse a backward edge from cset2 to cset1.
        visited2[cluster] = True
        clusters2 += [cluster]
        for match in matches2[cluster]:
            if not visited1[match]:
                traverse_matches(matches1, matches2, match, clusters1, \
                                 clusters2, visited1, visited2, not forward)


def match_components(matches1, matches2):
    """
    Finds the connected components in the union of two symmetric, directed,
    bipartite graphs representing the cluster matches from one individual (ID1)
    to another (ID2) and vice versa.

    :param matches1: a dict mapping cluster indices of ID1 to the indices of
    clusters in ID2 that they match with
    :param matches2: a dict mapping cluster indices of ID1 to the indices of
    clusters in ID2 that they match with
    :returns: a list of tuples representing the connected components; the first
    (resp., second) member of each tuple is a list of indices of the clusters in
    the ID1 (resp., ID2) in that component
    """
    visited1, visited2 = {}, {}  # Map cluster indices to True/False.
    for cluster in matches1:
        visited1[cluster] = False
    for cluster in matches2:
        visited2[cluster] = False

    components = []
    for cluster in visited1:
        clusters1, clusters2 = [], []
        if not visited1[cluster]:
            traverse_matches(matches1, matches2, cluster, clusters1, clusters2,\
                             visited1, visited2, forward=True)
            components.append((clusters1, clusters2))

    return components


def cluster_size(cset, clusters, relative=True):
    """
    Calculates the fraction (if relative) or raw number (otherwise) of amino
    acid sequences in the specified clusters.

    :param cset: a ClusterSet object
    :param clusters: a list of indices of cluters to measure
    :returns: the fraction/number of amino acid sequences in the given clusters
    """
    size = sum([len(cset.clusters[i]) for i in clusters])
    return size / len(cset.clustered_seqs) if relative else size


def cluster_size_pearson(dataset):
    """
    Calculates the minimum Pearson correlation between pairs of cluster sizes of
    all pairs of individuals in the given dataset.

    :param dataset: a string dataset name
    :returns: a dict mapping cluster file names to dicts mapping cluster file
    names to the minimum Pearson correlation between pairs of cluster sizes for
    this pair of cluster files
    """
    pearson = {}

    # Load the matches dict from file.
    with open(osp.join('matches', dataset + '.pkl'), 'rb') as f:
        matches = pickle.load(f)

    for f_i in matches:
        cset_i = load_obj(osp.join('clusters', dataset, f_i))
        pearson[f_i] = {f_i : 1.0}  # Perfectly correlated with itself.

        for f_j in matches[f_i]:
            cset_j = load_obj(osp.join('clusters', dataset, f_j))
            matches_ij = matches[f_i][f_j]['matches']
            matches_ji = matches[f_j][f_i]['matches']

            # Get the connected components in the matching graph.
            components = match_components(matches_ij, matches_ji)

            # Calculate the sizes of the i and j clusters in each component.
            sizes_i, sizes_j = [], []
            for component in components:
                sizes_i.append(cluster_size(cset_i, component[0]))
                sizes_j.append(cluster_size(cset_j, component[1]))

            # Record the minimum Pearson correlation between i and j clusters
            # of any component.
            pearson[f_i][f_j] = np.min(np.corrcoef(sizes_i, sizes_j))

    return pearson


def cluster_stats(dataset):
    """
    Compute cluster- and sequence-level statistics for the given dataset,
    including the number of clusters, number of (un)clustered sequences,
    and pairwise fractions of missing sequences and clusters.

    :param dataset: a string dataset name
    :returns: a dict mapping cluster file names to their numbers of clusters
    :returns: a dict mapping cluster file names to dicts containing their number
    of clustered and unclustered sequences
    :returns: a dict mapping cluster file names to dicts mapping cluster file
    names to the fraction of missing sequences for this pair
    :returns: a dict mapping cluster file names to dicts mapping cluster file
    names to the fraction of missing clusters for this pair
    """
    num_clusters = {}      # Map cluster files to # clusters.
    num_seqs = {}          # Map cluster files to # sequences (un)clustered.
    seqs_missing = {}      # Map pairs of cluster files to the fraction of
                           # sequences in clusters_i missing from clusters_j.
    clusters_missing = {}  # Map pairs of cluster files to the fraction of
                           # clusters in clusters_i missing from clusters_j.

    # Load the matches dict from file.
    with open(osp.join('matches', dataset + '.pkl'), 'rb') as f:
        matches = pickle.load(f)

    for f_i in matches:
        cset_i = load_obj(osp.join('clusters', dataset, f_i))
        num_clusters[f_i] = len(cset_i.clusters)
        num_seqs[f_i] = {'clustered': len(cset_i.clustered_seqs), \
                         'unclustered': len(cset_i.seq_df) \
                                        - len(cset_i.clustered_seqs)}
        seqs_missing[f_i] = {f_i : 0}
        clusters_missing[f_i] = {f_i : 0}

        for f_j in matches[f_i]:
            seqs_missing[f_i][f_j] = sum([len(cset_i.clusters[miss]) \
                for miss in matches[f_i][f_j]['missing']])
            seqs_missing[f_i][f_j] /= len(cset_i.clustered_seqs)

            clusters_missing[f_i][f_j] = \
                len(matches[f_i][f_j]['missing']) / len(cset_i.clusters)

    return num_clusters, num_seqs, seqs_missing, clusters_missing