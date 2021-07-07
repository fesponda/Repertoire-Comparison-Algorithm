# Project:  PublicClusters
# Filename: match.py
# Authors:  Fernando Esponda (fernando.esponda@itam.mx) and Joshua J. Daymude
#           (jdaymude@asu.edu).

"""
match: Match clusters from distinct individuals.
"""

from cluster import *

from collections import defaultdict
from itertools import combinations


def generate_nbr_seqs(amino_seq, dist_metric='Hamming'):
    """
    Generates all amino acid sequences within distance 1 of the given sequence.

    :param amino_seq: a string representing an amino acid sequence
    :param dist_metric: a string distance metric ('Hamming' or 'Levenshtein')
    :returns: a list of strings representing sequences within distance 1
    """
    nbr_seqs = []
    if dist_metric == 'Hamming':
        for i in range(len(amino_seq)):
            for a in amino_acids():
                # Substitute the i-th amino acid with amino acid a.
                nbr_seqs.append(amino_seq[0:i] + a + amino_seq[i+1:])
    elif dist_metric == 'Levenshtein':
        for i in range(len(amino_seq)):
            for a in amino_acids():
                # Substitute the i-th amino acid with amino acid a.
                nbr_seqs.append(amino_seq[0:i] + a + amino_seq[i+1:])
                # Insert amino acid a just before the i-th amino acid.
                nbr_seqs.append(amino_seq[0:i] + a + amino_seq[i:])
            # Delete the i-th amino acid.
            nbr_seqs.append(amino_seq[0:i] + amino_seq[i+1:])
        for a in amino_acids():
            # Insert amino acid a at the end of the sequence.
            nbr_seqs.append(amino_seq + a)

    return nbr_seqs


def match_clusters(cset1, cset2):
    """
    Computes which clusters of cset2 match those of cset1.

    TODO: Can be further optimized by ordering clusters by size and sequences
    by degree.

    :param cset1: a ClusterSet object containing the first set of clusters
    :param cset2: a ClusterSet object containing the second set of clusters
    :returns: a dict mapping cluster indices of cset1 to the indices of clusters
    in cset2 that they match with
    :returns: a list of indices of clusters in cset1 that cset2 is missing.
    """
    matches, missing = {}, []

    # For efficiency, retrieve and store the amino acid sequences in each
    # cluster of cset2 instead of repeatedly looking them up.
    cluster_seqs2 = {}
    for j, cluster_j in enumerate(cset2.clusters):
        cluster_seqs2[j] = set(cset2.seq_df.loc[cluster_j]['amino_seq'])

    for i, cluster_i in enumerate(cset1.clusters):
        matches[i] = []
        is_missing = True
        matched_idxs2 = set()  # Clusters from cset2 that matched cluster_i.

        # For efficiency, we generate the neighbors of each amino acid sequence
        # exactly once and then check for matching sequences in all clusters of
        # cset2 in one go.
        for seq_idx_i in cluster_i:
            seq_i = cset1.seq_df['amino_seq'][seq_idx_i]
            nbr_seqs_i = set(generate_nbr_seqs(seq_i))

            # For futher efficiency, clusters of cset2 that already matched with
            # an earlier sequence in cluster_i aren't checked again.
            for j in set(range(len(cset2.clusters))) - matched_idxs2:
                if len(nbr_seqs_i & cluster_seqs2[j]) > 0:
                    is_missing = False
                    matches[i].append(j)
                    matched_idxs2 |= set([j])

        if is_missing:
            del(matches[i])
            missing.append(i)

    return matches, missing


def invert_matches(matches, cluster_idxs2):
    """
    If matches is the result of match_clusters(cset1, cset2) and cluster_idxs2
    is a list of cluster indices in cset2, then this function returns a matches
    dict and missing list as if match_clusters(cset2, cset1) was called.

    :param matches: a dict mapping cluster indices to lists of cluster indices
    that they match with
    :param clusters: a list of cluster indices representing the range of matches
    :returns: a dict mapping cluster indices from the range of matches to the
    cluster indices fromthe domain of matches that matched with them
    :returns: a list of indices of clusters in the domain of matches that its
    range is missing
    """
    matches_inv = {}
    for cluster in matches:
        for match in matches[cluster]:
            if match in matches_inv.keys():
                matches_inv[match].append(cluster)
            else:
                matches_inv[match] = [cluster]
    for match in matches_inv:
        matches_inv[match] = list(set(matches_inv[match]))

    return matches_inv, list(set(cluster_idxs2) - set(matches_inv.keys()))


def match_pairs(dataset, file_pairs, queue):
    """
    Computes matching clusters for each of the given pairs of cluster files and
    enqueues the results in the given multiprocessing queue.

    :param dataset: a string dataset name
    :param file_pairs: a list of pairs of cluster file names to match
    :param queue: a multiprocessing.Queue for storing the resulting matches
    """
    matches = defaultdict(dict)
    for f_i, f_j in tqdm(file_pairs, desc='PID{} Matching'.format(os.getpid())):
        cset_i = load_obj(osp.join('clusters', dataset, f_i))
        cset_j = load_obj(osp.join('clusters', dataset, f_j))
        matches_ij, missing_ij = match_clusters(cset_i, cset_j)
        matches[f_i][f_j] = {'matches': matches_ij, 'missing': missing_ij}

    queue.put(matches)


def match_all(dataset, num_procs=1):
    """
    Computes matching clusters between all pairs of cluster files in the given
    dataset's clusters directory and writes the results to file.

    :param dataset: a string dataset name
    :param num_procs: an int number of processors to parallelize over
    :returns: a dict mapping cluster file names to dicts mapping cluster file
    names to dicts containing the matching and missing clusters for this pair
    """
    # Partition pairs of file indices to match on over the number of processors.
    files = os.listdir(osp.join('clusters', dataset))
    file_pairs = [pair for pair in combinations(files, r=2)]
    file_pair_chunks = np.array_split(file_pairs, num_procs)

    # Start all processes, performing matching over all unique pairs.
    procs = []
    queue = mp.Queue()
    for chunk in file_pair_chunks:
        proc = mp.Process(target=match_pairs, args=(dataset, chunk, queue,))
        proc.start()
        procs.append(proc)

    # Collect the returned matches as the processes finish.
    chunk_matches = []
    for proc in procs:
        chunk_matches.append(queue.get())
        proc.join()

    # Merge all returned matches into one dict.
    matches = defaultdict(dict)
    for cmatches in chunk_matches:
        for f_i in cmatches:
            for f_j in cmatches[f_i]:
                matches[f_i][f_j] = cmatches[f_i][f_j]

    # For efficiency, use symmetry and inversion to get the remaining matches.
    for i, f_i in enumerate(files):
        cset_i = load_obj(osp.join('clusters', dataset, files[i]))
        for j, f_j in enumerate(files):
            if i > j:
                matches_ji = matches[f_j][f_i]['matches']
                clusters_i = list(range(len(cset_i.clusters)))
                matches_ij, missing_ij = invert_matches(matches_ji, clusters_i)
                matches[f_i][f_j] = {'matches':matches_ij, 'missing':missing_ij}

    # Write the matches to file before returning them.
    with open(osp.join('matches', dataset + '.pkl'), 'wb') as f:
        pickle.dump(matches, f)

    return matches


if __name__ == "__main__":
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-D', '--dataset', type=str, required=True, \
                        help='Dataset to perform pairwise matching on')
    parser.add_argument('-P', '--num_procs', type=int, default=1, \
                        help='Number of processors to parallelize over')
    args = parser.parse_args()

    match_all(args.dataset, args.num_procs)
