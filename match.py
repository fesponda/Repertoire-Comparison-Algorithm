# Project:   PublicClusters
# Filename:  match.py
# Authors:   Fernando Esponda (fernando.esponda@itam.mx) and Joshua J. Daymude
#            (jdaymude@asu.edu).

from cluster import *


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


def clusterDistance(c1,c2,inter=''):
    """
    TODO
    """
    if inter=='':
        inter=set()
    ##Check for sequences not included
    #if (len(c2)<len(c1)):
    #    temp=c1
    #    c1=c2
    #    c2=temp
    s1=set()
    for amino1 in c1: #one away from c1
        s1 |= set(generate_nbr_seqs(amino1))
    s1=s1.difference(c1)
    s2=set()
    for amino1 in c2: #one away from c2
        s2 |= set(generate_nbr_seqs(amino1))
    s2 = s2.difference(c2)

    inter |= s1.intersection(s2) #number of sequences 1 away from both clusters

    return inter


def has_nbr_in_cluster(amino_seq, repertoire, dist_metric='Hamming'):
    """
    Returns True iff there is an amino acid sequence in the given repertoire
    that is within distance 1 of the given amino acid sequence.

    :param amino_seq: a string representing an amino acid sequence
    :param repertoire: a dict mapping amino acid sequences to their indices
    :param dist_metric: a string distance metric ('Hamming' or 'Levenshtein')
    :returns: True iff the cluster has a sequence within distance 1
    """
    if dist_metric == 'Hamming':
        for i in range(len(amino_seq)):
            for a in amino_acids():
                # Substitute the i-th amino acid with amino acid a.
                if amino_seq[0:i] + a + amino_seq[i+1:] in repertoire.keys():
                    return True
    elif dist_metric == 'Levenshtein':
        for i in range(len(amino_seq)):
            for a in amino_acids():
                # Substitute the i-th amino acid with amino acid a or insert
                # amino acid a just before the i-th amino acid.
                if amino_seq[0:i] + a + amino_seq[i+1:] in repertoire.keys() or\
                   amino_seq[0:i] + a + amino_seq[i:] in repertoire.keys():
                    return True
            # Delete the i-th amino acid.
            if amino_seq[0:i] + amino_seq[i+1:] in repertoire.keys():
                return True
        for a in amino_acids():
            # Insert amino acid a at the end of the sequence.
            if amino_seq + a in repertoire.keys():
                return True

    return False


def is_match(df1, df2, cluster1, cluster2, exact=True):
    """
    Determines whether there exists an amino acid sequence in cluster1 that
    matches a sequence in cluster2, where "matches" is either exact or within
    distance 1.

    :param df1: a DataFrame containing a column of amino acid sequences
    :param df2: a DataFrame containing a column of amino acid sequences
    :param cluster1: a list of amino acid sequence indices from df1
    :param cluster2: a list of amino acid sequence indices from df2
    :param exact: True iff a sequence must appear exactly to match
    :returns: True iff a sequence in cluster1 matches a sequence in cluster2
    """
    cluster2 = set(df2.loc[cluster2]['amino_seq'])

    for seq_idx in cluster1:
        if exact:
            if df1.amino_seq[seq_idx] in cluster2:
                return True
        elif has_nbr_in_cluster(df1.amino_seq[seq_idx], cluster2):
            return True

    return False


def match_clusters_slow(exp1, exp2):
    """
    [DO NOT USE] Computes which clusters of exp2 match those of exp1. Very
    inefficient due to repeatedly calculating whether a sequence has a neighbor
    in a given cluster. This information was already calculated in close_to
    during clustering and should be reused.

    :param exp1: an Experiment object containing the first set of clusters
    :param exp2: an Experiment object containing the second set of clusters
    :returns: a dict mapping TODO to TODO and a list of clusters present in exp1
    that are missing in exp2.
    """
    matches = {}
    missing = []
    df1 = read_data(exp1.fname)
    df2 = read_data(exp2.fname)
    clusters1 = exp1.experiment_log.keys()
    clusters2 = exp2.experiment_log.keys()

    if len(set(exp1.explored_set) - set(df1.index)) > 0 or \
       len(set(exp2.explored_set) - set(df2.index)) > 0:
        print('Index mismatch')
        return matches, missing

    for c1 in clusters1:
        matches[c1] = []
        is_missing = True
        for c2 in clusters2:
            if is_match(df1, df2, exp1.experiment_log[c1].group_indexes, \
                        exp2.experiment_log[c2].group_indexes, exact=False):
                is_missing = False
                matches[c1] += [c2]
        if is_missing:
            missing += [c1]

    return matches, missing


def match_clusters(exp1, exp2):
    """
    Computes which clusters of exp2 match those of exp1.

    TODO: Can be further optimized by ordering clusters by size and sequences
    by degree.

    :param exp1: an Experiment object containing the first set of clusters
    :param exp2: an Experiment object containing the second set of clusters
    :returns: a dict mapping TODO to TODO and a list of clusters present in exp1
    that are missing in exp2.
    """
    matches = {}
    missing = []
    df1 = read_data(exp1.fname)
    df2 = read_data(exp2.fname)    
    clusters1 = exp1.experiment_log.keys()
    clusters2 = exp2.experiment_log.keys()

    if len(set(exp1.explored_set) - set(df1.index)) > 0 or \
       len(set(exp2.explored_set) - set(df2.index)) > 0:
        print('Index mismatch')
        return matches, missing

    # Also not sure about this?
    clusterDict = {}
    for c2 in clusters2:
        clusterDict[c2] = set(df2.loc[exp2.experiment_log[c2].group_indexes]['aminoAcid'])

    for c1 in clusters1:
        matches[c1] = []
        is_missing = True
        clusters2 = set(exp2.experiment_log.keys())
        out = set()

        for sequence1 in exp1.experiment_log[c1].group_indexes:
            if len(clusters2)==0:
                break
            aminoString = df1['aminoAcid'][sequence1]
            variations = set(generate_nbr_seqs(aminoString))
            for c2 in clusters2:
                common = variations.intersection(clusterDict[c2])
                if len(common) > 0:  # found
                    is_missing = False
                    matches[c1] += [c2]
                    out |= set([c2])
            clusters2 = clusters2.difference(out)

        if is_missing:
            missing.append(c1)

    return matches, missing


def matchExperimentsBack(files):
    clusterMatches={}
    ## There is some redundant computation (because of symetry) but I leave them now for a sanity check
    for file1 in files:
        experiment1=load_data_object(file1)
        #experiment1.fname = path + experiment1.fname.split('/')[-1]
        clusterMatches[file1]={}
        clusterMatches[file1]['size']=len(experiment1.explored_set)##assumes not repeated sequences
        for file2 in files:
            if file1 != file2:
                experiment2=load_data_object(file2)
               # experiment2.fname = path + experiment2.fname.split('/')[-1]
                clusterMatches[file1][file2]={}

                graph,missing=match_clusters(experiment1,experiment2)
                clusterMatches[file1][file2]['graph']=graph
                clusterMatches[file1][file2]['missing']=missing
                with open('matches.pickle', 'wb') as f:
                    pickle.dump(clusterMatches, f)

    return clusterMatches


def invertGraph(graph,groups):
    graph2={}
    for v in graph:
        for v2 in graph[v]:
            if v2 in graph2.keys():
                graph2[v2]+=[v]
            else:
                graph2[v2]=[v]
    for v in graph2:
        graph2[v]=list(set(graph2[v]))
    return graph2, list(set(groups).difference(set(graph2.keys())))


def matchExperiments(files,dataPath='',experimentPath=''):
    clusterMatches={}
    ## There is some redundant computation (because of symetry) but I leave them now for a sanity check
    for x,file1 in enumerate(files):

        experiment1=load_data_object(file1,dataPath=dataPath,experimentPath=experimentPath)
        #experiment1.fname = path + experiment1.fname.split('/')[-1]
        clusterMatches[file1]={}
        for y,file2 in enumerate(files):
            experiment2 = load_data_object(file2,dataPath=dataPath,experimentPath=experimentPath)
            #experiment2.fname = path + experiment2.fname.split('/')[-1]

            if x < y:
                clusterMatches[file1][file2] = {}
                graph,missing=match_clusters(experiment1,experiment2)
                clusterMatches[file1][file2]['graph'] = graph
                clusterMatches[file1][file2]['missing'] = missing
                with open(experimentPath+'matches.pickle', 'wb') as f:
                    pickle.dump(clusterMatches, f)
            elif x > y:
                clusterMatches[file1][file2] = {}
                graph,missing=invertGraph(clusterMatches[file2][file1]['graph'],experiment1.experiment_log.keys())
                clusterMatches[file1][file2]['graph'] = graph
                clusterMatches[file1][file2]['missing'] = missing
                with open(experimentPath+'matches.pickle', 'wb') as f:
                    pickle.dump(clusterMatches, f)


    return clusterMatches


def doubleTraverse(G1, G2, v, nodesG1, nodesG2, visitadoG1, visitadoG2, which):
    if which:
        visitadoG1[v] = True
        nodesG1 += [v]
        for v2 in G1[v]:
            if not visitadoG2[v2]:
                doubleTraverse(G1, G2, v2, nodesG1, nodesG2, visitadoG1, visitadoG2, not which)
    else:
        visitadoG2[v] = True
        nodesG2 += [v]
        for v2 in G2[v]:
            if not visitadoG1[v2]:
                doubleTraverse(G1, G2, v2, nodesG1, nodesG2, visitadoG1, visitadoG2, not which)


def consolidateGraph(g1,g2):
    visitadoG1={}
    for v in g1:
        visitadoG1[v]=False
    visitadoG2={}
    for v in g2:
        visitadoG2[v]=False

    multiGraph={}
    for v in visitadoG1:
        nodesG1=[]
        nodesG2=[]
        if not visitadoG1[v]:
            doubleTraverse(g1,g2,v,nodesG1,nodesG2,visitadoG1,visitadoG2,True)
            multiGraph[v]={'g1':nodesG1,'g2':nodesG2}
    return multiGraph


# We assume unrepeated strings. Other wise we would need to read the file
'''
a1=read_data(experiment1.fname)
s=set(a1.loc[experiment1.experiment_log[i].group_indexes]['aminoAcid'])
'''
def group_size(experiment1,groups,relative=True):
    tam=0
    for i in groups:
        tam+=len(experiment1.experiment_log[i].group_indexes)
    if relative:
        return 100*tam/len(experiment1.explored_set)
    return tam


if __name__ == "__main__":
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-E', '--exps', type=str, default='exps/', \
                        help='path to experiments directory')
    parser.add_argument('-D', '--data', type=str, default='data/', \
                        help='path to data directory')
    args = parser.parse_args()
    
    files = [f for f in listdir(args.exps) if 'Experiment' in f]
    clusterMatches = matchExperiments(files, args.data, args.exps)