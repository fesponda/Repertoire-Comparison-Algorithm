# Project:   PublicClusters
# Filename:  analysis.py
# Authors:   Fernando Esponda (fernando.esponda@itam.mx) and Joshua J. Daymude
#            (jdaymude@asu.edu).

from match import *

import numpy as np


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


def traverse_matches(matches1, matches2, cluster, clusters1, clusters2, \
                     visited1, visited2, first):
    """
    TODO

    :param matches1: a dict mapping clusters (i.e., lists of amino acid sequence
    indices) to lists of matching clusters
    :param matches2: a dict mapping clusters to lists of matching clusters
    :param cluster: a cluster from which to start the traversal
    :param clusters1: a list of clusters in matches1 that have been visited
    :param clusters2: a list of clusters in matches2 that have been visited
    :param visited1: a dict mapping clusters of matches1 to True iff visited
    :param visited2: a dict mapping clusters of matches2 to True iff visited
    :param first: True iff the next matches graph to traverse is matches1
    """
    if first:  # Traverse matches1.
        visited1[cluster] = True
        clusters1 += [cluster]
        for match in matches1[cluster]:
            if not visited2[match]:
                traverse_matches(matches1, matches2, match, clusters1, \
                                 clusters2, visited1, visited2, not first)
    else:  # Traverse matches2.
        visited2[cluster] = True
        clusters2 += [cluster]
        for match in matches2[cluster]:
            if not visited1[match]:
                traverse_matches(matches1, matches2, match, clusters1, \
                                 clusters2, visited1, visited2, not first)


def consolidate_matches(matches1, matches2):
    """
    TODO: Documentation.

    :param matches1: TODO
    :param matches2: TODO
    :returns: TODO
    """
    visited1, visited2 = {}, {}
    for cluster in matches1:
        visited1[cluster] = False
    for cluster in matches2:
        visited2[cluster] = False

    multiGraph={}
    for cluster in visited1:
        clusters1, clusters2 = [], []
        if not visited1[cluster]:
            traverse_matches(matches1, matches2, cluster, clusters1, \
                             clusters2, visited1, visited2, first=True)
            multiGraph[cluster] = {'g1': clusters1, 'g2': clusters2}
    return multiGraph


# Compute pairwise correlation in size
def pairwiseCorrelations(clusterMatches, dataPath='', experimentPath=''):
    result = {}

    minClusterSizePercentage = 0
    for f1 in clusterMatches.keys():

        # if('human'in f1):
        #     continue
        result[f1] = {}
        result[f1][f1] = 1.0
        # clusterMatches[f1].pop('size', None)
        experiment1 = load_data_object(f1, dataPath, experimentPath)

        for f2 in clusterMatches[f1].keys():
            #   if('human'in f2):
            #       continue
            g1 = clusterMatches[f1][f2]['graph']
            g2 = clusterMatches[f2][f1]['graph']

            multiGraph = consolidate_matches(g1, g2)
            x = []
            y = []
            experiment2 = load_data_object(f2, dataPath, experimentPath)
            for mNode in multiGraph:
                s1 = group_size(experiment1, multiGraph[mNode]['g1'])
                s2 = group_size(experiment2, multiGraph[mNode]['g2'])
                if s1 > minClusterSizePercentage and s2 > minClusterSizePercentage:
                    x += [s1]
                    y += [s2]
            result[f1][f2] = np.min(np.corrcoef(x, y))
            # row.append((f1,f2,np.min(np.corrcoef(x,y))))
    return result


# size of missing clusters
def missingClusterSizes(clusterMatches, dataPath='', experimentPath=''):
    totalMissing = {}
    numberClusterMissing = {}
    clusteredUnClustered = {}
    numberOfGroups = {}
    for f1 in clusterMatches.keys():

        experiment1 = load_data_object(f1, dataPath, experimentPath)
        totalMissing[f1] = {}
        numberClusterMissing[f1] = {}
        numberOfGroups[f1] = experiment1.exp_number
        numberClusterMissing[f1][f1] = 0
        totalMissing[f1][f1] = 0
        size=len(read_data(experiment1.file_name).index)
        clustered = len(experiment1.explored_set)
        unclustered =  size - clustered
        clusteredUnClustered[f1] = {'clustered': clustered, 'unclustered': unclustered,'size':size}

        for f2 in clusterMatches[f1].keys():
            # print(f1,f2)

            sum = 0.0
            for miss in clusterMatches[f1][f2]['missing']:
                sum += 100 * len(experiment1.experiment_log[miss].group_indexes) / len(experiment1.explored_set)
                # print(miss,100*len(experiment1.experiment_log[miss].group_indexes)/len(experiment1.explored_set))
            totalMissing[f1][f2] = sum
            if experiment1.exp_number!=0:
                numberClusterMissing[f1][f2] = 100 * len(clusterMatches[f1][f2]['missing']) / experiment1.exp_number
                #numberClusterMissing[f1][f2] = len(clusterMatches[f1][f2]['missing'])

            else:
                print('No clusters')
                numberClusterMissing[f1][f2]=0

            # print('Total missing>:',sum)
    return totalMissing, numberClusterMissing, clusteredUnClustered, numberOfGroups


