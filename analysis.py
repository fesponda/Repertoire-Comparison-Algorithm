from scipy.stats import spearmanr, stats
from scipy.spatial.distance import euclidean
from random import randint, shuffle,sample,random
from os import listdir
from os.path import isfile, join
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from time import time
import pickle
import networkx as nx
from math import isinf, exp
import random
import sys
import json
sys.path.append("./code/")
from cluster import *
from matches import *

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

            multiGraph = consolidateGraph(g1, g2)
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


