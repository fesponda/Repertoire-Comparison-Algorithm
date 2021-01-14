from scipy.stats import spearmanr, stats
from scipy.spatial.distance import euclidean
from random import randint, shuffle,sample,random
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import datetime
from time import time
import pickle
import networkx as nx
from math import isinf, exp
import random
import os
import sys
import json
sys.path.append("./code/")
from cluster import *
#path='/Users/fesponda/research/Immunology/Data/MissingVJ/'
#datapath='/Users/fesponda/research/Immunology/Data/MissingVJ/'

#experimentPath='/Users/fesponda/research/Immunology/Data/experiments/MissingVJ/'
#experimentPath='~/Immunology/Data/experiments/'


def generateVariations(aminoString, distanceMetric='Hamming'):
    current_matches = []
    if distanceMetric=='Hamming':
        for i in range(len(aminoString)):
            for j in dnaValues:
                amino = aminoString[0:i] + j + aminoString[i + 1:]
                current_matches += [amino]
    elif distanceMetric=='Edit':
        for i in range(len(aminoString)):
            for j in dnaValues:
                if  True or  min_dist_amino[aminoString[i] + j] <= 1: #True  or //strict aminoacid edit, no codon
                    amino = aminoString[0:i] + j + aminoString[i + 1:]
                    current_matches+=[amino]
                amino=aminoString[0:i] + j + aminoString[i:]
                current_matches+=[amino]
            amino=aminoString[0:i] + aminoString[i + 1:]
            current_matches+=[amino]

        for j in dnaValues:  # insertion in the last position
            amino=aminoString + j
            current_matches+=[amino]
    return current_matches

'''
# Only efficient when distance is 1 or 0
'''


def checkVariations(aminoString, cluster2, distanceMetric='Hamming'):
    found = False
    if distanceMetric=='Hamming':
        for i in range(len(aminoString)):
            for j in dnaValues:
                amino = aminoString[0:i] + j + aminoString[i + 1:]
                if amino in cluster2.keys():  ##substitutions
                    found = True
                    break
            if found:
                break
    elif distanceMetric=='Edit':

        for i in range(len(aminoString)):
            for j in dnaValues:
                if True or min_dist_amino[aminoString[i] + j] <= 1:
                    amino = aminoString[0:i] + j + aminoString[i + 1:]  # substitutions
                    if amino in cluster2:
                        found = True
                        break
                amino = aminoString[0:i] + j + aminoString[i:]
                if amino in cluster2.keys():  # insertion
                    found = True
                    break
            amino = aminoString[0:i] + aminoString[i + 1:]
            if not found and amino in cluster2.keys():  # deletion
                found = True

            if found:
                break
        if not found:
            for j in dnaValues:  # insertion in the last position
                amino = aminoString + j
                if amino in cluster2.keys():
                    found = True
                    break

    return found


def checkVariationsBack (aminoString, cluster2):
    found = False
    for i in range(len(aminoString)):
        for j in dnaValues:
            if min_dist_amino[aminoString[i] + j] <= 1:
                amino = aminoString[0:i] + j + aminoString[i + 1:]
                if amino in cluster2:
                    found = True
                    break
        if found:
            break

    return found


def isEquivalent(df1, df2, group1, group2, exact=True):
    equivalent = False
    #cluster2 = set(df2.loc[df2.index.intersection(group2)]['aminoAcid'])
    cluster2=set(df2.loc[group2]['aminoAcid'])

    for sequence1 in group1:
        if exact:
            if df1['aminoAcid'][sequence1] in cluster2:
                equivalent=True

                break
        elif checkVariations(df1['aminoAcid'][sequence1], cluster2):
            equivalent = True
            break

    return equivalent

def clusterDistance(c1,c2,inter=''):
    if inter=='':
        inter=set()
    ##Check for sequences not included
    #if (len(c2)<len(c1)):
    #    temp=c1
    #    c1=c2
    #    c2=temp
    s1=set()
    for amino1 in c1: #one away from c1
        s1 |= set(generateVariations(amino1))
    s1=s1.difference(c1)
    s2=set()
    for amino1 in c2: #one away from c2
        s2 |= set(generateVariations(amino1))
    s2 = s2.difference(c2)

    inter |= s1.intersection(s2) #number of sequences 1 away from both clusters

    return inter

def matchClusters(exp1, exp2):
    ##We first order clusters by size to better performance. Other optimizations (not implemented yet) will include:
    ##Sorting sequences by degree
    ## comparing samples of all groups to decide most likely matches and order accordingly
    # clusterOrder1,clusterOrder2=orderClusters(exp1,exp2)

    result = []
    missing = []
    graph = {}
    df1 = read_data(exp1.file_name)
    df2 = read_data(exp2.file_name)

    # groups1=sorted(exp1.experiment_log.keys(),key=lambda x: len(exp1.experiment_log[x].group_indexes))
    # groups2=sorted(exp2.experiment_log.keys(),key=lambda x: len(exp2.experiment_log[x].group_indexes))
    groups1 = exp1.experiment_log.keys()
    groups2 = exp2.experiment_log.keys()

    if len(set(exp1.explored_set).difference(set(df1.index))) > 0 or len(
            set(exp2.explored_set).difference(set(df2.index))) > 0:
        print('Index mismatch')
        return graph,[]

    for cl1 in groups1:
        missing_group = True
        equivalent = False
        graph[cl1] = []

    clusterDict={}
    for cl2 in groups2:
        clusterDict[cl2]=set(df2.loc[exp2.experiment_log[cl2].group_indexes]['aminoAcid'])

    for cl1 in groups1:
        inS=[cl1]
        groups2 = set(exp2.experiment_log.keys())
        out = set()

        for sequence1 in exp1.experiment_log[cl1].group_indexes:
            if len(groups2)==0:
                break
            aminoString=df1['aminoAcid'][sequence1]
            variations=set(generateVariations(aminoString))
            for cl2 in groups2:
                common=variations.intersection(clusterDict[cl2])
                if len(common)>0:
                    found = True
                    graph[cl1] += [cl2]
                    out |= set([cl2])
                    inS = []
            groups2 = groups2.difference(out)

        missing+=inS


    return graph, missing


def matchClustersBack2(exp1, exp2):
    ##We first order clusters by size to better performance. Other optimizations (not implemented yet) will include:
    ##Sorting sequences by degree
    ## comparing samples of all groups to decide most likely matches and order accordingly
    # clusterOrder1,clusterOrder2=orderClusters(exp1,exp2)

    result = []
    missing = []
    graph = {}
    df1 = read_data(exp1.file_name)
    df2 = read_data(exp2.file_name)

    # groups1=sorted(exp1.experiment_log.keys(),key=lambda x: len(exp1.experiment_log[x].group_indexes))
    # groups2=sorted(exp2.experiment_log.keys(),key=lambda x: len(exp2.experiment_log[x].group_indexes))
    groups1 = exp1.experiment_log.keys()
    groups2 = exp2.experiment_log.keys()

    if len(set(exp1.explored_set).difference(set(df1.index))) > 0 or len(
            set(exp2.explored_set).difference(set(df2.index))) > 0:
        print('Index mismatch')
        return graph,[]

    for cl1 in groups1:
        missing_group = True
        equivalent = False
        graph[cl1] = []

    clusterDict={}
    for cl2 in groups2:
        clusterDict[cl2]=set(df2.loc[exp2.experiment_log[cl2].group_indexes]['aminoAcid'])

    for cl1 in groups1:
        inS=[cl1]
        groups2 = set(exp2.experiment_log.keys())
        out = set()

        for sequence1 in exp1.experiment_log[cl1].group_indexes:
            if len(groups2)==0:
                break
            aminoString=df1['aminoAcid'][sequence1]

            for i in range(len(aminoString)):
                for j in dnaValues:

                    if min_dist_amino[aminoString[i] + j] <= 1:
                        amino = aminoString[0:i] + j + aminoString[i + 1:]
                        for cl2 in groups2:
                            if amino in clusterDict[cl2]:
                                found = True
                                graph[cl1]+=[cl2]
                                out|=set([cl2])
                                inS=[]
                        groups2 = groups2.difference(out)

        missing+=inS


    return graph, missing


def matchClustersBack(exp1, exp2):
    ##We first order clusters by size to better performance. Other optimizations (not implemented yet) will include:
    ##Sorting sequences by degree
    ## comparing samples of all groups to decide most likely matches and order accordingly
    # clusterOrder1,clusterOrder2=orderClusters(exp1,exp2)

    result = []
    missing = []
    graph = {}
    df1 = read_data(exp1.file_name)
    df2 = read_data(exp2.file_name)

    # groups1=sorted(exp1.experiment_log.keys(),key=lambda x: len(exp1.experiment_log[x].group_indexes))
    # groups2=sorted(exp2.experiment_log.keys(),key=lambda x: len(exp2.experiment_log[x].group_indexes))
    groups1 = exp1.experiment_log.keys()
    groups2 = exp2.experiment_log.keys()

    if len(set(exp1.explored_set).difference(set(df1.index))) > 0 or len(
            set(exp2.explored_set).difference(set(df2.index))) > 0:
        print('Index mismatch')
        return result, graph

    for cl1 in groups1:
        missing_group = True
        equivalent = False
        graph[cl1] = []
        for cl2 in groups2:
            equivalent = isEquivalent(df1, df2, exp1.experiment_log[cl1].group_indexes,
                                      exp2.experiment_log[cl2].group_indexes, exact=False)
            if (equivalent):
                missing_group = False
                result += [(cl1, cl2)]
                graph[cl1] += [cl2]
        if missing_group:
            missing += [cl1]

    return result, graph, missing


def matchExperimentsBack(files):
    clusterMatches={}
    ## There is some redundant computation (because of symetry) but I leave them now for a sanity check
    for file1 in files:
        experiment1=load_data_object(file1)
        #experiment1.file_name = path + experiment1.file_name.split('/')[-1]
        clusterMatches[file1]={}
        clusterMatches[file1]['size']=len(experiment1.explored_set)##assumes not repeated sequences
        for file2 in files:
            if file1 != file2:
                experiment2=load_data_object(file2)
               # experiment2.file_name = path + experiment2.file_name.split('/')[-1]
                clusterMatches[file1][file2]={}

                graph,missing=matchClusters(experiment1,experiment2)
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
        #experiment1.file_name = path + experiment1.file_name.split('/')[-1]
        clusterMatches[file1]={}
        for y,file2 in enumerate(files):
            experiment2 = load_data_object(file2,dataPath=dataPath,experimentPath=experimentPath)
            #experiment2.file_name = path + experiment2.file_name.split('/')[-1]

            if x < y:
                clusterMatches[file1][file2] = {}
                graph,missing=matchClusters(experiment1,experiment2)
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
a1=read_data(experiment1.file_name)
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

    if len(sys.argv) != 2:
        sys.exit(0)

    config=read_config(sys.path[0]+'/'+sys.argv[1])
    experimentPath=config['experimentPath']
    datapath=config['dataPath']
    perc=(float)(config['percentage_min_group_size'])
    file=config['file'].lower()
    files = [f for f in listdir(experimentPath) if 'Experiment' in f]
    clusterMatches = matchExperiments(files,dataPath=datapath,experimentPath=experimentPath)
      #  r = json.dumps(clusterMatches)
      #  with open('matches.json', 'w') as f:
      #      json.dump(r, f)