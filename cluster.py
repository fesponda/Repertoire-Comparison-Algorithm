from scipy.stats import spearmanr, stats
from scipy.spatial.distance import euclidean
from random import randint, shuffle,sample,random
from os import listdir
from os.path import isfile, join
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import datetime
from time import time
import pickle
import sys
#% matplotlib inline
#path='/Users/fesponda/research/Immunology/Data/'
experimentPath=''
#path='~/Immunology/Data/'
#path='/home/cesponda/Immunology/Data/'

mouse='Spleen_1'
min_group_size=5
perc=0.01
#mouse='S_7_AAA_2_norm_plus'

def create_codon_dict():
    dna_slc={}
    dna_slc['ATT']='I'
    dna_slc['ATC']='I'
    dna_slc['ATA']='I'
    dna_slc['CTT']='L'
    dna_slc['CTC']='L'
    dna_slc['CTA']='L'
    dna_slc['CTG']='L'
    dna_slc['TTA']='L'
    dna_slc['TTG']='L'
    dna_slc['GTT']='V'
    dna_slc['GTC']='V'
    dna_slc['GTA']='V'
    dna_slc['GTG']='V'
    dna_slc['TTT']='F'
    dna_slc['TTC']='F'
    dna_slc['ATG']='M'
    dna_slc['TGT']='C'
    dna_slc['TGC']='C'
    dna_slc['GCT']='A'
    dna_slc['GCC']='A'
    dna_slc['GCA']='A'
    dna_slc['GCG']='A'
    dna_slc['GGT']='G'
    dna_slc['GGC']='G'
    dna_slc['GGA']='G'
    dna_slc['GGG']='G'
    dna_slc['CCT']='P'
    dna_slc['CCC']='P'
    dna_slc['CCA']='P'
    dna_slc['CCG']='P'
    dna_slc['ACT']='T'
    dna_slc['ACC']='T'
    dna_slc['ACA']='T'
    dna_slc['ACG']='T'
    dna_slc['TCT']='S'
    dna_slc['TCC']='S'
    dna_slc['TCA']='S'
    dna_slc['TCG']='S'
    dna_slc['AGT']='S'
    dna_slc['AGC']='S'
    dna_slc['TAT']='Y'
    dna_slc['TAC']='Y'
    dna_slc['TGG']='W'
    dna_slc['CAA']='Q'
    dna_slc['CAG']='Q'
    dna_slc['AAT']='N'
    dna_slc['AAC']='N'
    dna_slc['CAT']='H'
    dna_slc['CAC']='H'
    dna_slc['GAA']='E'
    dna_slc['GAG']='E'
    dna_slc['GAT']='D'
    dna_slc['GAC']='D'
    dna_slc['AAA']='K'
    dna_slc['AAG']='K'
    dna_slc['CGT']='R'
    dna_slc['CGC']='R'
    dna_slc['CGA']='R'
    dna_slc['CGG']='R'
    dna_slc['AGA']='R'
    dna_slc['AGG']='R'
    dna_slc['TAA']='Stop'
    dna_slc['TAG']='Stop'
    dna_slc['TGA']='Stop'
    return dna_slc

# mindistance between amino acids
def min_dist_between_aminoacids(dna_slc):
    # revese codon dict first
    amino2codon={}
    for k in dna_slc.values():
        amino2codon[k]=[]
    for k, v in dna_slc.items():
        amino2codon[v].append(k)
    min_dist_amino={}
    for ami1 in amino2codon.keys():
        for ami2 in amino2codon.keys():
            #if ami1==ami2:
            #    continue
            min_dist=3
            for codon1 in amino2codon[ami1]:
                for codon2 in amino2codon[ami2]:
                    dist=0
                    for i in range(3):
                        if codon1[i]!=codon2[i]:
                            dist+=1
                    if dist < min_dist:
                        min_dist=dist
            min_dist_amino[ami1+ami2]=min_dist
            min_dist_amino['*'+ami1]=0
            min_dist_amino[ami1+'*']=0
            min_dist_amino['**']=0
    return min_dist_amino


dna_slc=create_codon_dict()
min_dist_amino=min_dist_between_aminoacids(dna_slc)
dnaValues=list(set(dna_slc.values()))
dnaValues+=['*']
dnaValues.remove('Stop')

class cluster_obj:

    def __init__(self, description, group_indexes):
        self.description = description
        self.group_indexes = group_indexes
        self.equivalent_groups = {}

    def compare_with_other_cluster(self, other, df, df2):
        ###aminoAcids
        in_self = df.loc[self.group_indexes]['aminoAcid'].unique()
        in_other = df2.loc[other.group_indexes]['aminoAcid'].unique()
        intersect = set(in_self).intersection(set(in_other))
        union = set(in_self).union(set(in_other))
        # return 100*len(intersect)/len(union)

        ##nucleotides
        in_self = df.loc[self.group_indexes]['nucleotide'].unique()
        in_other = df2.loc[other.group_indexes]['nucleotide'].unique()
        intersect = set(in_self).intersection(set(in_other))
        union = set(in_self).union(set(in_other))
        return 100 * len(intersect) / len(union)

        # in_self=df['aminoAcid'].unique()
        # in_other=df2['aminoAcid'].unique()
        # intersect=set(in_self).intersection(set(in_other))
        # union=set(in_self).union(set(in_other))
        # print(len(intersect)/len(union))

    def set_within_metric(self, number_of_amino_in, number_of_amino_out, number_of_amino_tot, number_of_amino_leaked):
        self.number_of_amino_in = number_of_amino_in
        self.number_of_amino_out = number_of_amino_out
        self.number_of_amino_tot = number_of_amino_tot
        self.number_of_amino_leaked = number_of_amino_leaked

    def set_equivalent_group(self, description, other_file_name, cluster_file):
        # mandar llamar
        grupos, df2, tot_templates_new, tot_templates_clus, templates_new, templates_clus = \
            extract_equivalent_cluster(other_file_name, cluster_file, self.group_indexes)
        exp = cluster_obj(description, grupos[0])

        exp.set_within_metric(*evaluate_within(df2, grupos[0])) #comented for other data
        self.equivalent_groups[other_file_name] = {'cluster_obj': exp, \
                                                   'sizes': (tot_templates_new, tot_templates_clus, templates_new,
                                                             templates_clus),
                                                   'relative_sizes': (1.0 * templates_new / tot_templates_new,
                                                                      1.0 * templates_clus / tot_templates_clus)}


class experiment:
    def __init__(self, description, file_name, source=True):
        self.description = description
        self.file_name = file_name
        self.source = source
        self.exp_number = 0
        self.experiment_log = {}
        self.explored_set = []
        self.consolidated={}
        self.consolidated_clusters={}

    def add_experiment(self, group_indexes, df=''):
        self.explored_set += group_indexes
        self.exp_number += 1
        exp = cluster_obj(self.description, group_indexes)
        if len(df) == 0:
            df = read_data(self.file_name)
        exp.set_within_metric(*evaluate_within(df, group_indexes))  # method outside class for now
        self.experiment_log[self.exp_number] = exp

    def set_equivalent_groups(self, description, other_file_name, groups=[]):
        for g in groups:
            if g in self.experiment_log.keys():
                self.experiment_log[g].set_equivalent_group(description, other_file_name, self.file_name)

    def compute_group_intersections(self):
        # To avoid possible futere mistakes Ill keep redundant info and do all against all (except self)
        inter = {}
        for other_individuals in list(self.experiment_log[1].equivalent_groups.keys()):
            inter[other_individuals] = {}
            for cluster in range(1, self.exp_number + 1):
                inter[other_individuals][cluster] = []
                c1 = self.experiment_log[cluster].equivalent_groups[other_individuals][
                    'cluster_obj'].group_indexes
                for cluster2 in range(1, self.exp_number + 1):
                    if cluster != cluster2:
                        c2 = self.experiment_log[cluster2].equivalent_groups[other_individuals][
                            'cluster_obj'].group_indexes
                        if len(set(c1).intersection(set(c2))) > 0:
                            inter[other_individuals][cluster].append(cluster2)

        return inter

    def traverse_consolidate(self,consolidate, visitado, llave, to_consolidate):
        visitado[llave] = True
        to_consolidate.append(llave)
        for val in consolidate[llave]:
            if not visitado[val]:
                self.traverse_consolidate(consolidate, visitado, val, to_consolidate)

    '''Compute intersection between groups within each data set and consolidate. It sufficies that there is an intersection
    in one dataset for the groups in all datasets to be merged'''

    def consolidate_groups(self):

        # for each dataset it computes which clusters have common sequences
        inter = self.compute_group_intersections()
        # If two clusters have common sequences we set them up for consolidation
        consolidate = {}
        for group in list(inter[list(inter.keys())[0]].keys()):
            consolidate[group] = [group]
            for key in inter.keys():
                consolidate[group] += inter[key][group]
            consolidate[group] = set(consolidate[group])
        # We compute the clusters to consolidate. The previous step creates a graphlike structure
        # in which it establishes, for example, that cluster 1 should be consolidated with cluster 2 and 2 with 3
        # We traverse it as a graph to consolidate all clusters in the path
        visitado = {}
        for llave in consolidate:
            visitado[llave] = False
        for llave in visitado:
            if not visitado[llave]:
                to_consolidate = []
                # Extracts one connected component in the graph. Clusters to consolidate
                self.traverse_consolidate(consolidate, visitado, llave, to_consolidate)
                r = set()
                for val in to_consolidate:

                    r |= consolidate[val]
                    if val != llave:
                        consolidate.pop(val, None)
                consolidate[llave] |= r
        self.consolidated=consolidate


    def consolidate(self):
        '''Compute the actual intersection. Just distinct sequences'''
        inter = {}
        ###this is for the original group. Perhaps change class strucuture to avoid
        inter[self.file_name] = {}
        count = 0
        for llave in self.consolidated:
            c1 = []
            for cluster in self.consolidated[llave]:
                c1 += self.experiment_log[cluster].group_indexes
            inter[self.file_name][llave] = set(c1)
            count += len(set(c1))

        for other_individuals in list(self.experiment_log[1].equivalent_groups.keys()):
            inter[other_individuals] = {}
            for llave in self.consolidated:
                c1 = []
                for cluster in self.consolidated[llave]:
                    c1 += self.experiment_log[cluster].equivalent_groups[other_individuals]['cluster_obj'].group_indexes
                inter[other_individuals][llave] = set(c1)

        self.consolidated_clusters=inter  #pd.DataFrame.from_dict(inter, orient='columns')

def dist_amino(s1,s2,max_dist): #assiming for now min_dist_amino as global so as not to carry around
    if len(s2)<len(s1):
        temp=s1
        s1=s2
        s2=temp
    dist=len(s2)-len(s1)
    for i in range(len(s1)):
        dist+=min_dist_amino[s1[i]+s2[i]]
        if dist > max_dist:
            break
    return dist


# Does not calculate the exact edit distance, it ends if it finds a transformation within min_dist
def edit_distance(s1,s2,max_dist,cont=0):
    if len(s1)==0 or len(s2)==0:
        return cont+len(s2)+len(s1)
    if cont>max_dist:
        return cont
    if s1[0]==s2[0]:
        return edit_distance(s1[1:],s2[1:],max_dist,cont)
    subs=edit_distance(s1[1:],s2[1:],max_dist,cont + 1)
    if subs <= max_dist:
        return subs
    remove=edit_distance(s1[1:],s2[0:],max_dist,cont + 1)
    if remove <= max_dist:
        return remove
    insert=edit_distance(s1[0:],s2[1:],max_dist,cont + 1)
    return min([subs,remove,insert])


def read_data(file1,sep=' '):
    if file1.split('.')[-1]=='csv':
        sep=','
    elif file1.split('.')[-1]=='tsv':
        sep='\t'
    df=pd.read_csv(file1,sep=sep)
    df.dropna(inplace=True,subset=['aminoAcid'])
    df=df.groupby('aminoAcid', as_index=False)['count (templates/reads)'].sum()
    df=df.sort_values(by=['count (templates/reads)'], ascending=False)
    return df


def save_clusters(df,grupos,mouse,path='',date='',other_description=''):
    names={}
    for key in grupos.keys():
        name=path+other_description+mouse+'group'+str(key)+date+'.csv'
        df.loc[grupos[key]].to_csv(name)
        names[key]=name
    return names


##Implements DBSCAN and returns a dictionary with the groups found. It contains the data frame indeces
def clusterOLD(df,num_groups=1,min_group_size=300):
    out_group=list(df.index.values)
    in_group=[]
    grupos={}
    max_dist=1 # for edit distance
    k=0
    while k< num_groups and len(out_group)>0:
        in_group.append(out_group.pop(0))
        cluster=[]
        while len(in_group) >0:
            sequence=in_group.pop(0)
            cluster.append(sequence)
            current_matches=[]
            for i in out_group:
                corr=dist_amino(df['aminoAcid'][sequence],df['aminoAcid'][i],max_dist)
                if(corr<=max_dist):
                    in_group.append(i)
                    current_matches.append(i)
            for sequence in current_matches:
                try:
                    out_group.remove(sequence)
                except:
                    print('not in')
            #if (len(cluster) % 50) ==0:
             #   print(len(out_group),len(in_group),len(cluster))
        if len(cluster)>= min_group_size:
            grupos[k]=cluster.copy()
            k+=1
        else:
            print('k=',k,' remaining=',len(out_group))

    return grupos

###New version

'''
# Only efficient when distance is 1 or 0
'''


def checkVariationsIndex(aminoString, cluster2, distanceMetric='Hamming'):
    current_matches = set()
    if distanceMetric=='Hamming':
        for i in range(len(aminoString)):
            for j in dnaValues:
                amino = aminoString[0:i] + j + aminoString[i + 1:]
                if amino in cluster2.keys():  ##substitutions
                    current_matches.add(cluster2[amino])
                    cluster2.pop(amino, None)
    elif distanceMetric=='Edit':

        for i in range(len(aminoString)):
            for j in dnaValues:
                if True or min_dist_amino[aminoString[i] + j] <= 1: #True so that all aminoacids are 1 away
                    amino = aminoString[0:i] + j + aminoString[i + 1:]
                    if amino in cluster2.keys():  ##substitutions
                        current_matches.add(cluster2[amino])
                        cluster2.pop(amino, None)
                amino=aminoString[0:i] + j + aminoString[i:]
                if amino in cluster2.keys():  # insertion
                    current_matches.add(cluster2[amino])
                    cluster2.pop(amino, None)
            amino=aminoString[0:i] + aminoString[i + 1:]
            if amino in cluster2.keys():  # deletion
                current_matches.add(cluster2[amino])
                cluster2.pop(amino, None)

        for j in dnaValues:  # insertion in the last position
            amino=aminoString + j
            if amino in cluster2.keys():
                current_matches.add(cluster2[amino])
                cluster2.pop(amino, None)

    return current_matches



'''
# Only efficient when distance is 1 or 0
'''
#Only substitutions
def checkVariationsIndexOld(aminoString, cluster2):
    current_matches = set()
    for i in range(len(aminoString)):
        for j in dnaValues:
            if min_dist_amino[aminoString[i] + j] <= 1:
                amino = aminoString[0:i] + j + aminoString[i + 1:]
                if amino in cluster2.keys():
                    current_matches.add(cluster2[amino])
                    cluster2.pop(amino, None)
    for j in dnaValues:  ##one long
        if aminoString + j in cluster2.keys():
            current_matches.add(cluster2[aminoString + j])
            cluster2.pop(aminoString + j, None)
    if aminoString[0:-1] in cluster2.keys():  # one short
        current_matches.add(cluster2[aminoString[0:-1]])
        cluster2.pop(aminoString[0:-1], None)
    return current_matches


##new for cluster


##Implements DBSCAN and returns a dictionary with the groups found. It contains the data frame indeces
def cluster(df,num_groups=1,min_group_size=30):
    out_group=set(df.index.values)
    aminoSet=dict(zip(df.aminoAcid,df.index))
    in_group=set()
    grupos={}
    max_dist=1 # for edit distance
    k=0
    while k< num_groups and len(out_group)>0:
        in_group.add(out_group.pop())
        cluster=[]
        while len(in_group) >0:
            sequence=in_group.pop()
            cluster.append(sequence)
            current_matches=checkVariationsIndex(df['aminoAcid'][sequence],aminoSet)
            in_group|=(current_matches-set([sequence]))
           # print('before',len(out_group),len(current_matches))
            out_group-=current_matches
           # print(len(out_group),len(current_matches))
        #    if (len(cluster) % 50) ==0:
        #        print(len(out_group),len(in_group),len(current_matches),len(cluster))


        if len(cluster)>= min_group_size:

            grupos[k]=cluster.copy()
            k+=1
    return grupos

###Calculate the  aminoacid 'leakage'
##
def evaluate_within(df, group_indexes):
    if 'vMaxResolved' not in df.columns:
        return 1, 1, 1, 1
    indexes_not_in_group = set(df.index.values).difference(set(group_indexes))

    all_data = df['aminoAcid'].unique()
    in_group = df.loc[group_indexes]['aminoAcid'].unique()
    ##filter v segments
    v_segments = df.loc[group_indexes]['vMaxResolved'].unique()
    ###filter j segments
    j_segments = df.loc[group_indexes]['jMaxResolved'].unique()
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
    if (len(group_indexes) > 0):
        print(1.0 * number_of_amino_leaked / number_of_amino_in)
    print(number_of_amino_in, number_of_amino_out, number_of_amino_tot, number_of_amino_leaked)
    return number_of_amino_in, number_of_amino_out, number_of_amino_tot, number_of_amino_leaked



## This is for extracting equivalent groups.
def extract_equivalent_cluster(new_data, cluster_file, group_indexes, field='aminoAcid'):
    df = read_data(new_data)[['nucleotide', 'aminoAcid', 'count (templates/reads)']]
    df2 = read_data(cluster_file)[['nucleotide', 'aminoAcid', 'count (templates/reads)']]
    tot_templates_new = df['count (templates/reads)'].sum()
    tot_templates_clus = df2['count (templates/reads)'].sum()

    df2 = df2.loc[group_indexes]
    df2.index = [str(i) + 'a' for i in range(len(df2))]
    df = pd.concat([df2, df])
    grupos2 = cluster(df[[field]], num_groups=1)
    grupos = {}
    grupos[0] = [i for i in grupos2[0] if str(i)[-1] != 'a']

    templates_new = df.loc[grupos[0]]['count (templates/reads)'].sum()
    templates_clus = df2['count (templates/reads)'].sum()

    return grupos, df[len(df2):], tot_templates_new, tot_templates_clus, templates_new, templates_clus


def save_rest_of_data(df, nombre, grupos, llaves=[0]):
    l_llaves = ''
    indexes = []
    for i in llaves:
        indexes += grupos[i]
        l_llaves = l_llaves + '-' + str(i)
    nombre = nombre + 'Rest' + l_llaves + '.csv'  # +l_llaves+datetime.datetime.fromtimestamp(time()).strftime('%Y-%m-%d')+'.csv'
    indexes_not_in_group = set(df.index.values).difference(set(indexes))
    df.loc[df.index.isin(indexes_not_in_group)].to_csv(nombre, index=False)

def save_data_object(name,data_object):
    with open (name,'wb') as fp:
        pickle.dump(data_object,fp)

def load_data_object(name,dataPath='',experimentPath=''):

    name=experimentPath+name.split('/')[-1]
    with open(name,'rb') as fp:
        data_object=pickle.load(fp)
    data_object.file_name = dataPath + data_object.file_name.split('/')[-1]
    return data_object

def read_config(file):
    config={}
    with open(file,'r') as f:
        for line in f:
            l=line.replace(" ",'').strip().split(':')
            config[l[0]]=l[1]
    if len(config)==0:
        print('No config file')
        sys.exit(0)

    return config

def select_data_subset(df,explored_indexes):
    indexes_not_explored=set(df.index.values).difference(set(explored_indexes))
    return df.loc[indexes_not_explored]

def run_one_group_all_files(experiment2,num):
    #l = [f for f in listdir(path) if 'S_7' in f]
    l = [f for f in listdir(path) if experiment2.file_name.split('/')[-1] not in f and '.tsv' in f]
    for file in l:
        experiment2.set_equivalent_groups('Using nucleotide sequences, distance=1', path + file, [num])
    name = 'other-data-group' + str(num)+ '-' + datetime.datetime.fromtimestamp(time()).strftime('%Y-%m-%d-%H')
    save_data_object(experimentPath+name, experiment2)
        #print(experiment2.experiment_log[num].equivalent_groups)


def run_all_files_inpath(path):


    subjects = [f for f in listdir(path) if isfile(join(path, f)) and f[0] != '.' and ".tsv" in f]
    for mouse in subjects:
        file = path + mouse
        experiment_1 = experiment('Using amino acid sequences, distance=1', file)
        experiment_1.experiment_number = '1'
        run_all_groups(experiment_1)




def run_all_groups_new_files(experiment2,start,finish):
    """

    :param experiment2: Experiment object
    :return:
    """
    l = [f for f in listdir(path) if 'Synthetic_mouse' in f]

    for file in l:
        experiment2.set_equivalent_groups('Using aminoacid sequences, distance=1', path + file, range(start,finish+1))
    name = 'synthetic-data-groups' + '3_synthmice'+ '-' + datetime.datetime.fromtimestamp(time()).strftime('%Y-%m-%d-%H')
    save_data_object(experimentPath+name, experiment2)
        #print(experiment2.experiment_log[num].equivalent_groups)

def run_all_groups(experiment1):
    cluster_num = 1
    df = read_data(experiment1.file_name)
    min_size = (int)(perc*len(df))
    if  min_size < min_group_size:
        min_size=min_group_size

    while cluster_num < 5000: ##maximim number of clusters to find

        #df = read_data(experiment1.file_name)

        df = select_data_subset(df, experiment1.explored_set)
        #df = df.sort_values(by=['count (templates/reads)'], ascending=False)

        if len(df) < min_size: #should be min size
            print('done')
            break
        # print('menos')
        grupos = cluster(df[['aminoAcid']],min_group_size=min_size)
        if len(grupos) == 0:
            break
        # print('cero')
        experiment1.add_experiment(grupos[0],df)#//avoid reading I think
        # print('uno')
        ''' 
        experiment1.set_equivalent_groups('Using nucleotide sequences, distance=1' + experiment_number,
                                          path + 'Spleen_2.tsv', [experiment1.exp_number])
        # print('dos')
        experiment1.set_equivalent_groups('Using nucleotide sequences, distance=1' + experiment_number,
                                          path + 'Spleen_3.tsv', [experiment1.exp_number])
        '''
#        name = mouse + 'Experiment' + experiment_number + '-' + datetime.datetime.fromtimestamp(time()).strftime(
#            '%Y-%m-%d-%H')
#       save_data_object(name, experiment1)

        cluster_num += 1
        #print('Iteration = ', cluster_num, ' Clusters processed = ', experiment1.exp_number)
    name = experiment1.file_name.split('/')[-1].split('.')[0] + 'Experiment' + experiment1.experiment_number + '-' + datetime.datetime.fromtimestamp(time()).strftime(
            '%Y-%m-%d-%H')
    # save_data_object(path + 'experiments/' + name, experiment1)
    save_data_object(experimentPath+name, experiment1)
if __name__ == "__main__":

    #print ('Number of arguments:', len(sys.argv), 'arguments.')
    #print('Argument List:', str(sys.argv))
    #print(path)
    #sys.exit(0)

    ###IF old experimet, read
    '''this lines where for when we ran partial experiments
    if len(sys.argv) > 1:
        experiment1 = load_data_object(str(sys.argv[1]))
        if len(sys.argv) == 2:
            run_all_groups(experiment1)
        elif len(sys.argv) == 3:
            run_one_group_all_files(experiment1,int(sys.argv[2]))
        else:
            run_all_groups_new_files(experiment1,int(sys.argv[2]),int(sys.argv[3]))
    else: # if new,create new
        experiment1=experiment('Using amino acid sequences, distance=1', file)
        run_all_groups(experiment1)
    '''
    ##Determine if we are runing one file or a whole directory
    if len(sys.argv) != 2:
        sys.exit(0)

    config=read_config(sys.argv[1])
    experimentPath=config['experimentPath']
    datapath=config['dataPath']
    perc=(float)(config['percentage_min_group_size'])
    file=config['file']

    if file != 'all': #one file

        experiment1 = experiment('Using amino acid sequences, distance=1', datapath+file)
        experiment1.experiment_number = '1'
        run_all_groups(experiment1)
    else: #all files in path
        run_all_files_inpath(datapath)
