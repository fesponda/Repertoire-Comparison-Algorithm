# -*- coding: utf-8 -*-
# Project: PublicClusters
# Filename: sample.py
# Authors: Fabio Calo (fcalodiz@itam.mx) and Joshua J. Daymunde (jdaymude@asu.edu).

"""
sample: Generates subsamples of given size from a data file.
"""

import pandas as pd
import os
import csv
import random as rand
import numpy as np
from datetime import date


def read_data(fname, datafmt):
    """
    Reads data from the given file as a DataFrame, removing all columns except
    the amino acid sequences and counts.
    
    :param fname: a string file name to read from (should be .tsv)
    :param datafmt: a string data format ('esponda2020' or 'emerson2017')
    :returns: a DataFrame containing amino acid sequences and their counts
    """
    
    # Read the data as a DataFrame, renaming columns according to the dataset.
    df = pd.read_csv(fname, sep = '\t')
    
    if datafmt == 'esponda2020':
        col_map = {'aminoAcid': 'amino_seq', 'count (templates/reads)': 'count'}
    elif datafmt == 'emerson2017':
        col_map = {'amino_acid': 'amino_seq', 'reads': 'count'}
    else:
        raise ValueError('unrecognized data format \"' + datafmt + '\"')
    
    df = df.rename(columns = col_map)
    
    # Group by amino acid sequence.
    df = df.groupby('amino_seq', as_index = False)['count'].sum()
    
    return df


def create_sample_dir(sampling_path, df, sample_type, i):
    """
    Creates a new directory and output file for the i-th sample, and copies the
    column headers to the new file.
    
    :param sampling_path: path used to create directories for the sample files
    :param df: dataframe that stores the file's data
    :param sample_type: 0 (repeated with replacement), 1 (unique with replacement) or
                        2 (without replacement)
    :param i: index of the current sample
    """
    
    try:
        # Creates the i-th sample's directory
        if not os.path.exists(sampling_path + str(date.today()) + "\\Sample" + str(i)):
            os.mkdir(sampling_path + str(date.today()) + "\\Sample" + str(i))
    
        # Creates the i-th sample file
        sample_file = open(sampling_path + str(date.today()) + "\\Sample" + str(i) + "\\Sample.tsv", 'a')
        
        # Copies the column headers to the new file
        with open(sampling_path + str(date.today()) + "\\Sample" + str(i) + "\\Sample.tsv", 'w') as out_file:
            tsv_writer = csv.writer(out_file, delimiter = '\t')
            
            if sample_type == 0 or sample_type == 1:
                tsv_writer.writerow([df.columns.values[0], df.columns.values[1]])
            else:
                tsv_writer.writerow([df.columns.values[0]])
                    
    except OSError:
        print('Could not make sample directory number ', i)


def roulette_selection(seq_dict, selected_dict):
    """
    Uses the roulette-selection method to randomly select a sequence. Each sequence has a
    range determined by its accumulated count, and the accumulated count of the sequence
    that precedes it. This method generates a random number between 1 and the total number
    of sequences, counting repetitions (i.e. the sum of all values in the count column), 
    and uses binary search to find the sequence corresponding to said number.
    
    :param seq_dict: dictionary that stores sequence indexes and each sequence's
                         accumulated count
    :param selected_dict: dictionary that keeps track of which sequences have already 
                          been selected
    :returns: the index of the sequence corresponding to the random number chosen.
    """
    
    # Total number of unique sequences.
    tot_uniq_seq = len(seq_dict)
    # Random number between 1 and the total number of sequences (counting repetitions).
    selected_index = rand.randint(1, seq_dict[tot_uniq_seq - 1])
    lo = 0
    hi = tot_uniq_seq
    seq_index = int(tot_uniq_seq / 2)
    seq_found = False
    
    # Uses binary search over the sequence indexes, until one is found whose range
    # includes the chosen random number.
    while not seq_found:
        if seq_index == 0 or (selected_index <= seq_dict[seq_index] and selected_index > seq_dict[max(0, seq_index - 1)]):
            seq_found = True
        elif selected_index > seq_dict[seq_index]:
            lo = seq_index
            seq_index += int((hi - lo + 1) / 2)
        else:
            hi = seq_index
            seq_index -= int((hi - lo + 1) / 2)
    
    return seq_index


def sample(fname, size, sample_type, datafmt, num_samples = 1):
    """
    Selects a subsample of given size, corresponding to the number of total sequences.
    The size may or may not count repetitions, depending on the sample_type parameter.
    
    :param fname: a string file name to read from (should be .tsv)
    :param size: number of sequences (repeated or unique) to select in each sample 
    :param sample_type: 0 (repeated with replacement), 1 (unique with replacement) or
                        2 (without replacement)
    :param datafmt: a string data format ('esponda2020' or 'emerson2017')
    :param num_samples: total number of samples desired
    """
    
    # Path used to create directories for the sample files.
    sampling_path = ""
    # Dictionary that stores sequence indexes and each sequence's accumulated count.
    seq_dict = {}
    # Dictionary that keeps track of which sequences have already been selected.
    selected_dict = {}
    # Dataframe that stores the file's data.
    df = read_data(fname, datafmt)
    
    # Populates seq_dict with the acumulated count of each sequence (sequence
    # indexes are used instead of the string names).
    for index, rows in df.iterrows():
        if index == 0:
            seq_dict[index] = rows.iloc[1]
        else:
            seq_dict[index] = seq_dict[index - 1] + rows.iloc[1]
    
    try:
        # Creates a new parent directory for the requested samples.
        if not os.path.exists(sampling_path + str(date.today())):
            os.mkdir(sampling_path + str(date.today()))
    except OSError:
        print('Could not make sample parent directory')
    
    for i in range(0, num_samples):
        create_sample_dir(sampling_path, df, sample_type, i)
        
        if sample_type == 0:
            sample_wr_repeated(seq_dict, selected_dict, size, i)
        elif sample_type == 1:
            sample_wr_unique(seq_dict, selected_dict, size, i)
        elif sample_type == 2:
            selected_list = sample_wor(df, seq_dict, size, i)
        else:
            raise ValueError('Illegal argument: sample_type = ' + sample_type +
                         ', expected 0 or 1')
        
        # Writes the selected data to an output file.
        with open(sampling_path + str(date.today()) + "\\Sample" + str(i) + "\\Sample.tsv", 'a') as out_file:
            if sample_type == 0 or sample_type == 1:
                for k in selected_dict.keys():
                    selected_data = [df.loc[k].iloc[0], selected_dict[k]]
                    tsv_writer = csv.writer(out_file, delimiter = '\t')
                    tsv_writer.writerow(selected_data)
            else:
                for k in selected_list:
                    selected_data = [df.loc[k].iloc[0]]
                    tsv_writer = csv.writer(out_file, delimiter = '\t')
                    tsv_writer.writerow(selected_data)
        
        selected_dict.clear()


def sample_wr_repeated(seq_dict, selected_dict, size, i):
    """
    Selects a subsample of given size (counting repetitions), with replacement.
    
    :param seq_dict: dictionary that stores sequence indexes and each sequence's
                         accumulated count
    :param selected_dict: dictionary that keeps track of which sequences have already 
                          been selected
    :param size: number of sequences to select in each sample (counting repetitions) 
    :param i: index of the current sample
    """
    
    # Selects the desired number of sequences into each sample file using the
    # roulette-wheel selection  method.
    for j in range(0, size):
        seq_index = roulette_selection(seq_dict, selected_dict)
        
        if seq_index not in selected_dict.keys():
            selected_dict[seq_index] = 1;
        else:
            selected_dict[seq_index] += 1;


def sample_wr_unique(seq_dict, selected_dict,  size, i):
    """
    Selects a subsample of given size (not counting repetitions).
    
    :param seq_dict: dictionary that stores sequence indexes and each sequence's
                         accumulated count
    :param selected_dict: dictionary that keeps track of which sequences have already 
                          been selected
    :param size: number of unique sequences to select in each sample
    :param i: index of the current sample
    """
    
    # Total number of unique sequences.
    tot_uniq_seq = len(seq_dict)
    
    # Selects the desired number of sequences into each sample file using the
    # roulette-wheel selection  method.
    while len(selected_dict) < min(size, tot_uniq_seq):
        seq_index = roulette_selection(seq_dict, selected_dict)
        
        if seq_index not in selected_dict.keys():
            selected_dict[seq_index] = 1;
        else:
            selected_dict[seq_index] += 1;


def sample_wor(df, seq_dict, size, i):
    """
    Selects a subsample of given size without replacement.
    
    :param df: dataframe that stores the file's data
    :param seq_dict: dictionary that stores sequence indexes and each sequence's
                         accumulated count
    :param selected_dict: dictionary that keeps track of which sequences have already 
                          been selected
    :param size: number of unique sequences to select in each sample
    :param i: index of the current sample
    """
    
    # Total number of unique sequences.
    tot_uniq_seq = len(seq_dict)
    # Stores each sequence's probability of being picked.
    P = []
    
    for k in seq_dict.keys():
        P.append(df.loc[k].iloc[1] / seq_dict[tot_uniq_seq - 1])
    
    return np.random.choice(range(0, tot_uniq_seq), size, replace = False, p = P)

