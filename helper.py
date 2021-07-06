# Project:  PublicClusters
# Filename: helper.py
# Authors:  Fernando Esponda (fernando.esponda@itam.mx) and Joshua J. Daymude
#           (jdaymude@asu.edu).

"""
helper: Helper functions for the RCA pipeline.
"""

import pandas as pd
import pickle
import os
import os.path as osp


def amino_acids():
    """
    :returns: a list of the 20 amino acids and the wildcard '*'
    """
    return ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', \
            'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']


def read_data(fname, datafmt):
    """
    Reads data from the given file as a DataFrame, removing all columns except
    the amino acid sequences and counts, grouping and sorting rows by count.

    :param fname: a string file name to read from (should be .tsv)
    :param datafmt: a string data format ('esponda2020' or 'emerson2017')
    :returns: a DataFrame containing amino acid sequences and their counts
    """
    # Read the data as a DataFrame, renaming columns according to the dataset.
    df = pd.read_csv(fname, sep='\t')
    if datafmt == 'esponda2020':
        col_map = {'aminoAcid': 'amino_seq', 'count (templates/reads)': 'count'}
    elif datafmt == 'emerson2017':
        col_map = {'amino_acid': 'amino_seq', 'reads': 'count'}
    else:
        raise ValueError('unrecognized data format \"' + datafmt + '\"')
    df = df.rename(columns=col_map)

    # Group by amino acid sequence and sort from most to least reads.
    df = df.groupby('amino_seq', as_index=False)['count'].sum()
    df = df.sort_values(by=['count'], ascending=False)

    return df


def dump_obj(fname, obj):
    """
    Writes a python object to file, creating any missing directories.

    :param fname: a string file name
    :param obj: a python object to write to file
    """
    os.makedirs(osp.split(fname)[0], exist_ok=True)
    with open(fname, 'wb') as f:
        pickle.dump(obj, f)


def load_obj(fname):
    """
    Reads a python object from file.

    :param fname: a string file name
    :returns: a python object
    """
    with open(fname, 'rb') as f:
        return pickle.load(f)