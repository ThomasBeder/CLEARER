#####################################
# Code to extract features from TMHMM server output
#####################################
# This code might not be necessary since we can import the file to excel and use tab as delimiter
# Three fields (ExpAA, First60, PredHel) are extracted from the input file
from statistics import mean
import pandas as pd
import re

def tmhmm(inputFile, path):
    # loop through file and compute count, mean etc
    id = ''
    lock = False
    res = {}
    for line in inputFile:
        if 'Length:' in line:
            # prevId = id
            id = line.split()[1]
            res[id] = {}
            res[id]['tm_len'] = line.split()[3]
            lock = True
            continue
        if lock and 'Number of predicted TMHs' in line: 
            res[id]['TMHs'] = line.split(':')[1].strip()
            continue
        if lock and 'Exp number of AAs in TMHs' in line: 
            res[id]['ExpAA'] = line.split(':')[1].strip()
            continue
        if lock and 'Exp number, first 60 AAs' in line: 
            res[id]['First60AA'] = line.split(':')[1].strip()
            continue
        if lock and 'Total prob of N-in' in line: 
            res[id]['prob_N_in'] = line.split(':')[1].strip()
            continue
        lock = False
    
    result = pd.DataFrame(res)
    result = result.transpose()
    result.index.names = ['realNames']
    # print(result.head())
    # convert ID to orgID
    if result.index[0] is '1':
        result.index = result.index.astype('int64')
        path = path.replace('\\', '/')
        pattern = '\/\w{2}\_'
        org = re.findall(pattern, path)
        org = org[0].strip('/').strip('_')
        map_url = 'nameMapping/' + org + '.csv'
        print(map_url)
        mapper = pd.read_csv(map_url, index_col=1)
        result = pd.merge(mapper, result, left_index=True, right_index=True, how='inner')
        url = path.replace('.fa.out', '.csv')
        result.to_csv(url, index=False)
    else:
        url = path.replace('.fa.out', '.csv')
        result.to_csv(url)

import os, glob
path = 'outdata_latest/tmhmm'
if os.path.exists(path):
    files = glob.glob(os.path.join(path, '*.fa.out'))
    for path in files:    
        # readin input file
        file = open(path, 'r')
        inFile = file.readlines()
        file.close()
        # print(inFile[8])
        tmhmm(inFile, path)
        print(path + ' completed successfully!')