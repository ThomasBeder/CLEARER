#####################################
# Code to extract features from NetNglyc server output
#####################################
from statistics import mean
import pandas as pd
import re


def netn(inputFile, path):
    # loop through file and compute count, mean etc
    id = ''
    iscore = []
    lock = False
    comment = False
    res = {}
    for line in inFile:
        if '##' in line and not comment:
            comment = True
            continue
        elif '##' in line and comment:
            comment = False
            continue
        
        if 'Name:' in line:
            # prevId = id
            id = line.split()[1]
            res[id] = {}
            continue
        if 'No sites predicted in this sequence' in line:
            res[id]['ngly_poscount'] = 0
            res[id]['ngly_posmean'] = 0
            res[id]['ogly_negcount'] = 0
            res[id]['ogly_negmean'] = 0
            continue
            
        if len(line) > 1 and not comment:
            # print(line)
            if id is not '' and id == line.split()[0]:
                lock = True
                score = line.split()[3]
                iscore.append(float(score))
                continue
            if id != line.split()[0] and lock:
                # get the count, mean and fraction of significant scores 
                sig_score = [x for x in iscore if float(x) >= 0.5]
                nonsig_score = [x for x in iscore if float(x) < 0.5]
                res[id]['ngly_poscount'] = len(sig_score)
                res[id]['ngly_posmean'] = mean(sig_score) if len(sig_score) > 0 else 0
                
                res[id]['ogly_negcount'] = len(nonsig_score)
                res[id]['ogly_negmean'] = mean(nonsig_score) if len(nonsig_score) > 0 else 0
        
                iscore = []
                lock = False
                
    # print(res)
    result = pd.DataFrame(res)
    result = result.transpose()
    result.index.names = ['realNames']
    # result.index = result.index.astype('int64')
    # convert ID to orgID
    # mapper = pd.read_csv('nameMapping/Dr.csv', index_col= 1)
    # out = pd.merge(mapper, result, left_index=True, right_index=True, how='inner')
    # url = inputFile.replace('outdata/netnglyc/Dr.fa','features/netnglyc/Dr_prot.csv')
    # out.to_csv(url, index=None) 
    if result.index[0] is '1':
        result.index = result.index.astype('int64')
        path = path.replace('\\', '/')
        print(path)
        pattern = '\/\w{2}\.'
        org = re.findall(pattern, path)
        org = org[0].strip('/').strip('.')
        map_url = 'nameMapping/' + org + '.csv'
        print(map_url)
        mapper = pd.read_csv(map_url, index_col=1)
        result = pd.merge(mapper, result, left_index=True, right_index=True, how='inner')
        url = path.replace('.fa', '.csv')
        result.to_csv(url, index=False)
    else:
        url = path.replace('.fa', '.csv')
        result.to_csv(url)

import os, glob
path = 'outdata_latest/netnglyc'
if os.path.exists(path):
    files = glob.glob(os.path.join(path, '*.fa'))
    for path in files:    
        # readin input file
        file = open(path, 'r')
        inFile = file.readlines()
        file.close()
        netn(inFile, path)
        print(path + ' completed successfully!')