#####################################
# Code to extract features from NetCglyc server output
#####################################
from statistics import mean
import pandas as pd
import re
        
def netc(inputFile, path):
    # loop through file and compute count, mean etc
    id = ''
    iscore = []
    lock = False
    res = {}
    for line in inFile:
        if 'C-manno' in line and not lock:
            # prevId = id
            id = line.split()[0]
            # print(id)
            res[id] = {}
            lock = True
            continue
    
        if id == line.split()[0] and lock:
            score = line.split()[5]
            iscore.append(float(score))
            continue
        if id != line.split()[0] and lock:
            # get the count, mean and fraction of significant scores
            sig_iscore = [x for x in iscore if float(x) >= 0.5]
            res[id]['cgly_sigcount'] = len(sig_iscore)
            res[id]['cgly_count'] = len(iscore)
            # print(sig_gscore)
            res[id]['cgly_mean'] = mean(sig_iscore) if len(sig_iscore) > 0 else 0
            if len(sig_iscore) > 0 and len(iscore) > 0:
                res[id]['cgly_sigfrac'] = round(len(sig_iscore) / len(iscore), 3)
            else:
                res[id]['cgly_sigfrac'] = 0
            iscore = []
            lock = False
    # print(res)
    result = pd.DataFrame(res)
    result = result.transpose()
    result.index.names = ['realNames']
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
        url = path.replace('.fa', '.csv')
        result.to_csv(url, index=False)
    else:
        url = path.replace('.fa', '.csv')
        result.to_csv(url)

# path = 'outdata_latest/netcglyc/Dr_prot.fa'

# print(org)
# file = open(path, 'r')
# inFile = file.readlines()
# file.close()
# netc(inFile, path)
import os, glob
path = 'outdata_latest/netcglyc'
if os.path.exists(path):
    files = glob.glob(os.path.join(path, '*.fa'))
    for path in files:    
        # readin input file
        file = open(path, 'r')
        inFile = file.readlines()
        file.close()
        netc(inFile, path)
        print(path + ' completed successfully!')
    


