#####################################
# Code to extract features from NetOglyc server output
#####################################
from statistics import mean
import pandas as pd
import re


def neto(inputFile, path):
    id = ''
    gscore = []
    iscore = []
    lock = False
    res = {}
    for line in inFile:
        if 'Name:' in line:
            # prevId = id
            id = line.split()[1]
            res[id] = {}
            continue
        if len(line) > 1:
            if id is not '' and id == line.split()[0]:
                lock = True
                score = line.split()
                gscore.append(float(score[3]))
                iscore.append(float(score[4]))
                continue
            if id != line.split()[0] and lock:
                # get the count, mean and fraction of significant scores 
                sig_gscore = [x for x in gscore if float(x) >= 0.5]
                sig_iscore = [x for x in iscore if float(x) >= 0.5]
                res[id]['ogly_gcount'] = len(sig_gscore)
                # print(sig_gscore)
                res[id]['ogly_gmean'] = mean(sig_gscore) if len(sig_gscore) > 0 else 0
                res[id]['ogly_gfrac'] = round(len(sig_gscore)/len(gscore),3)
                
                res[id]['ogly_icount'] = len(sig_iscore)
                res[id]['ogly_imean'] = np.mean(sig_iscore) if len(sig_iscore) > 0 else 0
                res[id]['ogly_ifrac'] = round(len(sig_iscore)/len(iscore), 3)
        
                gscore = []
                iscore = []
                lock = False
    result = pd.DataFrame(res)
    result = result.transpose()
    result.index.names = ['realNames']
    # convert ID to orgID
    if result.index[0] is '1':
        # print(result.index[0])
        result.index = result.index.astype('int64')
        path = path.replace('\\', '/')
        # print(path)
        pattern = '\/\w{2}\.'
        org = re.findall(pattern, path)
        org = org[0].strip('/').strip('.')
        map_url = 'nameMapping/' + org + '.csv' # remove '_renamed' and '_additional'
        print(map_url)
        mapper = pd.read_csv(map_url, index_col=1) # confirm if index_col is 0 0r 1
        result = pd.merge(mapper, result, left_index=True, right_index=True, how='inner')
        url = path.replace('.fa', '.csv')
        result.to_csv(url, index=False)
    else:
        url = path.replace('.fa', '.csv')
        result.to_csv(url)


path = 'outdata_latest/raw_data/netoglyc/Ag.fa'
with open(path,'r') as inFile:    # line used for single file
    neto(inFile, path)

# code section used for multiple files
import os, glob
path = 'outdata_latest/raw_data/netoglyc'
if os.path.exists(path):
    files = glob.glob(os.path.join(path, '*.fa'))
    for path in files:    
        # readin input file
        file = open(path, 'r')
        inFile = file.readlines()
        file.close()
        neto(inFile, path)
        print(path + ' completed successfully!')