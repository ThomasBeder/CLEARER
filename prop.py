#####################################
# Code to extract features from Prop server output
#####################################
from statistics import mean
import pandas as pd
import re

def prop(inputFile, path):
    # loop through file and compute count, mean etc
    id = ''
    iscore = []
    lock = False
    res = {}
    pattern = '\s+\d+\s\w{5,18}'
    # line = '  467 ENSDARG00000001057 '
    # id = re.search(pattern, line, re.I)
    # print(id.group().split()[1])
    for line in inputFile:
        getid = re.search(pattern, line, re.I)
        # print(line)
        if getid is not None:
            id = getid.group().split()[1]
            # print(id)
            res[id] = {}
            continue
        
        if len(line) > 1:
            # print(line)
            if id is not '' and id in line:
                lock = True
                score = line.split( )[3]
                # print('The score is %s' %score)
                iscore.append(float(score))
                continue
            if id not in line and lock:
                # get the count, mean and fraction of significant scores 
                sig_score = [x for x in iscore if float(x) >= 0.5]
                # nonsig_score = [x for x in iscore if float(x) < 0.5]
                res[id]['ngly_poscount'] = len(sig_score)
                res[id]['ngly_posmean'] = mean(sig_score) if len(sig_score) > 0 else 0
                # res[id]['ogly_negcount'] = len(nonsig_score)
                # res[id]['ogly_negmean'] = mean(nonsig_score) if len(nonsig_score) > 0 else 0
                iscore = []
                lock = False

    result = pd.DataFrame(res)
    result = result.transpose()
    result.index.names = ['ID']
    # convert ID to orgID
    if result.index[0] is '1':
        # print(result.index[0])
        result.index = result.index.astype('int64')
        path = path.replace('\\', '/')
        # print(path)
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
path = 'outdata_latest/prop'
if os.path.exists(path):
    files = glob.glob(os.path.join(path, '*.fa'))
    for path in files:    
        # readin input file
        file = open(path, 'r')
        inFile = file.readlines()
        file.close()
        # print(inFile[8])
        prop(inFile, path)
        print(path + ' completed successfully!')