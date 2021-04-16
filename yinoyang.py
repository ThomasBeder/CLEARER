#####################################
# Code to extract features from Yin-Yang  server output
#####################################
from statistics import mean
import pandas as pd
import re

def yinoyang(inputFile, path, TrueID = True):
    # loop through file and compute count, mean etc
    id = ''
    iscore, count = [], []
    lock = False
    res = {}
    pattern = '\+'
    
    for line in inputFile:
        # getid = re.search(pattern, line, re.I)
        # if getid is not None:
        #     id = getid.group().split()[1]
        #     res[id] = {}
        #     continue
        if 'Name:' in line:
            # prevId = id
            id = line.split()[1]
            res[id] = {}
            continue
        
        if len(line) > 1:
            # print(line)
            if id is not '' and id in line:
                lock = True
                getsig = re.findall(pattern, line, re.I)
                sig = len(getsig)
                score = line.split( )[4]
                count.append(sig)
                iscore.append(float(score))
                continue
            # if len(iscore) < 1 and id not in line and lock: # assign zero score to genes without prediction
            #     # get the count, mean and fraction of significant scores 
            #     # sig_score = [x for x in iscore if float(x) >= 0.5]
            #     # nonsig_score = [x for x in iscore if float(x) < 0.5]
            #     res[id]['oyang_count'] = 0
            #     res[id]['oyang_max'] = 0
            #     res[id]['oyang_mostsig'] = 0
            #     # reset lock and variable objects
            #     lock = False
            if id not in line and lock:
                # get the count, mean and fraction of significant scores 
                # sig_score = [x for x in iscore if float(x) >= 0.5]
                # nonsig_score = [x for x in iscore if float(x) < 0.5]
                res[id]['oyang_count'] = len(count)
                res[id]['oyang_max'] = max(iscore) if len(iscore) > 0 else 0
                res[id]['oyang_mostsig'] = max(count) if len(count) > 0 else 0
                res[id]['oyang_mean'] = mean(iscore) if len(iscore) > 0 else 0
                # reset lock and variable objects
                iscore, count = [], []
                lock = False
    
                
    result = pd.DataFrame(res)
    result = result.transpose()
    result.index.names = ['realNames']
    
    # convert ID to orgID
    # if result.index[0] is '1':
    if not TrueID:
        # print(result.index[0])
        # result.index = result.index.astype('int64')
        path = path.replace('\\', '/')
        pattern = '\/\w{2}\.'
        print(path)
        org = re.findall(pattern, path)
        org = org[0].strip('/').strip('.')
        map_url = 'nameMapping_additional/yinoyang/' + org + '.csv'
        
        mapper = pd.read_csv(map_url, index_col=0)
        result = pd.merge(mapper, result, left_index=True, right_index=True, how='inner')
        url = path.replace('.fa', '.csv')
        result.to_csv(url, index=False)
        print(result.head())
    else:
        url = path.replace('.fa', '.csv')
        result.to_csv(url)

path = 'outdata_latest/additional/yinyang/Sc.fa'
with open(path,'r') as inFile:    # line used for single file
    yinoyang(inFile, path, False)
    
    
# import os, glob
# path = 'outdata_latest/yinoyang'
# if os.path.exists(path):
#     files = glob.glob(os.path.join(path, '*.fa'))
#     for path in files:    
#         # readin input file
#         file = open(path, 'r')
#         inFile = file.readlines()
#         file.close()
#         yinoyang(inFile, path)
#         print(path + ' completed successfully!')