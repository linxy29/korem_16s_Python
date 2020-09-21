## Get the otus ratio table and a dictionary with co-occurrence rate(with similarity rate 90%).
import timeit
import numpy as np
import pandas as pd

start = timeit.default_timer()

similarity = float(input('The cutoff of similarity is: '))
data_used = input('Use log of count(log) or count(count): ')

while data_used not in ['log', 'count']:
    data_used = input('Please enter <log> or <count> ')

# 1. get copies numbers file
## counts of otus with minimal processing (basically just uniting the ends)
raw_otu_l10_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/raw_otu_l10.df')

## get log(add 1 to each value and take log)

logOtusRatio_df = raw_otu_l10_df.apply(lambda x : np.log(x + 1), axis = 0)
#print(logOtusRatio_df.head())

if data_used == 'log':
    data_df = logOtusRatio_df
else:
    data_df = raw_otu_l10_df


# 2. function of calculating copies ratios

# function to calculate ratio
def getRatio(list1, list2):
    res = []
    for i in range(len(list1)):
        if list1[i]==0 or list2[i]==0:
            res.append(0)
        else:
            res.append(list1[i]/list2[i])
    return res

# 3. function to select high coocurrence otus
def getCoocur(list1, list2):
    a = (list1 > 0)
    b =(list2 > 0)
    return (a&b).sum() / (a|b).sum()

# 4. get ratio table and coocur_rate
index_names = data_df.index.values
otusRatio_df = pd.DataFrame(index = index_names)

Coocur_rate = {}

## get pairs_list and add column
with open('D:/Codes/korem_16s/Data/usearch_allpairsGlobal/results' + str(int(similarity*100)) + '.useout') as read_file:
    #pairs_list = []
    for line in read_file:
        pair = line[:-1].split("\t")
        #print(line)
        #print(pair)
        if getCoocur(data_df[pair[0]], data_df[pair[1]]) >= 0.7:
            # print(getCoocur(data_df[item[0]], data_df[item[1]]))
            otusRatio_df[pair[0] + "-" + pair[1]] = getRatio(data_df[pair[0]], data_df[pair[1]])
            Coocur_rate[pair[0] + "-" + pair[1]] = getCoocur(data_df[pair[0]], data_df[pair[1]])

otusRatio_df.to_pickle('D:/Codes/korem_16s/Data/Real_Data/otusRatio' + str(int(similarity*100)) + '_' + data_used + '.df')

first2pairs = {k: Coocur_rate[k] for k in list(Coocur_rate)[:2]}
print(first2pairs)

# Save
np.save('D:/Codes/korem_16s/Data/Real_Data/Coocur_rate' + str(int(similarity*100)) + '_' + data_used + '.npy', Coocur_rate)

# Load
Coocur_rate_dic = np.load('D:/Codes/korem_16s/Data/Real_Data/Coocur_rate' + str(int(similarity*100)) + '_' + data_used + '.npy',allow_pickle='TRUE').item()
first2pairs = {k: Coocur_rate_dic[k] for k in list(Coocur_rate_dic)[:2]}
print(first2pairs)

stop = timeit.default_timer()
print('Time: ', stop - start)