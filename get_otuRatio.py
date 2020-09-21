## Get the otus ratio table and a dictionary with co-occurrence rate(with similarity rate 90%).

import pandas as pd
import numpy as np

# 1. get pairs_list
with open('D:/Codes/korem_16s/Data/usearch_allpairsGlobal/results90.useout') as read_file:
    pairs_list = []
    for line in read_file:
        pair = line[:-1].split("\t")
        #print(line)
        #print(pair)
        pairs_list.append(pair)

# 2. get copies numbers file
raw_otu_l10_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/raw_otu_l10.df')
raw_otu_l10_df # counts of otus with minimal processing (basically just uniting the ends)

# 4. function of calculating copies ratios
#print(raw_otu_l10_df['RO1'])

# function to calculate ratio
def getRatio(list1, list2):
    res = []
    for i in range(len(list1)):
        if list1[i]==0 or list2[i]==0:
            res.append(0)
        else:
            res.append(list1[i]/list2[i])
    return res

# 5. function to select high coocurrence otus
def getCoocur(list1, list2):
    a = (list1 > 0)
    b =(list2 > 0)
    return (a&b).sum() / (a|b).sum()

# 6. get ratio table
index_names = raw_otu_l10_df.index.values
otusRatio_df = pd.DataFrame(index = index_names)

## add column
for item in pairs_list:
    if getCoocur(raw_otu_l10_df[item[0]], raw_otu_l10_df[item[1]]) >= 0.7:
        #print(getCoocur(raw_otu_l10_df[item[0]], raw_otu_l10_df[item[1]]))
        otusRatio_df[item[0] + "-" + item[1]] = getRatio(raw_otu_l10_df[item[0]], raw_otu_l10_df[item[1]])

#print(otusRatio_df.head())
otusRatio_df.to_pickle('D:/Codes/korem_16s/Data/Real_Data/otusRatio90.df')

#test = pd.read_pickle("D:/Codes/korem_16s/Data/Real_Data/otusRatio.df")
#print(test.head())


# otusRatio_df = pd.read_pickle("D:/Codes/korem_16s/Data/Real_Data/otusRatio.df")
#
# for item in pairs_list[20000001:]:
#     if getCoocur(raw_otu_l10_df[item[0]], raw_otu_l10_df[item[1]]) >= 0.7:
#         #print(getCoocur(raw_otu_l10_df[item[0]], raw_otu_l10_df[item[1]]))
#         otusRatio_df[item[0] + "-" + item[1]] = getRatio(raw_otu_l10_df[item[0]], raw_otu_l10_df[item[1]])
#
# otusRatio_df.to_pickle('D:/Codes/korem_16s/Data/Real_Data/otusRatio2.df')

# 7. get co-occurrence dictionary
Coocur_rate = {}
for item in pairs_70Coocur_list:
    Coocur_rate[item[0] + "-" + item[1]] = getCoocur(raw_otu_l10_df[item[0]], raw_otu_l10_df[item[1]])

first2pairs = {k: Coocur_rate[k] for k in list(Coocur_rate)[:2]}
print(first2pairs)

# Save
np.save('D:/Codes/korem_16s/Data/Real_Data/Coocur_rate90.npy', Coocur_rate)

# Load
Coocur_rate_dic = np.load('D:/Codes/korem_16s/Data/Real_Data/Coocur_rate90.npy',allow_pickle='TRUE').item()
first2pairs = {k: Coocur_rate_dic[k] for k in list(Coocur_rate_dic)[:2]}
print(first2pairs)


