# calculate ratio

##！！test version(test_results.useout)

# 1. get pairs_list
with open('D:/Codes/korem_16s/Data/usearch_allpairsGlobal/test_results.useout') as read_file:
    pairs_list = []
    for line in read_file:
        pair = line[:-1].split("\t")
        pairs_list.append(pair)

#print(pairs_list[:10])
#print(len(pairs_list))

# 2. get copies numbers file
import pandas as pd
raw_otu_l10_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/raw_otu_l10.df')
raw_otu_l10_df # counts of otus with minimal processing (basically just uniting the ends)

# distribution of numbers of 0 in the raw_otu table
print(raw_otu_l10_df.isin([0]).sum().describe())

# 3. get PTR file
P2T_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/P2T.df')
P2T_df # growth estimates based on metagenomics for the same samples

# distribution of numbers of NaN in the P2T table
print(P2T_df.isna().sum().describe())

# 4. map otus_ratio with P2T through sample name
otu_name = raw_otu_l10_df.index.values
P2T_name = P2T_df.index.values
otu_name.sort()
P2T_name.sort()
#print(otu_name[:100])
#print(P2T_name[:100])

# 5. function of calculating copies ratios
list1 = raw_otu_l10_df['RO1']
list2 = raw_otu_l10_df['RO2']
res = []
for i in range(len(list1)):
    #print(list2[i]==0 and list2[i]==0)
    if list1[i]==0 or list2[i]==0:
        res.append(0)
    else:
        res.append(list1[i]/list2[i])

#print(list1[:5])
#print(list2[:5])
#print(res[:5])

def getRatio(list1, list2):
    res = []
    for i in range(len(list1)):
        if list1[i]==0 or list2[i]==0:
            res.append(0)
        else:
            res.append(list1[i]/list2[i])
    return res

test = getRatio(raw_otu_l10_df['RO1'], raw_otu_l10_df['RO2'])
#print(test[:5])
#print(raw_otu_l10_df['RO1'][:5])
#print(raw_otu_l10_df['RO2'][:5])

# 6. get ratio table
#print(len(pairs_list))
index_names = raw_otu_l10_df.index.values
otusRatio_df = pd.DataFrame(index = index_names)
## add column
for item in pairs_list:
    otusRatio_df[item[0] + "-" + item[1]] = getRatio(raw_otu_l10_df[item[0]], raw_otu_l10_df[item[1]])

print(otusRatio_df.head())