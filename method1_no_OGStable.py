import numpy as np
import pandas as pd

## 100% matched otus
# 1. get queries group
## 100% matched dataframe(not split)
map100_res_df = pd.read_csv("D:/Codes/korem_16s/Data/method1/map100_res.csv")
dup_target_df = pd.read_csv("D:/Codes/korem_16s/Data/method1/dup_target.csv")
# counts of otus with minimal processing (basically just uniting the ends)
raw_otu_l10_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/raw_otu_l10.df')
# growth estimates based on metagenomics for the same samples
P2T_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/P2T.df')
# sample information
genomes_metadata = pd.read_csv("D:/Codes/korem_16s/Data/Real_Data/genomes_metadata.tsv", delimiter='\t')

# 2. group otus and get numbers
## should introduce flexible start point: CB1, CB2
otus_com_dict = {}
i = 0
for index, row in dup_target_df.iterrows():
    query_group = map100_res_df.loc[(map100_res_df['target'] == dup_target_df.loc[index,'target']) & (map100_res_df['target_start'] < dup_target_df.loc[index,'target_start'] + 3) & (map100_res_df['target_start'] > dup_target_df.loc[index,'target_start'] - 3),'query'].tolist()
    if query_group not in otus_com_dict.values():
        i += 1
        otus_com_dict['CB'+str(i)] = query_group
print(otus_com_dict)

for key, value in otus_com_dict.items():
    raw_otu_l10_df[key] = raw_otu_l10_df[value].sum(axis=1)
print(raw_otu_l10_df.head())

# 3. replace otus with 'CB'
rev_otusCom_dict = {k: oldk for oldk, oldv in otus_com_dict.items() for k in oldv}
map100_res_df["query"].replace(rev_otusCom_dict, inplace=True)
print(map100_res_df.head())

map100_res_df.to_csv (r'D:/Codes/korem_16s/Data/method1/map100_replace_res.csv', index = False, header=True)

# 4. split target
map100_res_df[["target_strain", "target_recordId", "target_refId", "target_Chromosome", "target_other"]] = map100_res_df.target.str.split("|",expand=True)
map100_res_df[["target_Wholstart", 'rm1', "target_Wholend"]] = map100_res_df.target_other.str.split('.',expand=True)
map100_res_df['target_Wholend'] = map100_res_df['target_Wholend'].str.replace(r'\D', '')
map100_res_df.drop(['target_refId', 'rm1', 'target_other'], axis=1, inplace=True)

# 5. get ratio
def getRatio(list1, list2):
    res = []
    for i in range(len(list1)):
        if list1[i]==0 or list2[i]==0:
            res.append(0)
        else:
            if list1[i] > list2[i]:
                res.append(list1[i]/list2[i])
            else:
                res.append(list2[i]/list1[i])
    return res

index_names = raw_otu_l10_df.index.values
queriesRatio_df = pd.DataFrame(index = index_names)
print(queriesRatio_df.head())

map100_ID = map100_res_df['target_recordId'].unique().tolist()
for item in map100_ID:
    query_list = map100_res_df.loc[map100_res_df['target_recordId'] == item]['query'].tolist()
    unique_query = map100_res_df.loc[map100_res_df['target_recordId'] == item]['query'].unique().tolist()
    if len(unique_query)> 1:
        #print(len(unique_query))
        A_num = query_list.count(unique_query[0])
        B_num = query_list.count(unique_query[1])
        A_list = raw_otu_l10_df[unique_query[0]]/A_num
        B_list = raw_otu_l10_df[unique_query[1]]/B_num
        queriesRatio_df[item] = getRatio(A_list, B_list)

print(queriesRatio_df.head())

# 5. map otus and P2T with species
otusRatio_IDvsStrain = {}
ID_list = queriesRatio_df.columns.tolist()
for item in ID_list:
    strain_list = map100_res_df.loc[map100_res_df['target_recordId'] == item]['target_strain'].tolist()
    if strain_list[0] in otusRatio_IDvsStrain:
        otusRatio_IDvsStrain[strain_list[0]].append(item)
    else:
        otusRatio_IDvsStrain[strain_list[0]] = []
        otusRatio_IDvsStrain[strain_list[0]].append(item)

P2T_IDvsStrain = {}
P2TID_list = P2T_df.columns.tolist()
for item in P2TID_list:
    strain_list = genomes_metadata.loc[genomes_metadata['SegalLabID'] == int(item), 'Species'].tolist()
    if strain_list[0] in P2T_IDvsStrain:
        P2T_IDvsStrain[strain_list[0]].append(item)
    else:
        P2T_IDvsStrain[strain_list[0]] = []
        P2T_IDvsStrain[strain_list[0]].append(item)

P2TSet = set(P2T_IDvsStrain)
RatioSet = set(otusRatio_IDvsStrain)

for name in RatioSet.intersection(P2TSet):
    print(name, P2T_IDvsStrain[name], otusRatio_IDvsStrain[name])

# 6. change sample names
## raw_otu_l10_df

## sample name modified dictionary
raw_otu_samplename = {'E15_1646_972761': 'E15_972761_FD443',
 'E15_2018_461108': 'E15_461108_FD451',
 'E15_2032_969484': 'E15_969484_FD445',
 'E15_2060_90317': 'E15_90317_FD444',
 'E15_2066_898593': 'E15_898593_FD447',
 'E15_2130_574949': 'E15_574949_FD448',
 'E15_2236_505209': 'E15_505209_FD452',
 'E15_2502_521450': 'E15_521450_FD465',
 'E15_2570_924164': 'E15_924164_FD471',
 'E15_2596_704040': 'E15_704040_FD500',
 'E15_2864_70890': 'E15_70890_FD493',
 'E15_2980_864274': 'E15_864274_FD560',
 'E9_1366_858809': 'E9_858809_FD446',
 'E9_1712_206581': 'E9_206581_FD441',
 'E9_1714_547936': 'E9_547936_FD442',
 'E9_2060_90317': 'E9_90317_FD444',
 'E9_2132_216164': 'E9_216164_FD449',
 'E9_2134_715320': 'E9_715320_FD450',
 'E9_2156_924456': 'E9_924456_FD455',
 'E9_2180_812549': 'E9_812549_FD454',
 'E9_2184_944720': 'E9_944720_FD453',
 'E9_2314_979887': 'E9_979887_FD462',
 'E9_2316_465475': 'E9_465475_FD457',
 'E9_2340_433060': 'E9_433060_FD464',
 'E9_2342_951397': 'E9_951397_FD460',
 'E9_2344_653308': 'E9_653308_FD456',
 'E9_2372_600497': 'E9_600497_FD458',
 'E9_2376_21331': 'E9_21331_FD461',
 'E9_2400_523506': 'E9_523506_FD478',
 'E9_2454_105043': 'E9_105043_FD463',
 'E9_2490_234223': 'E9_234223_FD472',
 'E9_2492_935630': 'E9_935630_FD473',
 'E9_2498_656056': 'E9_656056_FD474',
 'E9_2506_304669': 'E9_304669_FD469',
 'E9_2508_43542': 'E9_43542_FD468',
 'E9_2554_709998': 'E9_709998_FD459',
 'E9_2556_189010': 'E9_189010_FD466',
 'E9_2562_431980': 'E9_431980_FD479',
 'E9_2566_211610': 'E9_211610_FD470',
 'E9_2574_174522': 'E9_174522_FD467',
 'E9_2576_711574': 'E9_711574_FD477',
 'E9_2580_975914': 'E9_975914_FD476',
 'E9_2582_795066': 'E9_795066_FD480',
 'E9_2594_381126': 'E9_381126_FD499',
 'E9_2614_143104': 'E9_143104_FD482',
 'E9_2632_310843': 'E9_310843_FD481',
 'E9_2648_748945': 'E9_748945_FD475',
 'E9_2660_212550': 'E9_212550_FD484',
 'E9_2674_573873': 'E9_573873_FD485',
 'E9_2678_519352': 'E9_519352_FD492',
 'E9_2724_88269': 'E9_88269_FD491',
 'E9_2726_671730': 'E9_671730_FD486',
 'E9_2728_529861': 'E9_529861_FD488',
 'E9_2734_93410': 'E9_93410_FD490',
 'E9_2736_663673': 'E9_663673_FD494',
 'E9_2744_377008': 'E9_377008_FD483',
 'E9_2806_203727': 'E9_203727_FD487',
 'E9_2822_536212': 'E9_536212_FD497',
 'E9_2824_917390': 'E9_917390_FD496',
 'E9_2828_355044': 'E9_355044_FD489',
 'E9_2866_428224': 'E9_428224_FD502',
 'E9_2884_132895': 'E9_132895_FD495',
 'E9_2894_939235': 'E9_939235_FD498',
 'E9_2896_672284': 'E9_672284_FD501',
 'E9_3122_66206': 'E9_66206_FD504'}

queriesRatio_df2 = queriesRatio_df.rename(index=raw_otu_samplename)

## select sample with 'FD'
index = queriesRatio_df2.index.tolist()
print("length of index: " + str(len(index)))
sub = 'FD'
with_FD = [s for s in index if sub in s]
print(with_FD[:10])
print("length of with_FD: " + str(len(with_FD)))
queriesRatio_df3 = queriesRatio_df2.loc[with_FD,:]
print(queriesRatio_df3.head())
print("number of rows: " + str(queriesRatio_df3.shape[0]))

import re
queriesRatio_index = queriesRatio_df3.index.tolist()

otus_FD_name = []
for name in queriesRatio_index:
    substr = re.findall(r"FD\d+",name)
    otus_FD_name.append(substr[0])

print(otus_FD_name[:20])
print(queriesRatio_index[:10])

## check duplicates
duplicates = []
for item in otus_FD_name:
    if otus_FD_name.count(item) > 1:
        duplicates.append(item)
print(duplicates)
## there are duplicate sample names

queriesRatio_df3.index = otus_FD_name

P2T_FD_name = list(P2T_df.index.to_series().str.split("_").str[0])
P2T_df.index = P2T_FD_name

## change NaN into 0
P2T_df.fillna(0,inplace = True)

## select samples in both two dataframe
## use 'otusRatio_df3' and 'P2T_df'
ind = [item for item in otus_FD_name if item in P2T_FD_name]
print("length of ind:" + str(len(ind)))

## remove those are duplicate
ind2 = [item for item in ind if item not in duplicates]
print(ind2[:20])
print("length of ind2:" + str(len(ind2)))

com_otusRatio_df = queriesRatio_df3.loc[ind2,:]
print("number of rows: " + str(com_otusRatio_df.shape[0]))
com_P2T_df = P2T_df.loc[ind2,:]
print("number of rows: " + str(com_P2T_df.shape[0]))

# 6. calculate correlation(need to modified)
from scipy import stats
index = ['2121-GCF_002005185.1', '3001-GCF_003312465.1', '1727-GCF_003312465.1', '147-GCF_003443395.1', '2661-GCF_003443395.1', '1300-GCF_003443395.1', '906-GCF_003443395.1', '2252-GCF_003443395.1', '2875-GCF_003443395.1']
correlation_df = pd.DataFrame(index = index, columns = ['pearsonCor', 'pearsonPvalue', 'spearmanCor', 'spearmanPvalue'])

#for item in ['3001', '1727']:
for index, row in correlation_df.iterrows():
    index_list = index.split('-')
    P2T = com_P2T_df[index_list[0]]
    otusR = com_otusRatio_df[index_list[1]]
    correlation_df.loc[index,'pearsonCor'] = stats.pearsonr(otusR,P2T)[0]
    correlation_df.loc[index,'pearsonPvalue'] = stats.pearsonr(otusR,P2T)[1]
    correlation_df.loc[index,'spearmanCor'] = stats.spearmanr(otusR,P2T)[0]
    correlation_df.loc[index,'spearmanPvalue'] = stats.spearmanr(otusR,P2T)[1]

print(correlation_df)