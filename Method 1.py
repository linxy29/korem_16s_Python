import numpy as np
import pandas as pd

## Import Data
# data with map rate greater or equal to 97%
mapRate97_res = pd.read_csv("D:/Codes/korem_16s/Data/method1/mapRate97_res.csv")
# counts of otus with minimal processing (basically just uniting the ends)
raw_otu_l10_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/raw_otu_l10.df')
# growth estimates based on metagenomics for the same samples
P2T_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/P2T.df')
# sample information
genomes_metadata = pd.read_csv("D:/Codes/korem_16s/Data/Real_Data/genomes_metadata.tsv", delimiter='\t')

## 1. filter otus mapped to different species
onespe_query_list = []
query_list = mapRate97_res['query'].unique()
for item in query_list:
    spe_num = mapRate97_res.loc[mapRate97_res['query'] == item]['target_strain'].unique()
    if len(spe_num) == 1:
        onespe_query_list.append(item)

print(onespe_query_list[:10])
print('length of onespe query:' + str(len(onespe_query_list)))

onespe_mapRate97_df = mapRate97_res.loc[mapRate97_res['query'].isin(onespe_query_list)].sort_values('query')

## 2. combine otus(with otus in more than one records)
unique_recordID = onespe_mapRate97_df['target_recordId'].unique().tolist()
otusX_dic = {}
i = 1
for item in unique_recordID:
    sort_df = onespe_mapRate97_df.loc[onespe_mapRate97_df['target_recordId'] == item].sort_values('map_seqStart')
    sort_df['seq_diff'] = sort_df['map_seqStart'].diff()
    sort_df.fillna(0, inplace=True)
    # i += 1
    list = [item]
    for index, row in sort_df.iterrows():
        if row['seq_diff'] < 100:
            list.append(row['query'])
        else:
            if sorted(list) not in otusX_dic.values():
                otusX_dic['X' + str(i)] = sorted(list)
                i += 1
            list = [item]
            list.append(row['query'])

print("lenghth of outsX_dic:" + str(len(otusX_dic)))

# change to dataframe

otusX_df = pd.DataFrame(index = otusX_dic.keys(), columns = ['recordID', 'otus'])
for key, value in otusX_dic.items():
    otusX_df.loc[key, 'recordID'] = value[0]
    otusX_df.loc[key, 'otus'] = value[1:]

otusX_df.reset_index(inplace=True)
print(otusX_df.head())

## 4. get XsRatio table
## Xs numbers
for index, row in otusX_df.iterrows():
    raw_otu_l10_df[row['index']] = raw_otu_l10_df[row['otus']].sum(axis=1)

print(raw_otu_l10_df.head())

ID_list = otusX_df['recordID'].unique().tolist()
print(ID_list[:10])
print(len(ID_list))

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
XsRatio_df = pd.DataFrame(index = index_names)
print(XsRatio_df.head())

for item in ID_list:
    unique_query = otusX_df.loc[otusX_df['recordID'] == item]['index'].tolist()
    if len(unique_query)> 1:
        #print(len(unique_query))
        A_otus = otusX_df.loc[otusX_df['index'] == unique_query[0]]['otus'].tolist()
        B_otus = otusX_df.loc[otusX_df['index'] == unique_query[1]]['otus'].tolist()
        A_num = len(A_otus[0])
        B_num = len(B_otus[0])
        A_list = raw_otu_l10_df[unique_query[0]]/A_num
        B_list = raw_otu_l10_df[unique_query[1]]/B_num
        XsRatio_df[item] = getRatio(A_list, B_list)

print(XsRatio_df.head())

## 5. change sample name

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

XsRatio_df2 = XsRatio_df.rename(index=raw_otu_samplename)

## select sample with 'FD'
index = XsRatio_df2.index.tolist()
print("length of index: " + str(len(index)))
sub = 'FD'
with_FD = [s for s in index if sub in s]
print(with_FD[:10])
print("length of with_FD: " + str(len(with_FD)))
XsRatio_df3 = XsRatio_df2.loc[with_FD,:]
print("number of rows: " + str(XsRatio_df3.shape[0]))

# duplicate sample names
import re
XsRatio_index = XsRatio_df3.index.tolist()

otus_FD_name = []
for name in XsRatio_index:
    substr = re.findall(r"FD\d+",name)
    otus_FD_name.append(substr[0])

print(otus_FD_name[:20])
print(XsRatio_index[:10])

## check duplicates
duplicates = []
for item in otus_FD_name:
    if otus_FD_name.count(item) > 1:
        duplicates.append(item)
print(duplicates)

XsRatio_df3.index = otus_FD_name

P2T_FD_name = P2T_df.index.to_series().str.split("_").str[0].tolist()
P2T_df.index = P2T_FD_name

## change NaN into 0
P2T_df.fillna(0,inplace = True)

## select samples in both two dataframe
## use 'otusRatio_df3' and 'P2T_df'
ind = [item for item in otus_FD_name if item in P2T_FD_name]
print("length of ind:" + str(len(ind)))

## remove those are duplicate
ind2 = [item for item in ind if item not in duplicates]
print("length of ind2:" + str(len(ind2)))

com_otusRatio_df = XsRatio_df3.loc[ind2,:]
print("number of rows: " + str(com_otusRatio_df.shape[0]))
com_P2T_df = P2T_df.loc[ind2,:]
print("number of rows: " + str(com_P2T_df.shape[0]))

## 6. map two dataframe based on their species
otusRatio_IDvsStrain = {}
ID_list = com_otusRatio_df.columns.tolist()
print('the length of ID_list:' + str(len(ID_list)))
for item in ID_list:
    strain_list = onespe_mapRate97_df.loc[onespe_mapRate97_df['target_recordId'] == item]['target_strain'].tolist()
    if strain_list[0] in otusRatio_IDvsStrain:
        otusRatio_IDvsStrain[strain_list[0]].append(item)
    else:
        otusRatio_IDvsStrain[strain_list[0]] = []
        otusRatio_IDvsStrain[strain_list[0]].append(item)

print(otusRatio_IDvsStrain)

P2T_IDvsStrain = {}
P2TID_list = com_P2T_df.columns.tolist()
for item in P2TID_list:
    strain_list = genomes_metadata.loc[genomes_metadata['SegalLabID'] == int(item), 'Species'].tolist()
    if strain_list[0] in P2T_IDvsStrain:
        P2T_IDvsStrain[strain_list[0]].append(item)
    else:
        P2T_IDvsStrain[strain_list[0]] = []
        P2T_IDvsStrain[strain_list[0]].append(item)

print(P2T_IDvsStrain)

P2TSet = set(P2T_IDvsStrain)
RatioSet = set(otusRatio_IDvsStrain)

for name in RatioSet.intersection(P2TSet):
    print(name, P2T_IDvsStrain[name], otusRatio_IDvsStrain[name])

from scipy import stats
correlation_df = pd.DataFrame(index = ['3001', '1727'], columns = ['pearsonCor', 'pearsonPvalue', 'spearmanCor', 'spearmanPvalue'])

otusR = com_otusRatio_df['GCF_003312465.1']
for item in ['3001', '1727']:
    P2T = com_P2T_df[item]
    correlation_df.loc[item,'pearsonCor'] = stats.pearsonr(otusR,P2T)[0]
    correlation_df.loc[item,'pearsonPvalue'] = stats.pearsonr(otusR,P2T)[1]
    correlation_df.loc[item,'spearmanCor'] = stats.spearmanr(otusR,P2T)[0]
    correlation_df.loc[item,'spearmanPvalue'] = stats.spearmanr(otusR,P2T)[1]

print(correlation_df)