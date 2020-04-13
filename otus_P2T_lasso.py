## Get the dictionary include otus pairs selected by lasso

import numpy as np
import pandas as pd
import re
from sklearn.linear_model import LassoCV

# 1. import ratio table and coocurrence rate dictionary
otusRatio_df = pd.read_pickle("D:/Codes/korem_16s/Data/Real_Data/otusRatio2.df")
Coocur_rate = np.load('D:/Codes/korem_16s/Data/Real_Data/Coocur_rate.npy',allow_pickle='TRUE').item()

# 2. get PTR file
P2T_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/P2T.df')

# 3. change sample names
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

otusRatio_df2 = otusRatio_df.rename(index=raw_otu_samplename)
index = otusRatio_df2.index.tolist()

## select sample with 'FD'
sub = 'FD'
with_FD = [s for s in index if sub in s]
otusRatio_df3 = otusRatio_df2.loc[with_FD,:]

otusRatio_index = otusRatio_df3.index.tolist()
otus_FD_name = []
for name in otusRatio_index:
    substr = re.findall(r"FD\d+",name)
    otus_FD_name.append(substr[0])

## check duplicates
duplicates = []
for item in otus_FD_name:
    if otus_FD_name.count(item) > 1:
        duplicates.append(item)

otusRatio_df3.index = otus_FD_name

P2T_FD_name = list(P2T_df.index.to_series().str.split("_").str[0])
P2T_df.index = P2T_FD_name

## change NaN into 0
P2T_df.fillna(0,inplace = True)

## select samples in both two dataframe
## use 'otusRatio_df3' and 'P2T_df'
ind = [item for item in otus_FD_name if item in P2T_FD_name]

## remove those are duplicate
ind2 = [item for item in ind if item not in duplicates]

com_otusRatio_df = otusRatio_df3.loc[ind2,:]
com_P2T_df = P2T_df.loc[ind2,:]

# 4. lasso
## select otusRatio
cutoff = 0.90
selected_otus = []
for key, value in Coocur_rate.items():
 if value >= cutoff:
  selected_otus.append(key)

selected_otusRatio_df = com_otusRatio_df.loc[:, selected_otus]

# Save
selected_otusRatio_df.to_pickle('D:/Codes/korem_16s/Data/Real_Data/com_otusRatio90_df.npy')
com_P2T_df.to_pickle('D:/Codes/korem_16s/Data/Real_Data/com_P2T_df.npy')

selected_otusRatio_paris = selected_otusRatio_df.columns
lasso_otusPairs = {}
for column in com_P2T_df:
    model = LassoCV(cv=10,max_iter=1000).fit(selected_otusRatio_df, com_P2T_df[column])
    coef = model.coef_
    coef_notZero = coef!=0
    lasso_otusPairs[column] = selected_otusRatio_paris[coef_notZero].tolist()

# Save
np.save('D:/Codes/korem_16s/Data/Real_Data/lasso_otusRatio.npy', lasso_otusPairs)