# using taxonomic assignment and dot plots to validate method 2's results

import numpy as np
import pandas as pd
import random
import re
import itertools

# 1. taxonomic assignment
## rdp database
with open('D:/Codes/korem_16s/Data/lasso_validation/sintax_rdp_res.sintax') as read_file:
    sintax_rdp_dict = {}
    for line in read_file:
        sin_line = line[:-1].split("\t")
        #print(sin_line)
        taxa_split = sin_line[-1].split(",")
        #print(taxa_split)
        sintax_rdp_dict[sin_line[0]] = taxa_split

## ltp database
with open('D:/Codes/korem_16s/Data/lasso_validation/sintax_ltp_res.sintax') as read_file:
    sintax_ltp_dict = {}
    for line in read_file:
        sin_line = line[:-1].split("\t")
        #print(sin_line)
        taxa_split = sin_line[-1].split(",")
        #print(taxa_split)
        sintax_ltp_dict[sin_line[0]] = taxa_split

lasso_otusRatio_dic = np.load('D:/Codes/korem_16s/Data/Lasso/lasso_otusRatio.npy',allow_pickle='TRUE').item()

## rdp
taxa_per_dict = {}
for key, value in lasso_otusRatio_dic.items():
    # get all otus
    if len(value)!= 0:
        otus_list = []
        for item in value:
            each_split = item.split("-")
            otus_list.append(each_split[0])
            otus_list.append(each_split[1])
        unique_otus = list(set(otus_list))
        # calculate numbers of otus belongs to same strain
        same_num = 0
        diff_num = 0
        all_num = 0
        for pair in itertools.combinations(unique_otus, 2):
            all_num += 1
            lenA = len(sintax_rdp_dict[pair[0]])
            lenB = len(sintax_rdp_dict[pair[1]])
            if lenA == lenB:
                if sintax_rdp_dict[pair[0]][lenA-1] == sintax_rdp_dict[pair[1]][lenA-1]:
                    #taxa = sintax_rdp_dict[pair[0]][lenA-1]
                    same_num += 1
                else:
                    diff_num += 1
            else:
                diff_num += 1
        #taxa_per_dict[key] = [same_num, diff_num, all_num, taxa]
        taxa_per_dict[key] = [same_num, diff_num, all_num, diff_num/all_num]
print(taxa_per_dict)

##ltp
taxa_per_dict = {}
for key, value in lasso_otusRatio_dic.items():
    # get all otus
    if len(value)!= 0:
        otus_list = []
        for item in value:
            each_split = item.split("-")
            otus_list.append(each_split[0])
            otus_list.append(each_split[1])
        unique_otus = list(set(otus_list))
        # calculate numbers of otus belongs to same strain
        same_num = 0
        diff_num = 0
        all_num = 0
        for pair in itertools.combinations(unique_otus, 2):
            all_num += 1
            lenA = len(sintax_ltp_dict[pair[0]])
            lenB = len(sintax_ltp_dict[pair[1]])
            if lenA == lenB:
                if sintax_ltp_dict[pair[0]][lenA-1] == sintax_ltp_dict[pair[1]][lenA-1]:
                    #taxa = sintax_ltp_dict[pair[0]][lenA-1]
                    same_num += 1
                else:
                    diff_num += 1
            else:
                diff_num += 1
        #taxa_per_dict[key] = [same_num, diff_num, all_num, taxa]
        taxa_per_dict[key] = [same_num, diff_num, all_num, diff_num/all_num]
print(taxa_per_dict)

# 2. draw dot plots
P2T_pred_df = pd.read_csv("D:/Codes/korem_16s/Data/Lasso/P2T_pred_df.csv")
P2T_true_df = pd.read_csv("D:/Codes/korem_16s/Data/Lasso/P2T_true_df.csv")

import matplotlib.pyplot as plt

#for column in P2T_pred_df:
for column in P2T_pred_df.iloc[:,1:49]:
    pred_rank = P2T_pred_df[column][P2T_pred_df[column]!=0].rank(ascending = 1)
    true_rank = P2T_true_df[column][P2T_true_df[column]!=0].rank(ascending = 1)
    plt.scatter(pred_rank, true_rank)
    plt.title(column)
    plt.xlabel('Predicted P2T')
    plt.ylabel('True P2T')
    plt.show()