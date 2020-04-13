## using 10-fold cross validation to get predicted P2T and calculate pearson correlation, spearman correlation

import numpy as np
import pandas as pd
import random
from sklearn.linear_model import LinearRegression

# Load data
com_otusRatio90_df = pd.read_pickle('D:/Codes/korem_16s/Data/Lasso/com_otusRatio90_df.npy')
lasso_otusRatio_dic = np.load('D:/Codes/korem_16s/Data/Lasso/lasso_otusRatio.npy',allow_pickle='TRUE').item()
com_P2T_df = pd.read_pickle('D:/Codes/korem_16s/Data/Lasso/com_P2T_df.npy')

## make sure have at least 48 sample
zero_num = (com_P2T_df>0).sum()
selected_strain1 = zero_num[zero_num>957*0.05].index.tolist()
print('length of selected_strain:' + str(len(selected_strain1)))

## select strain that have corresponding otus pairs
selected_strain2 = [s for s in selected_strain1 if len(lasso_otusRatio_dic[s])>0]
print('length of selected_strain:' + str(len(selected_strain2)))

# P2T_pred table
index_names = com_P2T_df.index.values
P2T_pred_df = pd.DataFrame(index = index_names, columns = selected_strain2)

## add column
for column in com_P2T_df[selected_strain2]:
#for column in com_P2T_df[['882', '2861', '310', '3001', '1256']]:
    ## 10-fold-CV sampling
    sample = com_P2T_df[com_P2T_df[column]>0].index.tolist()
    fold_size = int(len(sample)//10)
    random.shuffle(sample)
    tenFold_sample = [sample[:fold_size],sample[fold_size:fold_size*2],sample[fold_size*2:fold_size*3],sample[fold_size*3:fold_size*4],sample[fold_size*4:fold_size*5],sample[fold_size*5:fold_size*6],sample[fold_size*6:fold_size*7],sample[fold_size*7:fold_size*8],sample[fold_size*8:fold_size*9],sample[fold_size*9:]]
    for fold in tenFold_sample:
        test = fold
        train = [s for s in sample if s not in test]
        y_train = com_P2T_df.loc[train,:][column]
        x_train = com_otusRatio90_df.loc[train,:][lasso_otusRatio_dic[column]]
        model = LinearRegression().fit(x_train, y_train)
        x_test = com_otusRatio90_df.loc[test,:][lasso_otusRatio_dic[column]]
        P2T_pred_df.loc[test, column] = model.predict(x_test)

P2T_pred_df.fillna(0,inplace = True)

# pearson and spearman correlation
from scipy import stats

index_names = P2T_pred_df.columns.values
correlation_df = pd.DataFrame(index = index_names, columns = ['pearsonCor', 'pearsonPvalue', 'spearmanCor', 'spearmanPvalue', 'coef_num', 'zeroCoef_num'])
#print(correlation_df.head())

for column in P2T_pred_df:
    pred = P2T_pred_df[column]
    true = com_P2T_df[column]
    correlation_df.loc[column,'pearsonCor'] = stats.pearsonr(pred,true)[0]
    correlation_df.loc[column,'pearsonPvalue'] = stats.pearsonr(pred,true)[1]
    correlation_df.loc[column,'spearmanCor'] = stats.spearmanr(pred,true)[0]
    correlation_df.loc[column,'spearmanPvalue'] = stats.spearmanr(pred,true)[1]
    correlation_df.loc[column,'coef_num'] = len(lasso_otusRatio_dic[column])
    correlation_df.loc[column,'zeroCoef_num'] = 7455 - len(lasso_otusRatio_dic[column])

print(correlation_df.sort_values('spearmanCor'))

P2T_pred_df.to_csv('D:/Codes/korem_16s/Data/Lasso/P2T_pred_df.csv')
com_P2T_df[selected_strain2].to_csv('D:/Codes/korem_16s/Data/Lasso/P2T_true_df.csv')

