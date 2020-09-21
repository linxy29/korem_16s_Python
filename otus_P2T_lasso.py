## Get the dictionary include otus pairs selected by lasso

import numpy as np
import pandas as pd
import re
from sklearn.linear_model import LassoCV
from korem_16s_Python import Functions

similarity = float(input('The cutoff of similarity is: '))
data_used = input('Use log of count(log) or count(count): ')
concurr = float(input('The cutoff of concurrence rate is: '))

while data_used not in ['log', 'count']:
    data_used = input('Please enter <log> or <count> ')

# 1. import ratio table and coocurrence rate dictionary
otusRatio_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/otusRatio' + str(int(similarity*100)) + '_' + data_used + '.df')  # the cutoff of similarity is 95
Coocur_rate = np.load('D:/Codes/korem_16s/Data/Real_Data/Coocur_rate' + str(int(similarity*100)) + '_' + data_used + '.npy',allow_pickle='TRUE').item()

# 2. get PTR file
P2T_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/P2T.df')

# 3. change sample names
## sample name modified dictionary
com_otusRatio_df, com_P2T_df = Functions.ChangeSelect_FD(otusRatio_df, P2T_df)

# 4. lasso
## select otusRatio
selected_otus = []
for key, value in Coocur_rate.items():
 if value >= concurr:
  selected_otus.append(key)

selected_otusRatio_df = com_otusRatio_df.loc[:, selected_otus]

# Save
selected_otusRatio_df.to_pickle('D:/Codes/korem_16s/Data/Lasso/S' + str(int(similarity*100)) + '_' + data_used + '_com_otusRatio_df.npy')
com_P2T_df.to_pickle('D:/Codes/korem_16s/Data/Lasso/S' + str(int(similarity*100)) + '_com_P2T_df.npy')

selected_otusRatio_paris = selected_otusRatio_df.columns
lasso_otusPairs = {}
for column in com_P2T_df:
    filt_P2T = com_P2T_df.loc[com_P2T_df[column] > 0][column]
    filt_otusRatio = selected_otusRatio_df[com_P2T_df[column] > 0]
    #print(len(filt_P2T))
    if len(filt_P2T) >= 10:
     model = LassoCV(cv=10,max_iter=1000).fit(filt_otusRatio, filt_P2T)
     coef = model.coef_
     coef_notZero = coef!=0
     lasso_otusPairs[column] = selected_otusRatio_paris[coef_notZero].tolist()

print(lasso_otusPairs)
# Save
np.save('D:/Codes/korem_16s/Data/Lasso/S' + str(int(similarity*100)) + '_' + data_used + '_lasso_otusRatio' + str(int(concurr*100)) + '.npy', lasso_otusPairs)