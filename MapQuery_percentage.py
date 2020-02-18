# calculate the percentage of 100% mapped query
import pandas as pd
raw_otu_l10_df = pd.read_pickle('D:/Codes/korem_16s/Data/Real_Data/raw_otu_l10.df')
raw_otu_l10_df # counts of otus with minimal processing (basically just uniting the ends)

sum_col = raw_otu_l10_df.sum(axis = 0, skipna = True)
#sum_col
#sum_col['RO1']
#sum_col.sum()

with open('D:/Codes/korem_16s/Data/100MapQuery.txt') as input_file:
    map_query = []
    for line in input_file:
        edline = line.rstrip('\n')
        map_query.append(edline)

#map_query

mapped_reads = 0
for item in map_query:
    mapped_reads = sum_col[item]

total_reads = sum_col.sum()
print(mapped_reads/total_reads)