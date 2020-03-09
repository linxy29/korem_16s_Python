# get all pairs of queries with 95% similarity
import re

with open('D:/Codes/korem_16s/Data/usearch_allpairsGlobal/raw_otus_sub1.fa') as fas_file:
    sequences = fas_file.read()
    sequences = re.split("^>", sequences, flags=re.MULTILINE)  # Only splits string at the start of a line.
    del sequences[0]
rawOtus_sub1_dic = {}
for fasta in sequences:
    header, sequence = fasta.split("\n", 1)  # Split each fasta into header and sequence.
    sequence = sequence.replace("\n", "")
    rawOtus_sub1_dic[header] = sequence

# first2pairs = {k: rawOtus_sub1_dic[k] for k in list(rawOtus_sub1_dic)[:2]}
# print(first2pairs)

with open('D:/Codes/korem_16s/Data/usearch_allpairsGlobal/raw_otus_sub2.fa') as fas_file:
    sequences = fas_file.read()
    sequences = re.split("^>", sequences, flags=re.MULTILINE)  # Only splits string at the start of a line.
    del sequences[0]
rawOtus_sub2_dic = {}
for fasta in sequences:
    header, sequence = fasta.split("\n", 1)  # Split each fasta into header and sequence.
    sequence = sequence.replace("\n", "")
    rawOtus_sub2_dic[header] = sequence

# first2pairs = {k: rawOtus_sub2_dic[k] for k in list(rawOtus_sub2_dic)[:2]}
# print(first2pairs)
#print(len(rawOtus_sub2_dic))

from Bio import pairwise2
selected_pairs = []
for key1, value1 in rawOtus_sub1_dic.items():
    #print(key1)
    #print(value1)
    for key2, value2 in rawOtus_sub2_dic.items():
        #length = len(value1)  # lenghth of all queries are all 214
        #print(length)
        score = pairwise2.align.globalxx(value1,value2, score_only=True)
        if score/214 >= 0.95:
            selected_pairs.append([key1, key2])

#print(selected_pairs)

with open('D:/Codes/korem_16s/Data/usearch_allpairsGlobal/results_sub3.useout',"w") as writefile:
    for pairs in selected_pairs:
        line = pairs[0] + "\t" + pairs[1] + "\n"
        writefile.write(line)