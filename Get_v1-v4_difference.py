# pairwise alignment for all combination of 16s copies
import re
from Bio.Seq import Seq
with open('D:/Codes/korem_16s/Data/rrnDB/rrnDB-5.6_16S_rRNA.fasta') as fas_file:
    sequences = fas_file.read()
    sequences = re.split("^>", sequences, flags=re.MULTILINE) # Only splits string at the start of a line.
    del sequences[0]
    #print(sequences[0:5])

fasta_dic = {}
for fasta in sequences:
    header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
    sequence = sequence.replace("\n","")
    header_split = header.split("|")
    NCBI_access = header_split[1]
    pos_line = re.split('\.\.|\s', header_split[-1])
    #print(pos_line)
    if NCBI_access not in fasta_dic:
        fasta_dic[NCBI_access] = []
    if header[-1].endswith('-'):
        seq = Seq(sequence)
        sequence = str(seq.reverse_complement())
    copyID = "c" + str(len(fasta_dic[NCBI_access])+1)
    sequence_info = pos_line[:2]
    sequence_info.append(copyID)
    sequence_info.append(sequence)
    #print(sequence_info)
    fasta_dic[NCBI_access].append(sequence_info)

#print(fasta_dic["GCF_000762265.1"][0])
#first2pairs = {k: fasta_dic[k] for k in list(fasta_dic)[:2]}
#print(first2pairs)

# filter records with more than one 16s copy and add copy ID
import itertools
from Bio import pairwise2

#first10000pairs = {k: fasta_dic[k] for k in list(fasta_dic)[:10000]}

res_list = []
for key, value in fasta_dic.items():
#for key, value in first10000pairs.items():
    if len(value) > 1:
        comb_seq = itertools.combinations(value, 2)
        for comb in comb_seq:
            #distance = abs(int(comb[1][0]) - int(comb[0][0]))
            distance = int(comb[1][0]) - int(comb[0][0])
            if distance > 150000:
                copies = comb[0][2] + "-" + comb[1][2]
                v1v2 = pairwise2.align.globalxx(comb[0][3][49:249],comb[1][3][49:249], score_only=True)
                v3v4 = pairwise2.align.globalxx(comb[0][3][449:799],comb[1][3][449:799], score_only=True)
                res = [key, copies, distance, 200-v1v2, 350-v3v4]
                res_list.append(res)
#print(res_list)

# write results
import csv
with open('D:/Codes/korem_16s/Data/pairwise_alig_res.csv',"w") as csvfile:
    writer = csv.writer(csvfile)
    # columns name
    writer.writerow(["record_id","copies","distance","v1v2_diff","v3v4_diff"])
    # positions and other information
    writer.writerows(res_list)