# split one fasta file into several files
import re
from Bio.Seq import Seq
with open('D:/Codes/korem_16s/Data/rrnDB/rrnDB-5.6_16S_rRNA.fasta') as fas_file:
    sequences = fas_file.read()
    sequences = re.split("^>", sequences, flags=re.MULTILINE) # Only splits string at the start of a line.
    del sequences[0]
    #print(sequences[0:5])

fasta_dic1 = {}
for fasta in sequences:
    header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
    sequence = sequence.replace("\n","")
    header_split = header.split("|")
    NCBI_access = header_split[1]
    if NCBI_access not in fasta_dic1:
        fasta_dic1[NCBI_access] = []
    if header[-1].endswith('-'):
        seq = Seq(sequence)
        sequence = str(seq.reverse_complement())
    fasta_dic1[NCBI_access].append(sequence)

#print(fasta_dic1["GCF_000762265.1"][0])
#first2pairs = {k: fasta_dic1[k] for k in list(fasta_dic1)[:2]}
#print(first2pairs)

for key, value in fasta_dic1.items():
    f = open('D:/Codes/korem_16s/Data/rrnDB/fasta_sub_ID/%s.fasta' % key,'w')
    i = 0
    for item in value:
        i = i + 1
        header = ">" + key + "_" + str(i) + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
        sequence = item + "\n" # Replace newlines in sequence, remember to add one to the end.
        f.write(header + sequence)
    f.close()

# not finished yet
# split based on strain
fasta_dic2 = {}
for fasta in sequences:
    header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
    sequence = sequence.replace("\n","")
    header_split = header.split("|")
    strain = header_split[0]
    if strain not in fasta_dic2:
        fasta_dic2[strain] = []
    if header[-1].endswith('-'):
        seq = Seq(sequence)
        sequence = str(seq.reverse_complement())
    info = [header_split[1], sequence]
    fasta_dic2[strain].append(info)

#print(fasta_dic2["GCF_000762265.1"][0])
#first2pairs = {k: fasta_dic2[k] for k in list(fasta_dic1)[:2]}
#print(first2pairs)

for key, value in fasta_dic2.items():
    f = open('D:/Codes/korem_16s/Data/rrnDB/fasta_sub_strain/%s.fasta' % key,'w')
    i = 0
    for item in value:
        i = i + 1
        header = ">" + item[0] + "_" + str(i) + "|" + key + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
        sequence = item[1] + "\n" # Replace newlines in sequence, remember to add one to the end.
        f.write(header + sequence)
    f.close()