# make rrnDB_fasta file suitable for blastdb
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
    if NCBI_access not in fasta_dic:
        fasta_dic[NCBI_access] = []
        fasta_dic[NCBI_access].append(header_split[0])
    if header[-1].endswith('-'):
        seq = Seq(sequence)
        sequence = str(seq.reverse_complement())
    fasta_dic[NCBI_access].append(sequence)

#print(fasta_dic["GCF_000762265.1"][0])
#first2pairs = {k: fasta_dic[k] for k in list(fasta_dic)[:2]}
#print(first2pairs)

with open('D:/Codes/korem_16s/Data/rrnDB_blastdb.fasta',"w") as writefile:
    for key, value in fasta_dic.items():
        i = 0
        for item in value[1:]:
            i = i + 1
            header = ">" + key + "_" + str(i) + "|" + value[0] + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
            sequence = item + "\n" # Replace newlines in sequence, remember to add one to the end.
            writefile.write(header + sequence)