import re

with open('D:/Codes/korem_16s/Data/rrnDB/rrnDB-5.6_16S_rRNA.fasta') as fa_file:
    ssp_list = []
    for line in fa_file:
        if line.startswith(">"):
            fsp_line = line[1:-1].split("|")
            #print(fsp_line)
            pos_line = re.split('\.\.|\s', fsp_line[-1])
            ssp_line = fsp_line[:-1]
            ssp_line.extend(pos_line)
            #print(ssp_line)
            ssp_list.append(ssp_line)
#print(ssp_list[0:4])

import csv
with open('D:/Codes/korem_16s/Data/rrnDB/rrnDB-5.6_16s_position.csv',"w") as csvfile:
    # columns name
    writer = csv.writer(csvfile)
    # positions and other information
    writer.writerow(["organism_name","record_id","RefSeq_sequence","chromosome","start","end","strand"])
    writer.writerows(ssp_list)

'''with open('D:/Codes/korem_16s/Data/rrnDB/rrnDB-5.6_16s_NCBIaa.txt',"w") as txtfile:
    for record in ssp_list:
        txtfile.write(record[1]+'\n')'''
