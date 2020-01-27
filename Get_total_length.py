import glob
import re
import csv

read_files = glob.glob("C:/Users/81983/Downloads/ncbi-genomes-2020-01-27-all/*.txt")

with open("C:/Users/81983/Downloads/total_length_comb.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["GenBank assembly accession","RefSeq assembly accession","whether identical","total length"])
    for f in read_files:
        with open(f, "r") as infile:
            info = []
            lines = []
            for line in infile:
                if "# GenBank assembly accession: " in line:
                    edline = line.rstrip('\n')
                    info.append(edline.split(" ")[4])
                if "# RefSeq assembly accession: " in line:
                    edline = line.rstrip('\n')
                    info.append(edline.split(" ")[4])
                if "# RefSeq assembly and GenBank assemblies identical" in line:
                    edline = line.rstrip('\n')
                    info.append(edline.split(" ")[7])
                if "all	all	all	all	total-length" in line:
                    edline = line.rstrip('\n')
                    info.append(edline.split("\t")[5])
            #print(info)
            #print(testline)
            writer.writerow(info)