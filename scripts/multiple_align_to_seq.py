from Bio import SeqIO
import sys
import os.path


if len(sys.argv) != 2:
    print("################## Specify the path of multiple sequence alignment FASTA file! ##################")

filename = sys.argv[1]

if not os.path.isfile(filename):
    print("################## File " + filename + " not found! ##################")

else:
    print("################## Loading file " + filename + "! ##################")

directory = "parsed_sequences_" + filename

if not os.path.exists(directory):
    print("################## Creating directory " + directory + " ##################")
    os.mkdir(directory)


if os.name is "posix":
    directory = directory + "/"

if os.name is "nt":
    directory = directory + "\\"


for record in SeqIO.parse(filename, "fasta"):
    output_file = open(directory + record.name + ".fasta", "w")
    record.seq = record.seq.ungap('-')
    SeqIO.write(record, output_file, "fasta")
    output_file.close()

print("################## FINISHED! ######################")