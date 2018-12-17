from Bio import SeqIO
import sys
import os.path


if len(sys.argv) != 2:
    print("################## Specify the path of multiple sequence alignment FASTA file! ##################")
    exit(0)

filename = sys.argv[1]

if not os.path.isfile(filename):
    print("################## File " + filename + " not found! ##################")
    exit(0)

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


count = 0
for record in SeqIO.parse(filename, "fasta"):
    record.seq = record.seq.ungap('-')
    valid_alphabet = True
    for i in range(len(record.seq)):
        if record.seq[i] != 'A' and record.seq[i] != 'G' and record.seq[i] != 'T' and record.seq[i] != 'C':
            valid_alphabet = False
            break
    if valid_alphabet:
        count = count + 1
        output_file = open(directory + record.name + ".fasta", "w")
        SeqIO.write(record, output_file, "fasta")
        output_file.close()

print("################## FINISHED! ######################")
