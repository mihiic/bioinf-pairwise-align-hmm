from Bio import SeqIO
import optparse
import os.path

option_parser = optparse.OptionParser()

option_parser.add_option('--mode', action="store", dest="mode",
                         help="Mode can be ungap|toupper", default=None)
option_parser.add_option('--filename', action="store", dest="filename",
                         help="Relative or full path to file containing multiple sequences", default=None)

options, args = option_parser.parse_args()

if options.mode != "ungap" and options.mode != "toupper":
    print("Invalid --mode parameters! Try --help!")
    exit(0)

if not os.path.isfile(options.filename):
    print("File " + options.filename + " not found! Try --help!")
    exit(0)

else:
    print("Loading file " + options.filename + "!")

directory = "parsed_sequences_" + options.filename

if not os.path.exists(directory):
    print("Creating directory " + directory)
    os.mkdir(directory)


if os.name is "posix":
    directory = directory + "/"

if os.name is "nt":
    directory = directory + "\\"


count = 0
for record in SeqIO.parse(options.filename, "fasta"):
    if options.mode == "ungap":
        record.seq = record.seq.ungap('-')
    if options.mode =="toupper":
        record.seq = record.seq.upper()
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
