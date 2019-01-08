from Bio import SeqIO
import optparse
import os


def is_valid(nucleotide):
    alphabet = ['A', 'G', 'C', 'T']
    for a in alphabet:
        if nucleotide == a:
            return True
    return False


option_parser = optparse.OptionParser()

option_parser.add_option('--dir', action="store", dest="dir",
                         help="path to directory containing sequences", default=None)
option_parser.add_option('--delete', action="store", dest="delete",
                         help="Set to true if you want to delete sequences that don't have allowed alphabet",
                         default="false")

options, args = option_parser.parse_args()

if options.dir is None:
    print("Specify --dir parameter!")
    exit(0)



if options.delete == "true":
    print("Delete option is set to true. All the files that don't have the allowed alphabet will be deleted.")
    choice = input("Are you sure you want to proceed? [y/n]: ")

    if choice != "y":
        print("Exiting script check_alphabet.py!")
        exit(0)


if not os.path.exists(options.dir):
    print("Specified directory " + options.dir + " does not exist!")
    exit(0)


all_files = os.listdir(options.dir)

for file in all_files:
    rec = None
    for record in SeqIO.parse(options.dir + "/" + file, 'fasta'):
        rec = record
    record_length = len(rec.seq)
    valid_alphabet = True
    for i in range(record_length):
        if not is_valid(rec.seq[i]):
            valid_alphabet = False
            break

    if not valid_alphabet:
        print("The file " + file + " does not contain the allowed alphabet!")
        if options.delete == "true":
            print("Removing file " + file + " !")
            os.remove(options.dir + "/" + file)






