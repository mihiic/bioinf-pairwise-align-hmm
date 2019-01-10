from Bio import SeqIO
import optparse
import os


def compare(path1, path2):
    records1 = []
    records2 = []
    for record in SeqIO.parse(path1, 'fasta'):
        records1.append(record)

    for record in SeqIO.parse(path2, 'fasta'):
        records2.append(record)

    length1 = len(records1[0].seq)
    length2 = len(records2[0].seq)

    length = length1
    rest = 0

    if length1 > length2:
        length = length2
        rest = length1 - length2

    if length1 < length2:
        length = length1
        rest = length2 - length1

    match = 0

    for i in range(length):
        if records1[0].seq[i] == records2[0].seq[i] and records1[1].seq[i] == records2[1].seq[i]:
            match = match + 1

    return (match / (length + rest)) * 100


option_parser = optparse.OptionParser()

option_parser.add_option('--mode', action="store", dest="mode",
                         help="file|dir", default=None)
option_parser.add_option('--path1', action="store", dest="path1",
                         help="path to directory or file", default=None)
option_parser.add_option('--path2', action="store", dest="path2",
                         help="path to directory or file", default=None)


options, args = option_parser.parse_args()

if options.mode is None or (options.mode != 'file' and options.mode != 'dir'):
    print("INVALID PARAMETERS SUPPLIED! Try --help!")
    exit(0)

if options.path1 is None or options.path2 is None:
    print("INVALID PATHS! Try --help!")
    exit(0)

if not os.path.exists(options.path1) and not os.path.exists(options.path2):
    print("SUPPLY CORRECT PATHS! Use --path1 and --path2! Try --help!")
    exit(0)


if options.mode == 'file':
    result = compare(options.path1, options.path2)
    print("COMPARISON RESULT: " + str(result))

if options.mode == 'dir':
    all_files1 = os.listdir(options.path1)
    all_files2 = os.listdir(options.path2)
    if len(all_files1) != len(all_files2):
        print("ERROR! Directories don't contain the same number of file!")
        exit(0)

    length = len(all_files1)
    result = 0
    for i in range(length):
        result = result + compare(options.path1 + '/' + all_files1[i], options.path2 + '/' + all_files2[i])

    print("COMPARISON RESULT: " + str(result / length))
