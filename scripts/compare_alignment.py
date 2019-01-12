from Bio import SeqIO
import optparse
import os

gap_open = -10
gap_extend = -0.5


def score(path):
    records = []

    for record in SeqIO.parse(path, 'fasta'):
        records.append(record)

    length = len(records[0].seq)
    score = 0
    gap = False

    for i in range(length):
        pair = (records[0].seq[i], records[1].seq[i])
        if not gap:
            if '-' in pair:
                gap = True
                score += gap_open
            else:
                if records[0].seq[i] == records[1].seq[i]:
                    score = score + 5
                else:
                    score = score - 4
        else:
            if '-' not in pair:
                gap = False
                if records[0].seq[i] == records[1].seq[i]:
                    score = score + 5
                else:
                    score = score - 4
            else:
                score = score + gap_extend

    return score



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
    result1 = score(options.path1)
    result2 = score(options.path2)
    print("SCORE OPTIMAL: " + str(result1))
    print("SCORE HMM: " + str(result2))

if options.mode == 'dir':
    all_files1 = os.listdir(options.path1)
    all_files2 = os.listdir(options.path2)
    '''
    if len(all_files1) != len(all_files2):
        print("ERROR! Directories don't contain the same number of file!")
        exit(0)
    '''
    length = len(all_files2)
    result = 0
    diff = 0

    res_optimal = []
    res_hmm = []

    for i in range(length):
        result1 = score(options.path1 + "/" + all_files1[i])
        res_optimal.append(result1)
        result2 = score(options.path2 + "/" + all_files2[i])
        res_hmm.append(result2)
        diff = diff + (result1 - result2)
    print("MAX OPTIMAL ALIGNMENT SCORE: " + str(max(res_optimal)))
    print("MIN OPTIMAL ALIGNMENT SCORE: " + str(min(res_optimal)))
    print("MAX HMM SCORE: " + str(max(res_hmm)))
    print("MIN HMM SCORE: " + str(min(res_hmm)))
    print("Average difference in alignment scores between OPTIMAL AND HMM: " + str(diff / length))
