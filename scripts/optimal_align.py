from Bio.Emboss.Applications import NeedleCommandline, StretcherCommandline
import optparse
import os
import time

option_parser = optparse.OptionParser()

option_parser.add_option('--mode', action="store", dest="mode",
                         help="Choose optimal alignment tool to use needle|stretcher", default=None)
option_parser.add_option('--dir', action="store", dest="dir",
                         help="path to directory containg sequences", default=None)
option_parser.add_option('--genomename', action="store", dest="genomename",
                         help="name of the genome", default=None)


options, args = option_parser.parse_args()

if options.dir is None:
    print("Directory must be specified!")
    exit(0)

directory_to_store = "pairwise_alignment_" + options.genomename

if not os.path.exists(directory_to_store):
    print("Creating directory " + directory_to_store)
    os.mkdir(directory_to_store)


if os.name is "posix":
    directory_to_store = directory_to_store + "/"
    sequences_directory = options.dir + "/"
if os.name is "nt":
    directory_to_store = directory_to_store + "\\"
    sequences_directory = options.dir + "\\"

all_sequences = os.listdir(options.dir)

sum = 0
cnt = 0

for i in range(len(all_sequences)):
    asequence = sequences_directory + all_sequences[i]
    print("Started for " + str(i))
    for j in range(i+1, len(all_sequences)):
        bsequence = sequences_directory + all_sequences[j]
        outfile = directory_to_store + "alignment_" + str(i) + "_and_" + str(j) + ".fasta"

        t1 = time.time()
        if options.mode == "needle":
            needle_cli = NeedleCommandline(asequence=asequence, bsequence=bsequence,
                                           gapopen=10.0, gapextend=0.5, outfile=outfile, aformat="fasta")
            stdout, stderr = needle_cli()

        if options.mode == "stretcher":
            stretcher_cli = StretcherCommandline(asequence=asequence, bsequence=bsequence,
                                           gapopen=16.0, gapextend=4.0, outfile=outfile, aformat="fasta")
            stdout, stderr = stretcher_cli()
        t2 = time.time()
        delta = t2 - t1
        sum = sum + delta
        cnt = cnt + 1
        print("Execution time in [s]: " + str(delta))
        print("Execution time in [min]: " + str(delta / 60))

    print("Finished for " + str(i))

print("Average execution time in [s]: " + str(sum / cnt))
print("Average execution time in [min]: " + str(sum / (cnt * 60)))
print("DONE!")
