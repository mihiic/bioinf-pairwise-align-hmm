import optparse
import os
import subprocess
import time

option_parser = optparse.OptionParser()

option_parser.add_option('--dir', action="store", dest="dir",
                         help="Directory containing sequences you want to align. Must be relative to C++ executable", default=None)
option_parser.add_option('--emission', action="store", dest="emission",
                         help="Path to emission matrix file. Must be relative to C++ executable!", default=None)
option_parser.add_option('--transition', action="store", dest="transition",
                         help="Path to emission matrix file. Must be relative to C++ executable!", default=None)
option_parser.add_option('--exe', action="store", dest="exe",
                         help="Full path to C++ executable!", default=None)


options, args = option_parser.parse_args()


if options.dir is None or options.emission is None or options.transition is None:
    print("INVALID PARAMETERS SUPPLIED! Try --help!")
    exit(0)

if not os.path.exists(options.exe):
    print("Specify full path to C++ executable!")
    exit(0)

all_sequences = os.listdir(options.dir)

directory_to_store = "scripts/pairwise_alignment_HMM_random10000"

if not os.path.exists(directory_to_store):
    print("Creating directory " + directory_to_store)
    os.mkdir(directory_to_store)

sum = 0
cnt = 0

for i in range(len(all_sequences)):
    seq_A = options.dir + "/" + all_sequences[i]
    for j in range(i + 1, len(all_sequences)):
        outfile = directory_to_store + "/" + "alignment_" + str(i) + "_and_" + str(j) + ".fasta"
        seq_B = options.dir + "/" + all_sequences[j]

        # zadnji argument je target sample size
        t1 = time.time()
        subprocess.call([options.exe, seq_A, seq_B, outfile, options.transition, options.emission, '250'])
        t2 = time.time()

        delta = t2 - t1
        sum = sum + delta
        cnt = cnt + 1
        print("Execution time in [s]: " + str(delta))
        print("Execution time in [min]: " + str(delta / 60))
        print("Finished: " + outfile)

    print("DONE FOR: " + str(i))

print("Average execution time in [s]: " + str(sum / cnt))
print("Average execution time in [min]: " + str(sum / (cnt * 60)))
