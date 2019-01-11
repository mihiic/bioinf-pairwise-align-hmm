import os
from shutil import copyfile

for filename in os.listdir('pairwise_alignment_HIV'):
    original = open('pairwise_alignment_HIV/' + filename, 'r')
    parts = []
    for line in original:
        if '>' in line:
            parts.append(line.strip())

    done = 0
    for hmm_filename in os.listdir('pairwise_alignment_HMM'):
        hmm = open('pairwise_alignment_HMM/' + hmm_filename, 'r')
        hmm_parts = []
        for line in hmm:
            if '>' in line:
                hmm_parts.append(line.strip())

        if hmm_parts[0] == parts[0] and hmm_parts[1] == parts[1]:
            copyfile('pairwise_alignment_HMM/' + hmm_filename, 'hmm_sorted/' + filename)

            done += 1
            print(filename + ' -> ' + hmm_filename)

        hmm.close()

    print('found {} matches for {}'.format(done, filename))
    original.close()
