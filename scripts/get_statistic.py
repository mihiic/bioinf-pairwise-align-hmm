from Bio import SeqIO
import os


def get_emission_matrix_index(base):
    if base is '-':
        return 0
    if base is 'A':
        return 1
    if base is 'T':
        return 2
    if base is 'C':
        return 3
    if base is 'G':
        return 4


def get_state(x, y):

    #state MATCH
    if (x == 'A' or x == "G" or x == "C" or x == "T") and (y == 'A' or y == "G" or y == "C" or y == "T"):
        return 0

    #state Ix
    if (x == 'A' or x == "G" or x == "C" or x == "T") and y == '-':
        return 1

    #state Iy
    if x == '-' and (y == 'A' or y == "G" or y == "C" or y == "T"):
        return 2


transitions_matrix = [[0, 0, 0],
                      [0, 0, 0],
                      [0, 0, 0]]

emissions_matrix = [[0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0]]

#states
match_state = 0
ix_state = 1
iy_state = 2

#symbols
gap = 0
A = 1
T = 2
C = 3
G = 4


directory = "pairwise_alignment_HIV"


all_files = os.listdir(directory)

for filename in all_files:

    records = []
    length = 0
    for record in SeqIO.parse(directory + "\\" + filename, 'fasta'):
        records.append(record)
        length = len(record.seq)


    last_state = None

    for i in range(length):
        nucleobase1 = records[0].seq[i]
        nucleobase2 = records[1].seq[i]

        l = get_emission_matrix_index(nucleobase1)
        k = get_emission_matrix_index(nucleobase2)

        emissions_matrix[l][k] += 1

        current_state = get_state(nucleobase1, nucleobase2)

        if last_state is None:
            last_state = current_state

        else:
            transitions_matrix[last_state][current_state] += 1
            last_state = current_state

