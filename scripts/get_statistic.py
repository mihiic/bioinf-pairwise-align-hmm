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
        return 1

    #state Ix
    if (x == 'A' or x == "G" or x == "C" or x == "T") and y == '-':
        return 2

    #state Iy
    if x == '-' and (y == 'A' or y == "G" or y == "C" or y == "T"):
        return 3


transitions_matrix = [[0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0]]

emissions_matrix = [[0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0]]


#states
begin_state = 0
match_state = 1
ix_state = 2
iy_state = 3
end_state = 4

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
            transitions_matrix[begin_state][current_state] += 1
            last_state = current_state

        else:
            transitions_matrix[last_state][current_state] += 1
            last_state = current_state
    transitions_matrix[last_state][end_state] += 1

sum_begin = 0
sum_match = 0
sum_ix = 0
sum_iy = 0

for i in range(5):
    sum_begin += transitions_matrix[begin_state][i]
    sum_match += transitions_matrix[match_state][i]
    sum_ix += transitions_matrix[ix_state][i]
    sum_iy += transitions_matrix[iy_state][i]

for i in range(5):
    transitions_matrix[begin_state][i] /= sum_begin
    transitions_matrix[match_state][i] /= sum_match
    transitions_matrix[ix_state][i] /= sum_ix
    transitions_matrix[iy_state][i] /= sum_iy

transition_matrix_file = open("transition_matrix.txt", "w")

for i in range(5):
    to_file = ""
    for j in range(5):
        if j == 4:
            to_file = to_file + " " + str(transitions_matrix[i][j])
            continue
        to_file = to_file + str(transitions_matrix[i][j]) + ","
    to_file = to_file + "\n"
    transition_matrix_file.write(to_file)

transition_matrix_file.close()


sum_emission_match = 0

for i in range(1, 5):
    for j in range(1, 5):
        sum_emission_match += emissions_matrix[i][j]

for i in range(1, 5):
    for j in range(1, 5):
        emissions_matrix[i][j] /= sum_emission_match

sum_emission_ix = 0
sum_emission_iy = 0

for i in range(5):
    sum_emission_ix += emissions_matrix[0][i]
    sum_emission_iy += emissions_matrix[i][0]

for i in range(5):
    emissions_matrix[0][i] /= sum_emission_ix
    emissions_matrix[i][0] /= sum_emission_iy

emission_matrix_file = open("emission_matrix.txt", "w")
for i in range(5):
    to_file = ""
    for j in range(5):
        if j == 4:
            to_file = to_file + " " + str(emissions_matrix[i][j])
            continue
        to_file = to_file + str(emissions_matrix[i][j]) + ","
    to_file = to_file + "\n"
    emission_matrix_file.write(to_file)

emission_matrix_file.close()

