# Bioinformatics course - Pairwise alignment using HMM

###Prerequisites

- **python3** (comes native with Linux distributions)
- **biopython** (https://biopython.org/)

To install biopython run:

```
sudo pip3 install biopython
```

### Documentation

The **/scripts** folder contains the following:  

- multiple_align_to_seq.py
- check_alphabet.py
- emboss_needle.py
- get_hmm_params.py


**multiple_align_to_seq.py**

Converts a multiple alignment fasta file to original sequences. Usage:

```
python3 multiple_align_to_seq.py <path_to_file>
```

Path can be relative to the script file or absolute.
The script will automatically remove sequences that have an alphabet different
from A, G, C and T. The results are stored in a directory that the script creates.


**check_alphabet.py**

Checks if all the sequences in a directory have the alphabet A, G, C and T. 
Usage:

```
python3 check_alphabet.py --dir <path_to_directory_containing_sequences> --delete true|false
```

For additional help run:

```
python3 check_alphabet.py --help
```

**emboss_needle.py**

Given a directory of sequences the script calculates optimal global alignments
between each pair using the online tool https://www.ebi.ac.uk/Tools/psa/emboss_needle/nucleotide.html

```
python3 emboss_needle.py --email <your_email> --dir <path_to_directory_containing_seqeunces> --genomename <name_of_the_organism>
```

The results are stored in a directory that the script will create (pairwise_alignment_ + genomename param).

For additional help run:
```
python3 emboss_needle.py --help
```

The emboss_needle.py script **won't accept files bigger than 1 MB**.


**get_hmm_params.py**

Calculates the emission and transition matrix for HMM given a directory with pairwise aligned sequences.
Usage:
```
python3 get_hmm_params.py --dir <path_to_directory_containing_pairwise_alginments>
```

The results are stored in emision_matrix.txt and transition_matrix.txt

For additional help run:
```
python3 get_hmm_params.py --help
```


**HMM pairwise alignment**

```
./bioinf <path_to_input_seq_A> <path_to_input_seq_B> <output_filename> <path_to_transition_matrix> <path_to_emission_matrix>
```

All paths must be relative to the **./bioinf** executable. If output filename is inside a folder, make sure that folder
exists before running the program.


**Faculty of Engineering and Computing (FER),  Zagreb 2018./2019. © Matej Jularić, Mihael Međan**
