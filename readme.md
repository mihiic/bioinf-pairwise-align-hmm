# Bioinformatics course - Pairwise alignment using HMM

### Prerequisites

- **python3** (comes native with Linux distributions)
- **biopython** (https://biopython.org/)
- **emboss tools need to be installed on system** (http://emboss.sourceforge.net/)

To install biopython run:

```
sudo pip3 install biopython
```

### Documentation

The **/scripts** folder contains the following:  

- multiple_to_seq.py
- check_alphabet.py
- emboss_needle.py
- get_hmm_params.py
- compare_alignment.py


**multiple_to_seq.py**

Splits a FASTA file containing more sequnces into multiple files. You can choose to remove 
gaps or to convert sequnces to upper case letters. Usage:

```
python3 multiple_align_to_seq.py --mode ungap|toupper --filename <pat>
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

**optimal_align.py**

Given a directory of sequences the script calculates optimal global alignments
between each pair using the needle (Needleman-Wunsch) or stretcher (Rapid Needleman-Wunsch) 
command line tool from EMBOSS (http://emboss.sourceforge.net/). Usage:

```
python3 optimal_align.py --mode needle|stretcher --dir <path_to_directory_containing_seqeunces> --genomename <name_of_the_organism>
```

The results are stored in a directory that the script will create (pairwise_alignment_ + genomename param).

For additional help run:
```
python3 optimal_align.py --help
```

**Files bigger than 1 MB are slowly aligned!**


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

**compare_alignment.py**

Gives comparison of alignments between two files containing alignments or between multiple files in different directories.
```
python3 compare_alignment.py --mode file|dir --path1 <path_to_file_or_dir> --path2 <path_to_file_or_dir>
```

For additional help run:

```
python3 compare_alignment.py --help
```


**hmm_wrapper.py**

Wraps the C++ code that runs the HMM alignment.

Run help for all the necessary parameters:

```
python3 hmm_wrapper.py --help
```


**HMM pairwise alignment**

```
./bioinf <path_to_input_seq_A> <path_to_input_seq_B> <output_filename> <path_to_transition_matrix> <path_to_emission_matrix>
```

All paths must be relative to the **./bioinf** executable. If output filename is inside a folder, make sure that folder
exists before running the program.


**Faculty of Engineering and Computing (FER),  Zagreb 2018./2019. © Matej Jularić, Mihael Međan**

![alt text](https://upload.wikimedia.org/wikipedia/en/2/27/Fakultet_elektrotehnike_i_ra%C4%8Dunarstva%2C_Sveu%C4%8Dili%C5%A1te_u_Zagrebu_%28logo%29.jpg)