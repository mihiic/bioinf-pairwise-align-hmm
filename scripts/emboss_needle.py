from urllib.request import urlopen
import urllib
import optparse
import requests
import os

option_parser = optparse.OptionParser()

option_parser.add_option('--email', action="store", dest="email",
                         help="email to send notification when done", default=None)
option_parser.add_option('--matrix', action="store", dest="matrix", default="EDNAFULL")
option_parser.add_option('--gapopen', action="store", dest="gapopen", default="10")
option_parser.add_option('--gapext', action="store", dest="gapext", default="0.5")
option_parser.add_option('--endweight', action="store", dest="endweight", default="false")
option_parser.add_option('--endopen', action="store", dest="endopen", default="10")
option_parser.add_option('--endextend', action="store", dest="endextend", default="0.5")
option_parser.add_option('--dir', action="store", dest="dir",
                         help="path to directory containg sequences", default=None)
option_parser.add_option('--genomename', action="store", dest="genomename",
                         help="name of the genome", default=None)


options, args = option_parser.parse_args()


if options.email is None:
    print("Email must be specified!")
    exit(0)

if options.dir is None:
    print("Directory must be specified!")
    exit(0)

base_url = "https://www.ebi.ac.uk/Tools/services/rest/emboss_needle"

params = {'email': options.email, 'title': options.genomename, 'matrix': options.matrix, 'gapopen': options.gapopen,
          'gapext': options.gapext, 'endweight': options.endweight, 'endopen': options.endopen,
          'endextend': options.endextend, 'format': "fasta", 'stype': "dna"}


def submit_request(parameters):
    request_data = urllib.parse.urlencode(parameters)
    job_submit_request = urlopen(base_url + "/run", bytes(request_data, 'utf-8'))
    jobID_bytes = job_submit_request.read()
    return jobID_bytes.decode()


def get_job_status(jobID):
    response = requests.get(base_url + '/status/' + jobID)
    if response.text == "FINISHED":
        return True
    return False


def get_result(jobID, filename):
    response = requests.get(base_url + '/result/' + jobID + '/aln')
    file = open(filename, "w")
    file.write(response.text)
    file.close()


directory_to_store = "pairwise_alignment_" + params['title']

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


for i in range(len(all_sequences)):
    fa = open(sequences_directory + all_sequences[i], 'r')
    asequence = fa.read()
    fa.close()
    print("Started for " + str(i))
    for j in range(i+1, len(all_sequences)):
        fb = open(sequences_directory + all_sequences[j], 'r')
        bsequence = fb.read()
        fb.close()
        params['asequence'] = asequence
        params['bsequence'] = bsequence

        job_ID = submit_request(params)
        while not get_job_status(job_ID):
            pass

        get_result(job_ID, directory_to_store + "alignment_" + str(i) + "_and_" + str(j) + ".fasta")

    print("Finished for " + str(i))

print("DONE!")











