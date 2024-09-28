"""
extract UMIs and reads 

+ve strand:  CS1_UMIF_F00phgback_READ_R00phgbackRC_UMIRRC_CS2RC
-ve strand:  CS2_UMIR_R00phgback_READRC_F00phgbackRC_UMIFRC_CS1RC

"""

###
# ignore pairwise2 deprecated warning
# use biopython 1.81
import warnings
from Bio import BiopythonDeprecationWarning
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
###

from Bio.Seq import Seq
from Bio import SeqIO

# need to use the older method since it provides start and end for all alignments
# whereas the newer one doesn't
from Bio import pairwise2

# missing start and end unforuntately
#from Bio import Align

# try substituion_matrices eventually and evaluate
#from Bio.Align import substitution_matrices

import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import csv
import os.path
import subprocess

# global constants
SEQ_CS1 = 'ACACTGACGACATGGTTCTACA' #ACACTgacgacatggttctaca
SEQ_CS2 = 'TACGGTAGCAGAGACTTGGTCT'
SEQ_CS1_RC = str(Seq(SEQ_CS1).reverse_complement()) #TGTAGAACCATGTCGTCAGTGT
SEQ_CS2_RC = str(Seq(SEQ_CS2).reverse_complement()) #AGACCAAGTCTCTGCTACCGTA


# these are now passed in as args so we can override them
#READ_CONSENSUS_F = 'GGTCTGCTGTTACTGGCGGC'
#     F00phgback = 'GGTCTGCTGTTACTGGCGGC'
#READ_CONSENSUS_R = 'ATGGTGATGATGATGTGCGG'
#     R00phgback = 'ATGGTGATGATGATGTGCGG'

# try this later
# https://community.nanoporetech.com/downloads
# GUPPY_MISMATCH_MATRIX = substitution_matrices.Array("ACGT", 2,
#                                         np.array(
#                                                 [[96.0, -316.0, -192.0, -369.0],
#                                                 [-316.0, 100.0, -352.0, -295.0],
#                                                 [-192.0, -352.0, 98.0, -329.0],
#                                                 [-369.0, -295.0, -329.0, 100.0]]
#                                             )
#                                       )

# TODO: potentially another args??
UMI_LEN = 25
CS_LEN = 22
MATCH_WINDOW = 100

# scores for pairwise2 matching (from porechop: 3,-6,-5,-2)
PAIRWISE2_MATCH = 1
PAIRWISE2_MISMATCH = -1
PAIRWISE2_GAP_OPEN = -0.5
PAIRWISE2_GAP_EXTENSION = -0.1

parser = argparse.ArgumentParser(description='Read UMI Extractor')
parser.add_argument('-i', '--input', type=str, required=True, help='Input fastq file')
parser.add_argument('-o', '--outputdir', type=str, required=True, help='Output directory for UMI binned fastq files')
parser.add_argument('-g1', '--graph1', type=str, required=True, help='Output svg graph file')
parser.add_argument('-g2', '--graph2', type=str, required=True, help='Output svg graph file')
parser.add_argument('-d', action='store_true', required=False, help='Display Graph')
parser.add_argument('-dc', '--discardconsensus', action='store_true', required=False, help='Discard a read if the read consensus does not match')
parser.add_argument('-t', '--table', type=str, required=True, help='Table of UMI read consensus mismatch')
parser.add_argument('-acs', '--alignmentscore', type=float, required=True, help='Alignment score in percentage (ranges from 1 to 99.0) for finding CS1, CS2')
parser.add_argument('-aum', '--alignmentscoreumi', type=float, required=True, help='Alignment score in percentage (ranges from 1 to 99.0) for finding read consensus')
parser.add_argument('-l', '--limit', type=int, default=0, required=False, help='Optional record limit - number of records to be processed')
parser.add_argument('-bbfs', '--backbone-forward', type=str, required=True, help='Backbone forward sequence')
parser.add_argument('-bbrs', '--backbone-reverse', type=str, required=True, help='Backbone reverse sequence')
args = parser.parse_args()

READ_CONSENSUS_F = args.backbone_forward
READ_CONSENSUS_R = args.backbone_reverse
READ_CONSENSUS_R_RC = str(Seq(READ_CONSENSUS_R).reverse_complement())


# newer aligner - but doesn't provide first match offset unlike the older pw2 aligner
# def create_aligner():
#     aligner = Align.PairwiseAligner()
#     aligner.mode = 'global'
#     aligner.match_score = 1
#     aligner.open_gap_score = -0.5
#     aligner.extend_gap_score = -0.1
#     aligner.substitution_matrix = GUPPY_MISMATCH_MATRIX
#     return aligner

def is_cs2_at_the_beginning(seq):
    alignments = pairwise2.align.localms(SEQ_CS2, seq, PAIRWISE2_MATCH, PAIRWISE2_MISMATCH, PAIRWISE2_GAP_OPEN, PAIRWISE2_GAP_EXTENSION)
    # alignments = create_aligner().align(SEQ_CS2, seq)
    if len(alignments) > 0:
        if alignments[0].score / len(SEQ_CS2) > args.alignmentscore:
            return (alignments[0].start < MATCH_WINDOW)
    return False

def beginning_match(cs_type, seq, percent):
    retval = -1
    alignments = pairwise2.align.localms(cs_type, seq, PAIRWISE2_MATCH, PAIRWISE2_MISMATCH, PAIRWISE2_GAP_OPEN, PAIRWISE2_GAP_EXTENSION)
    # alignments = create_aligner().align(cs_type, seq) #grab minimap/guppy ont penalties and use them here ... can be derived through probabilities of diff errors
    if len(alignments) > 0:
        if alignments[0].score / len(cs_type) > percent: #
            retval = alignments[0].end    #see if first alignment is optimal or not... choose the alignment that has the lowest start value... alignments should be a list... use python sort to get lowest a.start
    return retval

def end_match(cs_type, seq, percent):
    retval = -1
    alignments = pairwise2.align.localms(cs_type, seq, PAIRWISE2_MATCH, PAIRWISE2_MISMATCH, PAIRWISE2_GAP_OPEN, PAIRWISE2_GAP_EXTENSION)
    # alignments = create_aligner().align(cs_type, seq)
    if len(alignments) > 0:
        if alignments[0].score / len(cs_type) > percent:
            retval = alignments[0].start    #see if first alignment is optimal or not
    return retval

def append_to_dict(dict, key, value):
    if key not in dict.keys():
        dict[key] = []
    dict[key].append(value)

def print_dict(dict):
    for k, v in dict.items():
        print(k, " => ", len(v))#, " - ", v)

def output_dict_as_fastq(dict, output_folder):
    #output_folder = "Unique_Molecular_Identifier/binning_output_bc3_final"
    umi_counter = 1
    for k, values in dict.items():
        has_file_been_created = False

        id = ""
        seq = ""
        qual = ""
        rc_flag = False

        #print("*** output_dict k is ", k)
        #print(values)

        for param in values:
            id = param[0]
            seq = param[1]
            qual = param[2]
            rc_flag = param[3]

            discardRead = False
            if args.discardconsensus and not rc_flag:
                discardRead = True

            if not discardRead:
                if not has_file_been_created:
                    filename = output_folder + "/umi" + str(umi_counter) + "_" + str(k) + ".fastq"
                    umi_counter += 1
                    print("create new file ", filename)
                    f = open(filename, 'w')
                    has_file_been_created = True

                f.write("@" + id + "\n")
                f.write(str(seq) + "\n")
                f.write("+\n")
                f.write(qual + "\n")

        if has_file_been_created:
            f.close()

def cluster_umi_keys(dict, output_folder):
    # output all keys into a fastq file
    # run with meshclust

    meshclust = "/usr/local/bin/meshclust"

    if os.path.isfile(meshclust):
        meshclust_filename1 = output_folder + "_umis.fa"
        f = open(meshclust_filename1, 'w')
        for k, _ in dict.items():
            if len(k) == UMI_LEN*2 + 1:
                f.write(">" + k + "\n")
                # remove the _ in the middle
                f.write(k[0:UMI_LEN]+k[UMI_LEN+1:len(k)] + "\n")
        f.close()

        meshclust_filename2 = output_folder + "_umis_clustered"

        # assume success
        subprocess.call([meshclust, "-d", meshclust_filename1,
                         "-o", meshclust_filename2, "-t", "0.80"])


        # file format is here
        # https://github.com/HersheyTriesToCode/Identity/blob/master/README.md

        # collect a list of all the cluster centers
        cluster_centers_dict = {}
        with open(meshclust_filename2) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                if len(row) == 4:
                    # print(row[0]) # cluster identifier
                    # print(row[1]) # sequence identifier
                    # print(row[2]) # identity score
                    # print(row[3]) # One of four letters (C, M, E, O)
                    if row[3] == 'C':
                        cluster_id = row[0]
                        # remove > from the beginning
                        umi_center_key = row[1][1:len(row[1])]                        
                        append_to_dict(cluster_centers_dict, cluster_id, umi_center_key)

        # for each cluster center, we use the cluster identifier to know what cluster we are in
        # and we collapse the key mentioned as M to the same bin as C
        with open(meshclust_filename2) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                if len(row) == 4:
                    if row[3] == 'M':
                        cluster_id = row[0]
                        umi_center_key = cluster_centers_dict[cluster_id][0] # returns a list of one item, we want the item
                        umi_key_to_merge = row[1][1:len(row[1])]
                        values = dict[umi_key_to_merge]
                        del dict[umi_key_to_merge]
                        for value in values:
                            append_to_dict(dict, umi_center_key, value)

    else:
        print("Warning: meshclust not found - clustering UMI step skipped")

    return dict


# initialize variables
umi_dict = {}
records_successfully_processed = 0
read_consensus_match = 0 
cs2_beginning = 0

for index, record in enumerate(SeqIO.parse(args.input, "fastq")):

    if args.limit != 0:
        if records_successfully_processed >= args.limit:
            break

    if records_successfully_processed % 10000 == 0:
        print("records successfully processed", records_successfully_processed)
 
    id = record.id
    seq = record.seq
    seq_len = len(seq)

    is_reversed = False

    # check if cs2 is at the beginning, if it is then reverse complement the seq
    # feed it to the cs1 matching section
    if is_cs2_at_the_beginning(seq):
        cs2_beginning += 1
        #print("cs2 found at the beginning - rc it")
        seq = seq.reverse_complement()
        is_reversed = True

    offset = beginning_match(SEQ_CS1, seq, args.alignmentscore)
    if offset > -1 and offset < MATCH_WINDOW:
        # UMI = offset ... offset + len(UMI)
        # READ = offset + len(UMI) ... len(seq)-1-len(CS)? (but we will end up with junk in READ sometimes)
        # to improve READ we need to do a pairwise on cs2rc

        umi = ""
        umif = str(seq)[offset:offset+UMI_LEN]
        umir = ""

        # this is just an approximation - we can ignore this
        #read = str(seq)[offset+UMI_LEN:seq_len-1-CS_LEN]

        # look at the end of the seq - 50 bp
        end_offset = end_match(SEQ_CS2_RC, str(seq)[seq_len-1-MATCH_WINDOW:seq_len-1], args.alignmentscore)

        read_end = -1
        # our window is 50 and 28 is perfect (50-28==22)
        # i.e. 28 == (MATCH_WINDOW-CS_LEN)
        # therefore +-5
        if end_offset > (MATCH_WINDOW-CS_LEN)-5 and end_offset < (MATCH_WINDOW-CS_LEN)+5:
            read_end = seq_len - (MATCH_WINDOW - end_offset) - 1

            # read extraction
            read = str(seq)[offset+UMI_LEN:read_end-UMI_LEN]
            
            umir = str(seq)[read_end-UMI_LEN:read_end]
            umi = umif + '_' + umir

            rc_flag = False
           
            alignments = pairwise2.align.localms(READ_CONSENSUS_F, read[0:len(READ_CONSENSUS_F)], PAIRWISE2_MATCH, PAIRWISE2_MISMATCH, PAIRWISE2_GAP_OPEN, PAIRWISE2_GAP_EXTENSION)
            #print("*** rcf", alignments[0].score / len(READ_CONSENSUS_F) * 100.0)
            #print(pairwise2.format_alignment(*alignments[0]))

            if alignments[0].score / len(READ_CONSENSUS_F) > args.alignmentscoreumi:
                alignments = pairwise2.align.localms(READ_CONSENSUS_R_RC, read[len(read)-len(READ_CONSENSUS_R_RC):len(read)], PAIRWISE2_MATCH, PAIRWISE2_MISMATCH, PAIRWISE2_GAP_OPEN, PAIRWISE2_GAP_EXTENSION)
                #print("*** rcr", alignments[0].score / len(READ_CONSENSUS_R_RC) * 100.0)
                #print(pairwise2.format_alignment(*alignments[0]))

                if alignments[0].score / len(READ_CONSENSUS_R_RC) > args.alignmentscoreumi:
                    read_consensus_match += 1
                    rc_flag = True

            qual_lookup_ascii = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
            qual_string = ""
            qualities = record.letter_annotations["phred_quality"]
            for qual in qualities:
                qual_string += qual_lookup_ascii[qual]

            # we need to modify the qual_string to match
            # how we extracted read from seq
            # that is, reverse if required
            # and use offset
            if is_reversed:
                qual_string = qual_string[::-1]

            # chop the qual string using the same parameters as the read extraction
            qual_string = qual_string[offset+UMI_LEN:read_end-UMI_LEN]

            append_to_dict(umi_dict, umi, [id, read, qual_string, rc_flag])
            records_successfully_processed += 1


# do a mesh cluster on the umi_dict keys and collapse the keys that are clustered
umi_dict = cluster_umi_keys(umi_dict, args.outputdir)

output_dict_as_fastq(umi_dict, args.outputdir)
#print_dict(umi_dict)

# TODO: what about stats for when reject reads due to the alignment score being under the threshold
if records_successfully_processed>0:
    print("read consensus match % ", read_consensus_match/records_successfully_processed * 100)
    print("reads that were RCed % : ", cs2_beginning/records_successfully_processed*100)

max_k = ""
max_v = 0
values = []
for k, v in umi_dict.items():
    value = len(v)
    values.append(value)
    if value > max_v:
        max_v = value
        max_k = k

if max_k != "":
    print("max value is ", max_v)
    print("max umi is ", max_k)

labels = np.arange(0, len(umi_dict.keys()), 1)

# make a histogram 
# figure out a x/y graph
plt.plot(labels, values, marker='o', linestyle='-')

# Adding labels and title
plt.xlabel('UMI number')
plt.ylabel('No. of reads')
plt.title('UMI-read bar graph')
plt.savefig(args.graph1, format='svg')
if args.d:
    plt.show()

with open(args.table, 'w', newline='') as file:
    writer = csv.writer(file)
    field = ["UMI", "Num of reads", "% read consensus mismatch"]
    writer.writerow(field)
    values = []
    if len(umi_dict.items())>0:
        print("percentage of incorrect reads due to read consensus mismatch by bin")
        for k, v in umi_dict.items():
            bin_name = k
            number_of_reads_in_bin = len(v)
            number_of_reads_with_incorrect_read_consensus = 0
            rc_flag = False
            for param in v:
                rc_flag = param[3]
                if not rc_flag:
                    number_of_reads_with_incorrect_read_consensus += 1

            print("bin ", bin_name)
            print("number of reads in bin ", number_of_reads_in_bin)
            print("number of reads with incorrect read consensus ", number_of_reads_with_incorrect_read_consensus)
            percentage_of_read_consensus_mismatches_in_this_bin = number_of_reads_with_incorrect_read_consensus/number_of_reads_in_bin*100.0
            print("percentage of read consensus mismatches in this bin ", percentage_of_read_consensus_mismatches_in_this_bin)
            values.append(percentage_of_read_consensus_mismatches_in_this_bin)
            writer.writerow([bin_name, number_of_reads_in_bin, percentage_of_read_consensus_mismatches_in_this_bin])

labels = np.arange(0, len(umi_dict.keys()), 1)

# make a histogram 
# figure out a x/y graph
plt.clf()
plt.plot(labels, values, marker='o', linestyle='-')

# Adding labels and title
plt.xlabel('UMI number')
plt.ylabel('% of read consensus mismatches')
plt.title('UMI-read consensus bar graph')
plt.savefig(args.graph2, format='svg')
if args.d:
    plt.show()


# cd src/PythonJournals/Unique_Molecular_Identifier/binning_output 
# wc -l * | sort | less 