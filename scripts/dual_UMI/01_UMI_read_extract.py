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

# for umi extraction
import re

import difflib

# global constants
SEQ_CS1 = 'ACACTGACGACATGGTTCTACA'
SEQ_CS2 = 'TACGGTAGCAGAGACTTGGTCT'
SEQ_CS1_RC = str(Seq(SEQ_CS1).reverse_complement()) #TGTAGAACCATGTCGTCAGTGT
SEQ_CS2_RC = str(Seq(SEQ_CS2).reverse_complement()) #AGACCAAGTCTCTGCTACCGTA

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

UMI_LEN = 25

# scores for pairwise2 matching
PAIRWISE2_MATCH = 3
PAIRWISE2_MISMATCH = -6
PAIRWISE2_GAP_OPEN = -5
PAIRWISE2_GAP_EXTENSION = -2

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

# returns -1 if there is no match
# otherwise returns the start and end position of the match
def match(searchseq, seq, percent, debug=False):
    retval = [-1,-1]
    alignments = pairwise2.align.localms(searchseq, seq, PAIRWISE2_MATCH, PAIRWISE2_MISMATCH, PAIRWISE2_GAP_OPEN, PAIRWISE2_GAP_EXTENSION)
    if len(alignments) > 0:
        if debug:
            print("searchseq: " + searchseq)
            for i, alignment in enumerate(alignments):
                print(str(i) + " score: " + str(alignment.score) + " start:" + str(alignment.start) + " end:" + str(alignment.end) + " length: " + str(alignment.end-alignment.start))
                print(alignment)
        if alignments[0].score / len(searchseq) > percent:
            retval = [alignments[0].start,alignments[0].end]
    return retval

def fuzzy_search(needle, haystack):
    position = -1
    highest_index = -1
    highest_score = -1
    needle_length = len(needle)

    if needle_length <= len(haystack):
        for i in range(len(haystack)-needle_length):
            current_window = haystack[i:needle_length+i]
            sm = difflib.SequenceMatcher(None, needle, current_window)
            score = sm.ratio()
            if score > highest_score:
                highest_score = score
                highest_index = i

                # debug
                # print()
                # listOfOffsetsSubSeqLabels = []
                # listOfOffsetsSubSeqLabels.append((i, needle, "needle with score " + str(score)))
                # renderOffsetListCascade(haystack, listOfOffsetsSubSeqLabels)
                # print()

    if highest_score >= 0.95:
        position = highest_index

    return position

def found_cs2(seq):
    # this alignment match appears to produce large amounts of false positives
    # given the current parameters for the pairwise alignment 
    #return match(SEQ_CS2, seq, args.alignmentscore) != [-1,-1]

    cs2_string = str(SEQ_CS2)
    seq_string = str(seq)

    # exact string match
    #return seq_string.find(cs2_string) > -1

    # fuzzy match
    return fuzzy_search(cs2_string, seq_string) > -1

# seq is the seq which has offset from
# list is a list of [(offset, subseq, label),...]
# assumes the offsets are monotonicly increasing
def renderOffsetListCascade(seq, list):
    print(seq)
    for _, item in enumerate(list):
        (offset,subseq,label)=item
        print(' '*offset + subseq)
        print(' '*offset + label)

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
            f.write(">" + k + "\n")
            # remove the _ in the middle
            f.write(k.replace('_','') + "\n")
        f.close()

        meshclust_filename2 = output_folder + "_umis_clustered"

        ret = subprocess.run([meshclust, "-d", meshclust_filename1,
                         "-o", meshclust_filename2, "-t", "0.90"], capture_output=True)
        
        if len(ret.stderr) > 0:
            print("meshclust failed")
            print("error: ", ret.stderr)
            # ignore meshclust failure
            return dict
        
        if not os.path.isfile(meshclust_filename2):
            print("meshclust failed")
            print("error: output file not created ", meshclust_filename2)
            # ignore meshclust failure
            return dict

        # file format is here
        # https://github.com/HersheyTriesToCode/Identity/blob/master/README.md

        # collect a list of all the cluster centers
        print("cluster_umi_keys first pass - collect cluster centers")
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

        print("number of cluster centers ", len(cluster_centers_dict))

        # for each cluster center, we use the cluster identifier to know what cluster we are in
        # and we collapse the key mentioned as M to the same bin as C
        print("cluster_umi_keys second pass - merge bins")
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
                        print("center ", umi_center_key, " <- merging with key ", umi_key_to_merge)
                        for value in values:
                            append_to_dict(dict, umi_center_key, value)

    else:
        print("Warning: meshclust not found - clustering UMI step skipped")

    return dict


# initialize variables
umi_dict = {}
read_consensus_match = 0 
number_of_reversed_reads = 0
number_of_records = 0

umi_search_pattern = '.*(TTT[AGC]{3,7}TT[AGC]{3,7}TT[AGC]{3,7}TTT).*'

# create head and tail subsequences to use alignment to find offset into the sequence
UMI_EG = "TTTAAAAATTAAAAATTAAAAATTT"
head = SEQ_CS1 + UMI_EG + READ_CONSENSUS_F
tail = READ_CONSENSUS_R_RC + UMI_EG + SEQ_CS2_RC

for _, record in enumerate(SeqIO.parse(args.input, "fastq")):

    number_of_records += 1

    if args.limit != 0:
        if number_of_records >= args.limit:
            break

    if number_of_records % 10000 == 0:
        print("records processed", number_of_records)
 
    id = record.id
    seq = record.seq
    seq_len = len(seq)

    is_reverse_complement = False

    # if we have a cs2 then we know the seq needs to be reverse complemented to be the forward version
    if found_cs2(seq):
        number_of_reversed_reads += 1
        seq = seq.reverse_complement()
        is_reverse_complement = True

    print(seq)
    print(len(seq))
    
    listOfOffsetsSubSeqLabels = []
    
    start_of_read = -1
    end_of_read = -1
    read = ""
    
    umi = ""
    umif = ""
    umir = ""

    haveHead = False
    haveTail = False
    haveUmif = False
    haveUmir = False
    haveRead = False
    
    positions = match(head, seq, args.alignmentscore)   
    if positions != [-1,-1]:
        extracted_head = str(seq[positions[0]:positions[1]])
        haveHead = True
        listOfOffsetsSubSeqLabels.append((positions[0], extracted_head, "head"))
        start_of_read = positions[1]
        
        # regex find umif in extracted_head
        match_result = re.search(umi_search_pattern, extracted_head)
        if match_result:
            umif = match_result.group(1)
            haveUmif = True
            listOfOffsetsSubSeqLabels.append((positions[0]+extracted_head.find(umif), umif, "umif"))
        
    positions = match(tail, seq, args.alignmentscore)   
    if positions != [-1,-1]:
        extracted_tail = str(seq[positions[0]:positions[1]])
        haveTail = True
        listOfOffsetsSubSeqLabels.append((positions[0], extracted_tail, "tail"))
        
        if haveHead:
            end_of_read = positions[0]
            read = seq[start_of_read:end_of_read]
            haveRead = True
            listOfOffsetsSubSeqLabels.append((start_of_read, read, "read"))

        # regex find umir in extracted_tail
        match_result = re.search(umi_search_pattern, extracted_tail)
        if match_result:
            umir = match_result.group(1)
            haveUmir = True
            listOfOffsetsSubSeqLabels.append((positions[0]+extracted_tail.find(umir), umir, "umir"))

    if haveUmif and haveUmir:
        umi = umif + "_" + umir
    
    # troubleshoot the subseqs
    renderOffsetListCascade(seq, listOfOffsetsSubSeqLabels)

    if haveRead and haveUmif and haveUmir:
        print("complete record found")
    
        # rc_flag = False
        
        # alignments = pairwise2.align.localms(READ_CONSENSUS_F, read[0:len(READ_CONSENSUS_F)], PAIRWISE2_MATCH, PAIRWISE2_MISMATCH, PAIRWISE2_GAP_OPEN, PAIRWISE2_GAP_EXTENSION)
        # #print("*** rcf", alignments[0].score / len(READ_CONSENSUS_F) * 100.0)
        # #print(pairwise2.format_alignment(*alignments[0]))

        # if alignments[0].score / len(READ_CONSENSUS_F) > args.alignmentscoreumi:
        #     alignments = pairwise2.align.localms(READ_CONSENSUS_R_RC, read[len(read)-len(READ_CONSENSUS_R_RC):len(read)], PAIRWISE2_MATCH, PAIRWISE2_MISMATCH, PAIRWISE2_GAP_OPEN, PAIRWISE2_GAP_EXTENSION)
        #     #print("*** rcr", alignments[0].score / len(READ_CONSENSUS_R_RC) * 100.0)
        #     #print(pairwise2.format_alignment(*alignments[0]))

        #     if alignments[0].score / len(READ_CONSENSUS_R_RC) > args.alignmentscoreumi:
        #         read_consensus_match += 1
        #         rc_flag = True
        rc_flag = True # disable rc check for now

        # follow the same conventions as dorado fastq format
        # https://github.com/nanoporetech/dorado/blob/c3a2952356e2506ef1de73b0c9e14784ab9a974a/dorado/utils/fastq_reader.h#L20C68-L20C80
        qual_lookup_ascii = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
        qual_string = ""
        qualities = record.letter_annotations["phred_quality"]

        for qual in qualities:
            qual_string += qual_lookup_ascii[qual]

        # we need to modify the qual_string to match
        # how we extracted read from seq
        # that is, reverse if required
        # and use offset
        if is_reverse_complement:
            qual_string = qual_string[::-1]

        # chop the qual string using the same parameters as the read extraction
        qual_string = qual_string[start_of_read:end_of_read]

        append_to_dict(umi_dict, umi, [id, read, qual_string, rc_flag])


# do a mesh cluster on the umi_dict keys and collapse the keys that are clustered
# comment out this line to disable umi key clustering
umi_dict = cluster_umi_keys(umi_dict, args.outputdir)

output_dict_as_fastq(umi_dict, args.outputdir)
#print_dict(umi_dict)

# TODO: what about stats for when reject reads due to the alignment score being under the threshold
if number_of_records>0:
    print("number of records ", number_of_records)
    #disabled print("read consensus match % ", read_consensus_match/number_of_records * 100.0)
    print("reads that were RCed % : ", number_of_reversed_reads/number_of_records * 100.0)

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