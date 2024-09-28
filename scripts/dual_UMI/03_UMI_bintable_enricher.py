###
# ignore pairwise2 deprecated warning
# use biopython 1.81
import warnings
from Bio import BiopythonDeprecationWarning
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
###

import os
import sys
import csv
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

parser = argparse.ArgumentParser(description='UMI BinTable Enricher')
parser.add_argument('-o', '--outputdir', type=str, required=True, help='Output directory for UMI binned fastq files')
parser.add_argument('-b', '--bindir', type=str, required=True, help='Binning directory of UMI fastq files')
parser.add_argument('-da', action='store_true', required=False, help='Disable Alignment')
args = parser.parse_args()

BINTABLE_FILENAME = "bin_table.csv"
NEW_BINTABLE_FILENAME = "bin_table_enriched.csv"
SPOA_FILENAME = "bin_table_spoa.csv"
OUTPUT_DIR = args.outputdir
BIN_OUTPUT_DIR = args.bindir
BINTABLE_PATH = OUTPUT_DIR + "/" + BINTABLE_FILENAME
NEW_BINTABLE_PATH = OUTPUT_DIR + "/" + NEW_BINTABLE_FILENAME
SPOA_TABLE_PATH = OUTPUT_DIR + "/" + SPOA_FILENAME

# these are trimmed sequeneces
SEQ_H8 = "GGTCTGCTGTTACTGGCGGCCCAGCCGGCCATGGCCCAGGTGCAGCTGTTGGAGTCTGGAGCAGAGGTGAAAAAGCCCGGGGAGTCTCTGAAGATCTCCTGTAAGGGTTCTGGATACAGCTTTACCAGCTACTGGATCGGCTGGGTGCGCCAGATGCCCGGGAAAGGCCTGGAGTGGATGGGGATCATCTATCCTGGTGACTCTGATACCAGATACAGCCCGGCCTTCCAAGGCCAGGTCACCATCTCAGCCGACAAGTCCATCAGCACCGCCTACCTGCAGTGGAGCAGCCTGAAGGCCTCGGACACCGCCATGTATTACTGTGCGAGACGGGGGATTTTTGGAGTGGAAAATCTTGATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAGGTGGAGGCGGTTCAGGCGGAGGTGGCTCTGGCGGTGGCGCTGGCCAGTCTGCCCTGACTCAGCCTCGCTCAGTGTCCGGGTCTCCTGGACAGTCAGTCACCATCTCCTGCACTGGAACCAGCAGTGATGTTGGTGGTTATAACTATGTCTCCTGGTACCAACAGCACCCAGGCAAAGCCCCCAAACTCATGATTTATGATGTCAGTAAGCGGCCCTCAGGGGTCCCTGATCGCTTCTCTGGCTCCAAGTCTGGCAGCACGGCCTCCCTGACAATCTCTGGGCTCCAGGCTGAGGACGAGGCTGAATATTACTGCAGCTCATATACAACCAGCGGCACTTATGTCTTCGGAACTGGGACCAAGCTGACCGTCCTAGGTGCGGCCGCACATCATCATCACCAT"
SEQ_C12 = "GGTCTGCTGTTACTGGCGGCCCAGCCGGCCATGGCCCAGGTCCAGCTGGTACAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAGGCTTTCCTGCAAGACTTCTGGATACAACTTCACTAGTTATGCTATGCATTGGGTGCGCCGGGCCCCCGGACAAAGGCTTGAATGGATGGGATGGATCAACGCTGGCAATGGTAAGACAGAATATTCACCGGGGTTTCAGGGCAGAGTCACCATTACCACAGACACATCCGCGAGCACAGGCTTCATGGAACTGAGCAGCCTGAGATCTGAAGACACGGCTATGTATTACTGTGCGAGAGATGGCTTGGGTGGTCGCGCCTTCAACGGAATGGACGTCTGGGGCCACGGCACCCTGGTCACCGTCTCCTCAGGTGGAGGCGGTTCAGGCGGAGGTGGCTCTGGCGGTGGCGCTAGCTCCTATGAGCTGACACAGCCACCCTCGGTGTCAGTGTCCCCAGGACAGACGGCCAGGATCACCTGCTCTGGAGATGCATTGCCAAAGCAATATGCTTATTGGTACCAGCAGAAGCCAGGCCAGGCCCCTGTGCTGGTGATATATAAAGACAGTGAGAGGCCCTCAGGGATCCCTGAGCGATTCTCTGGCTCCAGCTCAGGGACAACAGTCACGTTGACCATCAGTGGAGTCCAGGCAGAAGACGAGGCTGACTATTACTGTCAATCAGCAGACAGCAGTGGTACTTGGGTGTTCGGCGGAGGGACCAAGCTGACCGTCCTAGGTGCGGCCGCACATCATCATCACCAT"
SEQ_2F5 = "GGTCTGCTGTTACTGGCGGCCCAGCCGGCCATGGCCGAGGTTCGCCTGCAACAGCCTGGGACTGAACTGGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGACTTCTGGCTACACCTTCACCAGGTACTGGATGCACTGGGTGAGGCAGAGGCCTGGACAAGGCCTTGAGTGGATTGGAAATATTCTTCCTTACGATGGTGGTACTAACTACAATGAGAGGTTCAAGAACAAGGCCACACTGACTGTAGACAGATCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCCCCCTACTATAGTGGGGACTTTGACTCCTGGGGCCAAGGCACCACTCTCACAGTCTCCTCGGGTGGTGGTGGTTCTGGCGGCGGCGGCTCCGGTGGAGGTGGATCCGATGTCCAGATGACACAGACTACATCCTCCCTGTCTGCCTCTCTGGGAGACAGAGTCACCATCAGTTGCAGGGCAAGTCAGGACATTAGCAATTATTTAAACTGGTATCAGCAGAAACCAGATGGAACTGTTAAACTCCTGATCTACTACACATCAAGATTACACTCAGGAATCCCATCAAGGTTCAGTGGCAGTGGGTCTGGAACAGATTATTCTCTCACCATTAGCAACCTGGAGCAAGAAGATGTTGCCACTTACTTTTGCCAACAGGGTAATACGCTTCCGTGGACGTTCGGTGGAGGCACCAAGCTGGAAATGAAACGCGCGGCCGCACATCATCATCACCAT"

# scores for pairwise2 matching (from porechop: 3,-6,-5,-2)
PAIRWISE2_MATCH = 1
PAIRWISE2_MISMATCH = -1
PAIRWISE2_GAP_OPEN = -0.5
PAIRWISE2_GAP_EXTENSION = -0.1
# PAIRWISE2_MATCH = 1
# PAIRWISE2_MISMATCH = -1
# PAIRWISE2_GAP_OPEN = -0.5
# PAIRWISE2_GAP_EXTENSION = -0.1

# if the spoa table exists, we read the whole thing
# into memory as a dictionary for fast lookups
spoa_dict = {}
if os.path.exists(SPOA_TABLE_PATH):
    with open(SPOA_TABLE_PATH, newline='') as input_spoa_table_file:
        reader = csv.reader(input_spoa_table_file)
        i = 0
        for row in reader:
            if i == 0:
                print() # ignore header row
            else:
                umi = row[0]
                consensus = row[1]
                spoa_dict[umi] = consensus

            i += 1

# open the input csv for enrichment
with open(BINTABLE_PATH, newline='') as input_csvfile:
    with open(NEW_BINTABLE_PATH, 'w', newline='') as output_csvfile:
        writer = csv.writer(output_csvfile)

        reader = csv.reader(input_csvfile)
        i = 0
        for row in reader:
            if i == 0:
                # expecting row is the header row
                if not len(row) == 3:
                    sys.stderr.write("Error: the bintable file had the wrong number of fields, expected 3 got " + str(len(row)))
                    sys.exit(1)
                else:
                    if args.da == False:
                        field = ["UMI", "Num of reads", "% read consensus mismatch", "consensus sequence", "consensus amino acid", "H8 match score", "C12 match score", "2F5 match score"]
                    else:
                        field = ["UMI", "Num of reads", "% read consensus mismatch", "consensus sequence", "consensus amino acid"]
                    writer.writerow(field)
            else:
                print(row) # debug

                umi = row[0]
                numOfReads = row[1]
                percentReadConsensusMismatch = row[2]
                consensus_seq = ""
                consensus_aa = ""
                seq_H8_match = 0.0
                seq_C12_match = 0.0
                seq_2F5_match = 0.0
               
                # if we have the spoa table in memory
                # then use it in preference
                if umi in spoa_dict:
                    consensus_seq = spoa_dict[umi]
                else:
                    # failback to running a medaka file
                    # we want to read the consensus seq if it exists
                    # binning/umi*_[UMI]/consensus_reference_0.fasta
                    consensus_reference_files = glob.glob(BIN_OUTPUT_DIR + "/" + r'umi*_' + umi + "/consensus_reference_0.fasta")
                    if len(consensus_reference_files) == 0:
                        # if we have no matches then we did not run a consensus program
                        x = BIN_OUTPUT_DIR + "/" + r'umi*_' + umi + ".fastq"
                        UMIs_with_one_read = glob.glob(x)
                        if len(UMIs_with_one_read) == 1:
                            for record in SeqIO.parse(UMIs_with_one_read[0], "fastq"):
                                consensus_seq = record.seq 
                                break

                        elif len(UMIs_with_one_read) == 0:
                            print("read probably discarded")
                            # file doesn't exist so... we don't have any seq to read
                            consensus_seq = ""
                        else:
                            sys.stderr.write("Error: UMI fastq file has more than one copy" + str(len(UMIs_with_one_read)))
                            sys.exit(1)
                    elif len(consensus_reference_files) == 1:
                        for record in SeqIO.parse(consensus_reference_files[0], "fasta"):
                            consensus_seq = record.seq
                            break
                    else:
                        sys.stderr.write("Error: umi binning dir returned more than 1 result, ", consensus_reference_files)
                        sys.exit(1)

                if len(consensus_seq) > 0:
                    consensus_aa = Seq(str(consensus_seq)).translate()

                    if args.da == False:
                        alignments = pairwise2.align.globalms(consensus_seq, SEQ_H8, PAIRWISE2_MATCH, PAIRWISE2_MISMATCH, PAIRWISE2_GAP_OPEN, PAIRWISE2_GAP_EXTENSION)
                        if len(alignments) > 0:
                            seq_H8_match = alignments[0].score / max(len(consensus_seq), len(SEQ_H8))

                        alignments = pairwise2.align.globalms(consensus_seq, SEQ_C12, PAIRWISE2_MATCH, PAIRWISE2_MISMATCH, PAIRWISE2_GAP_OPEN, PAIRWISE2_GAP_EXTENSION)
                        if len(alignments) > 0:
                            seq_C12_match = alignments[0].score / max(len(consensus_seq), len(SEQ_C12))

                        alignments = pairwise2.align.globalms(consensus_seq, SEQ_2F5, PAIRWISE2_MATCH, PAIRWISE2_MISMATCH, PAIRWISE2_GAP_OPEN, PAIRWISE2_GAP_EXTENSION)
                        if len(alignments) > 0:
                            seq_2F5_match = alignments[0].score / max(len(consensus_seq), len(SEQ_2F5))

                        writer.writerow([umi, numOfReads, percentReadConsensusMismatch, consensus_seq, consensus_aa, seq_H8_match, seq_C12_match, seq_2F5_match])
                    else:
                        writer.writerow([umi, numOfReads, percentReadConsensusMismatch, consensus_seq, consensus_aa])


            # debug - only process 3 rows from the input csv file
            #if i > 1000:
                #break
            i += 1

            