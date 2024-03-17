# Transformation for QC (pattern match UMIs to UMI template)
# UMI template = NNNYRNNNYRNNN
# N = A, T, G, C
# Y = C, T
# R = G, A
# look at https://www.w3schools.com/python/ref_string_translate.asp (sam to joe eg)
# look at https://www.w3schools.com/python/python_regex.asp 

import re 
import glob
import argparse

parser = argparse.ArgumentParser(description='UMI Pattern Matcher')
parser.add_argument('-b', '--bindir', type=str, required=True, help='Binning directory of UMI fastq files')
args = parser.parse_args()

def match_umi(umi):
    pattern = re.compile(r'^[ACGT]{3}[CT]{1}[AG]{1}[ACGT]{3}[CT]{1}[AG]{1}[ACGT]{3}$')
    return bool(pattern.match(umi))

umi_template = "NNNYRNNNYRNNN"

#binning_folder = "Unique_Molecular_Identifier/binning_output_3"

list_of_files = glob.glob(args.bindir + r'/*.fastq')
#print(list_of_files)
#print(len(list_of_files))

rejected_count = 0

for filename in list_of_files:

    pattern = re.compile(r'.*umi[0-9]+_([ACGT]+)\.fastq')
    match = pattern.match(filename)
    if match:
        #print(filename)
        #print(match.groups()[0]) 
        umi = match.groups()[0]

        if not match_umi(umi):
            print(f"UMI '{umi}' does not match.")
            rejected_count += 1

total = len(list_of_files)
print("total files examined ", total)
print("rejected ", rejected_count)
print("percentage matched ", 100.0 * (total - rejected_count) / total)



