# look at https://www.w3schools.com/python/ref_string_translate.asp (sam to joe eg)
# look at https://www.w3schools.com/python/python_regex.asp 

import re 
import glob
import argparse

parser = argparse.ArgumentParser(description='UMI Pattern Matcher')
parser.add_argument('-b', '--bindir', type=str, required=True, help='Binning directory of UMI fastq files')
args = parser.parse_args()

# TTTVVVVVTTVVVVVTTVVVVVTTT
# 3, 5,   2,5,   2,5,   3
# V = G, C, A
def match_umi1(umi):
    pattern = re.compile(r'^[T]{3}[GCA]{5}[T]{2}[GCA]{5}[T]{2}[GCA]{5}[T]{3}$')
    return bool(pattern.match(umi))

# AAABBBBBAABBBBBAABBBBBAAA
# 3, 5,   2,5,   2,5,   3
# B = G, C, T
def match_umi2(umi):
    pattern = re.compile(r'^[A]{3}[GCT]{5}[A]{2}[GCT]{5}[A]{2}[GCT]{5}[A]{3}$')
    return bool(pattern.match(umi))

#binning_folder = "Unique_Molecular_Identifier/binning_output_3"

list_of_files = glob.glob(args.bindir + r'/*.fastq')
#print(list_of_files)
#print(len(list_of_files))

rejected_count = 0
total_files_examined = 0

for filename in list_of_files:
    total_files_examined += 1
    pattern = re.compile(r'.*umi[0-9]+_([ACGT]+)_([ACGT]+)\.fastq')
    match = pattern.match(filename)
    if match:
        #print(filename)
        #print(match.groups()[0]) 
        umi1 = match.groups()[0]
        # print("UMI1 ", umi1)
        umi2 = match.groups()[1]
        # print("UMI2 ", umi2)

        if not match_umi1(umi1):
            print(f"UMI '{umi1}' does not match.")
            rejected_count += 1

        if not match_umi2(umi2):
            print(f"UMI '{umi2}' does not match.")
            rejected_count += 1

    else:
        total_files_examined -= 1
        print("warning: " + filename + " is not a UMI path")

total = total_files_examined * 2
print("total UMIs examined ", total)
print("rejected ", rejected_count)
print("percentage matched ", 100.0 * (total - rejected_count) / total)



