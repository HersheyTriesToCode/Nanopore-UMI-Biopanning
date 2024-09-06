# Purpose
To extract dual UMIs from start and end.

# Design
```
+ve strand:  CS1_UMIF_F00phgback_READ_R00phgbackRC_UMIRRC_CS2RC
-ve strand:  CS2_UMIR_R00phgback_READRC_F00phgbackRC_UMIFRC_CS1RC

UMIF - forward UMI
UMIR - reverse UMI
CS1 - conserved sequence 1
CS2 - conserved sequence 2
F00phgback - FORWARD PHAGE BACKBONE
R00phgback - REVERSE PHAGE BACKBONE
RC - reverse complement
```

# Install
1. install latest miniconda, use base environment is simplest
  * https://docs.anaconda.com/free/miniconda/index.html
    * at the time of writing:
      * `conda 24.3.0`
      * `Python 3.12.2`
2. use pip to install from docker/umi/requirements.txt
  * `pip install -r ../docker/umi/requirements.txt`
3. upgrade miniconda as needed later
  * https://docs.anaconda.com/free/anaconda/install/update-version/
    * `conda update -n base conda` 

# Tests

## Test using input_single_entry.fastq
```shell
# using dual_UMI/input/dual_umi_test by default
# all test output is in dual_UMI/output

cd dual_UMI

# modify dual_umi_pipeline.sh to use
# INPUT_FASTQ="input/dual_umi_test/input_single_entry.fastq"

# filter using test data - expect one read in the filtered.fastq
./dual_umi_pipeline.sh --stage 1

# single binning to binning/umi1_TTTAAAAATTAAAAATTAAAAATTT_AAAGGGGGAAGGGGGAAGGGGGAAA.fastq
./dual_umi_pipeline.sh --stage 2
```

## Test using input_double_entry.fastq

```shell
cd dual_UMI

# modify dual_umi_pipeline.sh to use
# INPUT_FASTQ="input/dual_umi_test/input_double_entry.fastq"

# filter using test data - expect two reads in the filtered.fastq
./dual_umi_pipeline.sh --stage 1

# single binning to binning/umi1_TTTAAAAATTAAAAATTAAAAATTT_AAAGGGGGAAGGGGGAAGGGGGAAA.fastq
# with both reads being identical
./dual_umi_pipeline.sh --stage 2

# check UMI_check.log is
# total UMIs examined  2
# rejected  0
# percentage matched  100.0
./dual_umi_pipeline.sh --stage 3

```
