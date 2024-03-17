# BioTools

The umi_read_extract_pipeline.sh within the pipeline folder can be run in stages 1,2,3,4,5 as `./umi_read_extract_pipeline.sh --stage 1,2,3,4,5`. There are various spice levels at which the pipeline can be run: where acs is alignment threshold for cs1 and cs2 sequences and aum is alignment threshold for the read consensus (phage backbone)

The repo consists of the following python scripts (dependencies can be satisfied by using umi docker file provided):
    
* Stage 1: 00_qual_len_filter.py which taken an input fastq file (processed by guppy for adaptor removal and barcode demultiplexing if dealing with multiple samples) and filters it based on a threshold quality score and a given length range. An output filtered file and log file along with graphs of quality and length distribution pre and post filtering are generated.  

* Stage 2: 01_UMI_read_extract.py extracts a 13bp single UMI situated between CS1 sequence and Read Consensus sequence (i.e. the phage backbone). Reads are reoriented so they all begin with CS1 sequence, extracted UMIs and scFv sequences are stored in a dictionary, and UMI binning/grouping occurs. Output fastq files are made for each UMI within the binning folder and a graph of UMI and frequency of reads is also produced. Output log is also generated.

* Stage 3: 02_UMI_check.py monitors if the extracted UMIs follow the template NNNYRNNNYRNNN pattern. Output log is produced. 

* Stage 4: SPOA or NGSpeciesID can be used to create consensus reads for each UMI fastq file (dependencies can be satisfied by using NGSpeciedID docker file provided). SPOA is quicker and NGSpeciesID is slower but includes polishing with Medaka

* Stage 5: 03_UMI_bintable_enricher.py creates an output csv file with UMI, number of reads within a UMI group, consensus sequence produced by SPOA/NGSpeciesID, corresponding amino acid consensus sequence, and match scores with known clones H8,C12,2F5 (this matching can be commented out to make the pipeline quicker)

    