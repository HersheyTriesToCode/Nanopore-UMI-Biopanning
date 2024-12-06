#!/bin/bash

# this is an example pipeline
# use this as a template, copy and modify to suit your data

# this is how relative paths work... on both mac and windows
# our directory is   blah/users/hershey/data
# our scripts are in blah/users/hershey/src/PythonJournals/UMI
# so the relative path from our directory to scripts is therefore ../src/PythonJournals/UMI

SCRIPT_DIR="../scripts/dual_UMI"
#SCRIPT_DIR="../src/nanopore-umi-biopanning/scripts/dual_UMI"

files=(
#    "../scripts/dual_UMI/input/dual_umi_test/input_single_entry.fastq"
    "../scripts/dual_UMI/input/dual_umi_test/input_double_entry.fastq"
)

OUT_DIR_PREFIX="output/dual_umi_test"

# filter length and quality values
# default qual_score <= 22.5
# for phage library use 820 to 950
# for aptamer library 45 use 162 to 202
FILTER_QUAL_THRESHOLD="-q 22.5"
FILTER_LENGTH_MIN="-lmin 820"
FILTER_LENGTH_MAX="-lmax 950"

# read extract
READ_CONSENSUS_F="-bbfs GGTCTGCTGTTACTGGCGGC" # F00phgback
READ_CONSENSUS_R="-bbrs ATGGTGATGATGATGTGCGG" # R00phgback
#READ_CONSENSUS_F="-bbfs GGAGGCTCTCGGGACGAC" #represents library 45 primer 1
#READ_CONSENSUS_R="-bbrs CTGTAAATCCTAAAGGCGGGACGAC" #represents library 45 primer 2

# bintable enricher
DISABLE_ALIGNMENT=""
#DISABLE_ALIGNMENT="-da" # uncomment to disable


# the higher the stricter the matching / rejection
# -acs Alignment score in percentage (ranges from 1 to 99.0) for finding CS1, CS2
# -aum Alignment score in percentage (ranges from 1 to 99.0) for finding read consensus
# -dc Discard a read if the read consensus does not match
SPICY_HIGH="-acs 0.95 -aum 0.95 -dc"
SPICY_MEDIUM="-acs 0.8 -aum 0.8 -dc"
SPICY_MILD="-acs 0.8 -aum 0.8"
SPICY_LOW="-acs 0.7 -aum 0.7"
SPICY_CUSTOM="-acs 0.x -aum 0.x -dc?"
DESIRED_SPICE=$SPICY_MEDIUM

function display_usage {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -h, --help        Display this help message."
    echo "  --dry-run         Do a dry run of the pipeline do not actually do anything."
    echo "  --stage <number>  Specify the stage numbers (1 to 6)."
    echo "      accepts multiple numbers e.g. 123456 runs all"
    exit 1
}

# Check for the number of arguments
if [ "$#" -eq 0 ]; then
    display_usage
fi

echo running pipeline...
date

# if 1 ie. true then we are only a dry run
DRY_RUN=0

# Process command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            display_usage
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --stage)
            if [ "$#" -gt 1 ]; then
                stage_number=$2
                shift 2
            else
                echo "Error: --stage requires a stage number."
                display_usage
            fi
            ;;
        *)
            echo "Error: Unknown option $1"
            display_usage
            ;;
    esac
done

# Check if stage_number is set
if [ -z "$stage_number" ]; then
    echo "Error: Please specify a stage number using --stage."
    display_usage
fi

for file in "${files[@]}"
do
    INPUT_FASTQ=$file
    filename=`basename $file`
    outdir="${filename%.*}"
    OUT_DIR=$OUT_DIR_PREFIX/$outdir
    BIN_OUT_DIR=$OUT_DIR/binning

    echo input fastq is $INPUT_FASTQ
    echo out dir is $OUT_DIR

    if [ $DRY_RUN -ne 1 ]
    then
        mkdir -p $OUT_DIR
        mkdir -p $BIN_OUT_DIR
    fi

    case $stage_number in
        *1*)
            echo running stage 1...
            if [ $DRY_RUN -ne 1 ]
            then
                if [ ! -e $INPUT_FASTQ ]
                then
                    echo "ERROR!!! ${INPUT_FASTQ} not found!!!"
                    exit -1
                else
                    python $SCRIPT_DIR/00_qual_len_filter.py \
                        -i $INPUT_FASTQ \
                        -o $OUT_DIR/filtered.fastq \
                        -g $OUT_DIR/filtered.svg \
                        $FILTER_QUAL_THRESHOLD \
                        $FILTER_LENGTH_MIN \
                        $FILTER_LENGTH_MAX \
                        2>&1 \
                        > $OUT_DIR/qual_len_filter.log

                    echo wc -l $INPUT_FASTQ
                    echo wc -l $OUT_DIR/filtered.fastq

                    tail $OUT_DIR/qual_len_filter.log
                fi
            fi
            ;;
    esac

    case $stage_number in
        *2*)
            echo running stage 2...

            if [ $DRY_RUN -ne 1 ]
            then
                if [ ! -s $OUT_DIR/filtered.fastq ]
                then
                    echo "ERROR!!! filtered.fastq is empty!!!"
                    exit -1
                else
                    python $SCRIPT_DIR/01_UMI_read_extract.py \
                        -i $OUT_DIR/filtered.fastq \
                        -o $BIN_OUT_DIR \
                        -g1 $OUT_DIR/umi_distribution.svg \
                        -g2 $OUT_DIR/bin_mismatch.svg \
                        -t $OUT_DIR/bin_table.csv \
                        $DESIRED_SPICE \
                        -l 0 \
                        $READ_CONSENSUS_F \
                        $READ_CONSENSUS_R \
                        2>&1 \
                        > $OUT_DIR/UMI_read_extract.log
                fi
            fi
            ;;
    esac

    case $stage_number in
        *3*)
            echo running stage 3...

            if [ $DRY_RUN -ne 1 ]
            then
                if [ ! `ls -1 ${BIN_OUT_DIR}/*.fastq 2>/dev/null` ]
                then
                    echo "ERROR!!! no fastq files in ${BIN_OUT_DIR}!!!"
                    exit -1
                else
                    python $SCRIPT_DIR/02_UMI_check.py -b $BIN_OUT_DIR \
                        2>&1 \
                        > $OUT_DIR/UMI_check.log
                fi
            fi
            ;;
    esac


    case $stage_number in
        *4*)
            echo running stage 4...

            if [ $DRY_RUN -ne 1 ]
            then
                if [ ! `ls -1 ${BIN_OUT_DIR}/*.fastq 2>/dev/null` ]
                then
                    echo "ERROR!!! no fastq files in ${BIN_OUT_DIR}!!!"
                    exit -1
                else
                    # pushd .
                    # cd $BIN_OUT_DIR
                    # for file in *.fastq; do
                    #     consensus_outdir=`basename $file .fastq`
                    #     test $( wc -l < $file ) -gt 4 && \
                    #         NGSpeciesID --ont \
                    #             --consensus --sample_size 500 \
                    #             --m 850 \
                    #             --s 100 \
                    #             --medaka \
                    #             --fastq $file \
                    #             --outfolder $consensus_outdir
                    # done
                    # popd

                    # spoa doco
                    # https://github.com/rvaser/spoa

                    SPOA_CSV_FILE=../bin_table_spoa.csv
                    pushd .
                    cd $BIN_OUT_DIR
                    echo "UMI,consensus sequence" > $SPOA_CSV_FILE
                    for file in *.fastq; do
                        SPOA_CMD="spoa $file -l 0 -r 0 -g -2"
                        UMI=`echo $file | cut -d '_' -f 2,3 | cut -d '.' -f 1`
                        SEQ=`$SPOA_CMD | sed '1d'`
                        echo "$UMI,$SEQ" >> $SPOA_CSV_FILE

                        # output option -r 4 to produce the graph used for the calculations
                        # both as a text output
                        # a dot file
                        # and a resultant svg file
                        # comment out when you are done
                        #echo $file.output4
                        #spoa $file -l 0 -r 4 -g -2 -d $file.dot > $file.output4
                        #dot -Tsvg -O $file.dot

                        # output option -r 1 to produce multiple sequence alignment (FASTA)
                        # we can compare SEQ (the centre i.e. the consensus) letter by letter
                        # against the separate reads in .output1
                        # where - means score -2
                        # where same means +1
                        # where mismatch means score -1
                        echo $file.output1
                        #echo $SEQ
                        spoa $file -l 0 -r 1 -g -2 > $file.output1

                        # both consensus (FASTA) and multiple sequence alignment (FASTA) is outputted
                        echo $file.output2
                        spoa $file -l 0 -r 2 -g -2 > $file.output2
                        

                        # with output1
                    done
                    popd
                fi

            fi
            ;;
    esac

    case $stage_number in
        *5*)
            echo running stage 5...

            if [ $DRY_RUN -ne 1 ]
            then
                if [ ! -e $OUT_DIR/bin_table.csv ]
                then
                    echo "ERROR!!! bin_table.csv not found in ${OUT_DIR}!!!"
                    exit -1
                else
                    python $SCRIPT_DIR/03_UMI_bintable_enricher.py \
                        -o $OUT_DIR \
                        -b $BIN_OUT_DIR \
                        $DISABLE_ALIGNMENT \
                        2>&1 \
                        > $OUT_DIR/UMI_bintable_enricher.log
                fi
            fi
            ;;
    esac


    case $stage_number in
        *6*)
            echo running stage 6...

            if [ $DRY_RUN -ne 1 ]
            then
                if [ ! -e $OUT_DIR/bin_table_enriched.fasta ]
                then
                    echo "ERROR!!! bin_table_enriched.fasta not found in ${OUT_DIR}!!!"
                    exit -1
                else
                    pushd .
                    cd $BIN_OUT_DIR
                    cd ..
                    /usr/local/bin/meshclust -d bin_table_enriched.fasta -o meshclust_output.txt -t 0.95
                    popd
                fi
            fi
            ;;
    esac

done


echo finished pipeline
date
