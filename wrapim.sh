#!/bin/bash

# This script wraps demultiplexing and mapping into a single command.
# In the setup file, the demultiplexed file paths should be given as demux/NAME_R1.fq.gz, 
# where NAME matches the FASTA headers used for demultiplexing.
# See https://github.com/gbedwell/intmap for more information.

set -euo pipefail

R2_FILE=""
BC_R1=""
BC_R2=""
NTHR=1
E=0.1
SINGLE_END=false
THREE_PRIME=false
R1_ONLY=false
R2_ONLY=false
NO_INDEX=false
FILE_MATCH=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -r1|--read1)
            R1_FILE="$2"
            shift 2
            ;;
        -r2|--read2)
            R2_FILE="$2"
            shift 2
            ;;
        -bc_r1|--barcode_r1)
            BC_R1="$2"
            shift 2
            ;;
        -bc_r2|--barcode_r2)
            BC_R2="$2"
            shift 2
            ;;
        -setup|--setup_file)
            SETUP_FILE="$2"
            shift 2
            ;;
        -json|--json_name)
            JSON_NAME="$2"
            shift 2
            ;;
        -out|--out_file)
            OUT_FILE="$2"
            shift 2
            ;;
        -nthr|--nthr)
            NTHR="$2"
            shift 2
            ;;
        -e|--error_rate)
            E="$2"
            shift 2
            ;;
        --r1_only)
            R1_ONLY=true
            shift
            ;;
        --r2_only)
            R2_ONLY=true
            shift
            ;;
        --no_index)
            NO_INDEX=true
            shift
            ;;
        -file_match|--file_match)
            FILE_MATCH="$2"
            shift 2
            ;;
        --three_prime)
            THREE_PRIME=true
            shift
            ;;
        --single_end)
            SINGLE_END=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Required options:"
            echo "  -r1, --read1               Multiplexed R1 FASTQ file"
            echo "  -setup, --setup_file       Analysis setup file"
            echo "  -json, --json_name         Anaysis JSON name"
            echo "  -out, --out_file           Analysis output file"
            echo "  -h, --help                 Show this help message"
            echo "Optional arguments (possibly required in certain circumstances):"
            echo "  -r2, --read2               Multiplexed R2 FASTQ file"
            echo "  -bc_r1, --barcode_r1       R1 barcode file"
            echo "  -bc_r2, --barcode_r2       R2 barcode file"
            echo "  -nthr, --nthr              Number of threads for demultiplexing (default: 1)"
            echo "  -e, --error_rate           Error rate for sequence matching (default: 0.1)"
            echo "  --r1_only                  Use only R1 barcode for demultiplexing"
            echo "  --r2_only                  Use only R2 barcode for demultiplexing"
            echo "  --no_index                 Do not build a barcode index"
            echo "  -file_match, --file_match  2-column csv/tsv for output file naming"
            echo "  --three_prime              Look for barcodes on the 3' end of R1"
            echo "  --single_end               Run mapping in single-end mode"
            echo "  -h, --help                 Show this help message"
            echo "Other mapping options should be defined in the setup file."
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 [OPTIONS]"
    echo "Required options:"
    echo "  -r1, --read1               Multiplexed R1 FASTQ file"
    echo "  -setup, --setup_file       Analysis setup file"
    echo "  -json, --json_name         Anaysis JSON name"
    echo "  -out, --out_file           Analysis output file"
    echo "  -h, --help                 Show this help message"
    echo "Optional arguments (possibly required in certain circumstances):"
    echo "  -r2, --read2               Multiplexed R2 FASTQ file"
    echo "  -bc_r1, --barcode_r1       R1 barcode file"
    echo "  -bc_r2, --barcode_r2       R2 barcode file"
    echo "  -nthr, --nthr              Number of threads for demultiplexing (default: 1)"
    echo "  -e, --error_rate           Error rate for sequence matching (default: 0.1)"
    echo "  --r1_only                  Use only R1 barcode for demultiplexing"
    echo "  --r2_only                  Use only R2 barcode for demultiplexing"
    echo "  --no_index                 Do not build a barcode index"
    echo "  -file_match, --file_match  2-column csv/tsv for output file naming"
    echo "  --three_prime              Look for barcodes on the 3' end of R1"
    echo "  --single_end               Run mapping in single-end mode"
    echo "  -h, --help                 Show this help message"
    echo "Other mapping options should be defined in the setup file."
    exit 0
fi

if [[ "$R1_ONLY" == true ]]; then
    if [[ -z "$R1_FILE" || -z "$BC_R1" || -z "$SETUP_FILE" || -z "$JSON_NAME" || -z "$OUT_FILE" ]]; then
        echo "Error: Some required parameters are missing."
        echo "Use -h or --help for usage information"
        exit 1
    fi
elif [[ "$R2_ONLY" == true ]]; then
    if [[ -z "$R1_FILE" || -z "$R2_FILE" || -z "$BC_R2" || -z "$SETUP_FILE" || -z "$JSON_NAME" || -z "$OUT_FILE" ]]; then
        echo "Error: Some required parameters are missing."
        echo "Use -h or --help for usage information"
        exit 1
    fi
else
    if [[ -z "$R1_FILE" || -z "$R2_FILE" || -z "$BC_R1" || -z "$BC_R2" || -z "$SETUP_FILE" || -z "$JSON_NAME" || -z "$OUT_FILE" ]]; then
        echo "Error: Some required parameters are missing."
        echo "Use -h or --help for usage information"
        exit 1
    fi
fi

echo
echo "Demultiplexing reads..."
echo

DEMUX_CMD="intmap_demux"

DEMUX_CMD+=" -r1 \"$R1_FILE\""

if [[ -n "$R2_FILE" ]]; then
    DEMUX_CMD+=" -r2 \"$R2_FILE\""
fi
if [[ -n "$BC_R1" ]]; then
    DEMUX_CMD+=" -bc_r1 \"$BC_R1\""
fi
if [[ -n "$BC_R2" ]]; then
    DEMUX_CMD+=" -bc_r2 \"$BC_R2\""
fi

DEMUX_CMD+=" -nthr \"$NTHR\""
DEMUX_CMD+=" -e \"$E\""
if [[ "$R1_ONLY" == true ]]; then DEMUX_CMD+=" --r1_only"; fi
if [[ "$R2_ONLY" == true ]]; then DEMUX_CMD+=" --r2_only"; fi
if [[ "$NO_INDEX" == true ]]; then DEMUX_CMD+=" --no_index"; fi
if [[ "$THREE_PRIME" == true ]]; then DEMUX_CMD+=" --three_prime"; fi
if [[ -n "$FILE_MATCH" ]]; then DEMUX_CMD+=" -file_match \"$FILE_MATCH\""; fi

# Run the command
eval $DEMUX_CMD

if [[ $? -eq 0 ]]; then
    echo
    echo "Demultiplexing completed successfully!"
else
    echo "Demultiplexing failed."
    exit 1
fi
echo
echo "Mapping integration sites..."
echo "Writing output to $OUT_FILE"

if [[ "$SINGLE_END" == true ]]; then
    intmap_multi \
        -setup_file "$SETUP_FILE" \
        -json_name "$JSON_NAME" \
        --single_end > "$OUT_FILE"
else
    intmap_multi \
        -setup_file "$SETUP_FILE" \
        -json_name "$JSON_NAME" > "$OUT_FILE"
fi

if [[ $? -eq 0 ]]; then
    echo
    echo "Mapping completed successfully!"
    echo "DONE!"
else
    echo "Mapping failed."
    exit 1
fi