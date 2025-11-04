#!/bin/bash

# This script aids in preparing the setup file for intmap_multi.
# It takes a tab-separated file containing sample names, LTR sequences, and linker sequences
# (by default columns 1, 2, and 3, respectively) and returns the following as comma-separated values:
# 1. Sample names
# 2. Demultiplexed file names (including path prefix)
# 3. ltr5, ltr3, linker5, linker3 (all 10 nucleotides by default)
# 4. ltr- and linker-end barcode fasta files.
#
# The script treats the 5' ends of the given LTR and linker sequences as barcode regions and
# the 3' ends as search regions. It therefore assumes that the entire sequenced portion
# of the LTR and linker sequences are provided (i.e., everything after the sequencing adapters).

NAME_COL=1
LTR_COL=2
LINKER_COL=3
PATH_PREFIX="demux/"
LTR5_LEN=10
LTR3_LEN=10
LINKER5_LEN=10
LINKER3_LEN=10
BC_LEN=10
LTRN=0
LINKN=0
FILE_EXT="fq.gz"
SE=false
MPLEX=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -name_col|--name_col)
            NAME_COL="$2"
            shift 2
            ;;
        -ltr_col|--ltr_col)
            LTR_COL="$2"
            shift 2
            ;;
        -linker_col|--linker_col)
            LINKER_COL="$2"
            shift 2
            ;;
        -path_prefix|--path_prefix)
            PATH_PREFIX="$2"
            shift 2
            ;;
        -ltr5_len|--ltr5_len)
            LTR5_LEN="$2"
            shift 2
            ;;
        -ltr3_len|--ltr3_len)
            LTR3_LEN="$2"
            shift 2
            ;;
        -linker5_len|--linker5_len)
            LINKER5_LEN="$2"
            shift 2
            ;;
        -linker3_len|--linker3_len)
            LINKER3_LEN="$2"
            shift 2
            ;;
        -bc_len|--bc_len)
            BC_LEN="$2"
            shift 2
            ;;
        -ltrn|--ltrn)
            LTRN="$2"
            shift 2
            ;;
        -linkn|--linkn)
            LINKN="$2"
            shift 2
            ;;
        -file_ext|--file_ext)
            FILE_EXT="$2"
            shift 2
            ;;
        --se)
            SE=true
            shift
            ;;
        -multiplex_prefix|--multiplex_prefix)
            MPLEX="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 <sample_info.txt> [OPTIONS]"
            echo
            echo "Required arguments:"
            echo "  <sample_info.txt>           Tab-separated file with sample info (names, LTR, linker)"
            echo "Optional arguments:"
            echo "  -name_col, --name_col N           Name column (default: 1)"
            echo "  -ltr_col, --ltr_col N             LTR column (default: 2)"
            echo "  -linker_col, --linker_col N       Linker column (default: 3)"
            echo "  -path_prefix, --path_prefix STR   Path prefix for output files (default: demux/)"
            echo "  -ltr5_len, --ltr5_len N           LTR5 barcode length (default: 10)"
            echo "  -ltr3_len, --ltr3_len N           LTR3 search region length (default: 10)"
            echo "  -linker5_len, --linker5_len N     Linker5 barcode length (default: 10)"
            echo "  -linker3_len, --linker3_len N     Linker3 search region length (default: 10)"
            echo "  -bc_len, --bc_len N               Barcode length (default: 10)"
            echo "  -file_ext, --file_ext STR         File extension (default: fq.gz)"
            echo "  --se                      Single-end data (default: false)"
            echo "  -ltrn, --ltrn N                   Replace left-most N bases of ltr_bc with N (default: 0)"
            echo "  -linkn, --linkn N                 Replace left-most N bases of linker_bc with N (default: 0)"
            echo "  -multiplex_prefix, --multiplex_prefix STR  Prefix for barcode fasta files"
            echo "  -h, --help                        Show this help message"
            echo
            echo "This script prepares setup files for intmap_multi. See https://github.com/gbedwell/intmap for more information."
            exit 0
            ;;
        *)
            TSV_FILE="$1"
            shift
            ;;
    esac
done

if [ -z "${TSV_FILE:-}" ]; then
    echo "Usage: $0 <sample_info.txt> [options]"
    echo "Options:"
    echo "  -name_col N        Name column (default: 1)"
    echo "  -ltr_col N         LTR column (default: 2)"
    echo "  -linker_col N      Linker column (default: 3)"
    echo "  -path_prefix STR   Path prefix (default: demux/)"
    echo "  -ltr5_len N        LTR5 length (default: 10)"
    echo "  -ltr3_len N        LTR3 length (default: 10)"
    echo "  -linker5_len N     Linker5 length (default: 10)"
    echo "  -linker3_len N     Linker3 length (default: 10)"
    echo "  -bc_len N          Barcode length (default: 10)"
    echo "  -file_ext STR      File extension (default: fq.gz)"
    echo "  --se               Single-end data (default: false)"
    echo "  -ltrn N            Replace left-most N bases of ltr_bc with N (default: 0)"
    echo "  -linkn N           Replace left-most N bases of linker_bc with N (default: 0)"
    exit 1
fi

names=()
ltr_seqs=()
linker_seqs=()

while IFS=$'\t' read -r -a cols || [ -n "${cols[*]}" ]; do
    # Remove all whitespace from sequences using parameter expansion
    name_clean="${cols[$((NAME_COL-1))]// /}"
    ltr_clean="${cols[$((LTR_COL-1))]// /}"
    linker_clean="${cols[$((LINKER_COL-1))]// /}"
    
    # Also trim leading/trailing whitespace
    name_clean="${name_clean#"${name_clean%%[![:space:]]*}"}"
    name_clean="${name_clean%"${name_clean##*[![:space:]]}"}"
    ltr_clean="${ltr_clean#"${ltr_clean%%[![:space:]]*}"}"
    ltr_clean="${ltr_clean%"${ltr_clean##*[![:space:]]}"}"
    linker_clean="${linker_clean#"${linker_clean%%[![:space:]]*}"}"
    linker_clean="${linker_clean%"${linker_clean##*[![:space:]]}"}"
    
    names+=("$name_clean")
    ltr_seqs+=("$ltr_clean")
    linker_seqs+=("$linker_clean")
done < "$TSV_FILE"

echo "Sample names:"
echo "${names[*]}" | tr ' ' ','

echo "Output files:"
R1_files=()
R2_files=()
for n in "${names[@]}"; do
    R1_files+=("${PATH_PREFIX}${n}_R1.${FILE_EXT}")
        if [[ "$SE" != true ]]; then
            R2_files+=("${PATH_PREFIX}${n}_R2.${FILE_EXT}")
    fi
done
echo "${R1_files[*]}" | tr ' ' ','
if [[ "$SE" != true ]]; then
    echo "${R2_files[*]}" | tr ' ' ','
fi

extract_segments() {
    local seq="$1" len5="$2" len3="$3"
    local l=${#seq}
    if (( l < len3 + 5 )); then
        echo "ERROR: Sequence '$seq' too short for 3' extraction" >&2
        return 1
    fi
    seg3="${seq: -$len3}"
    seg5="${seq: -$((len3+len5)):$len5}"
    if (( ${#seg5} < 5 )); then
        echo "ERROR: 5' segment too short in '$seq'" >&2
        return 1
    fi
    echo "$seg5 $seg3"
}

echo "ltr5:"
ltr5_arr=()
for s in "${ltr_seqs[@]}"; do
    ltr5_arr+=( "$(extract_segments "$s" "$LTR5_LEN" "$LTR3_LEN" | awk '{print $1}')" )
done
echo "${ltr5_arr[*]}" | tr ' ' ','

echo "ltr3:"
ltr3_arr=()
for s in "${ltr_seqs[@]}"; do
    ltr3_arr+=( "$(extract_segments "$s" "$LTR5_LEN" "$LTR3_LEN" | awk '{print $2}')" )
done
echo "${ltr3_arr[*]}" | tr ' ' ','

echo "linker5:"
linker5_arr=()
for s in "${linker_seqs[@]}"; do
    linker5_arr+=( "$(extract_segments "$s" "$LINKER5_LEN" "$LINKER3_LEN" | awk '{print $1}')" )
done
echo "${linker5_arr[*]}" | tr ' ' ','

echo "linker3:"
linker3_arr=()
for s in "${linker_seqs[@]}"; do
    linker3_arr+=( "$(extract_segments "$s" "$LINKER5_LEN" "$LINKER3_LEN" | awk '{print $2}')" )
done
echo "${linker3_arr[*]}" | tr ' ' ','

if [ -z "$MPLEX" ]; then
    ltr_file="ltr_bc.fa"
    linker_file="linker_bc.fa"
else
    ltr_file="${MPLEX}_ltr_bc.fa"
    linker_file="${MPLEX}_linker_bc.fa"
fi

echo "Writing ${ltr_file} and ${linker_file}..."

: > "$ltr_file"
: > "$linker_file"
for i in "${!names[@]}"; do
    ltr_bc="${ltr_seqs[$i]:0:$BC_LEN}"
    linker_bc="${linker_seqs[$i]:0:$BC_LEN}"
    if (( LTRN > 0 )); then
        ltr_bc="$(printf '%*s' "$LTRN" | tr ' ' 'N')${ltr_bc:$LTRN}"
    fi
    if (( LINKN > 0 )); then
        linker_bc="$(printf '%*s' "$LINKN" | tr ' ' 'N')${linker_bc:$LINKN}"
    fi
    echo ">${names[$i]}" >> "$ltr_file"
    echo "$ltr_bc" >> "$ltr_file"
    echo ">${names[$i]}" >> "$linker_file"
    echo "$linker_bc" >> "$linker_file"
done

echo "Done."