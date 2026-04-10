#!/bin/bash

# This script aids in preparing the setup file for intmap_multi.
# It takes a tab-separated file containing sample names, LTR sequences, and linker sequences
# (by default columns 1, 2, and 3, respectively) and returns the following as comma-separated values:
# 1. Sample names
# 2. Demultiplexed file names (including path prefix)
# 3. ltr and linker sequences (both 20 nucleotides by default)
# 4. ltr- and linker-end barcode fasta files for demultiplexing.
#
# The script treats the 5' ends of the given LTR and linker sequences as barcode regions and
# the 3' ends as search regions. It therefore assumes that the entire sequenced portion
# of the LTR and linker sequences are provided (i.e., everything after the sequencing adapters).

NAME_COL=1
LTR_COL=2
LINKER_COL=3
INDEX_COL1=0
INDEX_COL2=0
PATH_PREFIX="demux/"
LTR_LEN=20
LINKER_LEN=20
BC_LEN=10
LTRN=0
LINKN=0
FILE_EXT="fq.gz"
SE=false
MPLEX=""
NO_BC=false

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
        -index_col1|--index_col1)
            INDEX_COL1="$2"
            shift 2
            ;;
        -index_col2|--index_col2)
            INDEX_COL2="$2"
            shift 2
            ;;
        -path_prefix|--path_prefix)
            PATH_PREFIX="$2"
            shift 2
            ;;
        -ltr_len|--ltr_len)
            LTR_LEN="$2"
            shift 2
            ;;
        -linker_len|--linker_len)
            LINKER_LEN="$2"
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
        --no_bc)
            NO_BC=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 <sample_info.txt> [OPTIONS]"
            echo
            echo "Required arguments:"
            echo "  <sample_info.txt>                 Tab-separated file with sample info (names, LTR, linker)"
            echo "Optional arguments:"
            echo "  -name_col, --name_col N           Name column (default: 1)"
            echo "  -ltr_col, --ltr_col N             LTR column (default: 2)"
            echo "  -linker_col, --linker_col N       Linker column (default: 3)"
            echo "  -index_col1, --index_col1 N       I1 index sequence column (default: 0, disabled)
  -index_col2, --index_col2 N       I2 index sequence column (default: 0, disabled)"
            echo "  -path_prefix, --path_prefix STR   Path prefix for output files (default: demux/)"
            echo "  -ltr_len, --ltr_len N             LTR search sequence length (default: 20)"
            echo "  -linker_len, --linker_len N       Linker search sequence length (default: 20)"
            echo "  -bc_len, --bc_len N               Barcode length (default: 10)"
            echo "  -file_ext, --file_ext STR         File extension (default: fq.gz)"
            echo "  --se                              Single-end data (default: false)"
            echo "  -ltrn, --ltrn N                   Replace left-most N bases of ltr_bc with N (default: 0)"
            echo "  -linkn, --linkn N                 Replace left-most N bases of linker_bc with N (default: 0)"
            echo "  -multiplex_prefix, --multiplex_prefix STR  Prefix for barcode fasta files"
            echo "  --no_bc                           Skip writing barcode fasta files"
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
    echo "  -ltr_len N         LTR length (default: 10)"
    echo "  -linker_len N      Linker length (default: 10)"
    echo "  -bc_len N          Barcode length (default: 10)"
    echo "  -file_ext STR      File extension (default: fq.gz)"
    echo "  --se               Single-end data (default: false)"
    echo "  -ltrn N            Replace left-most N bases of ltr_bc with N (default: 0)"
    echo "  -linkn N           Replace left-most N bases of linker_bc with N (default: 0)"
    exit 1
fi

if (( LTR_LEN < 5 )); then
    echo "ERROR: --ltr_len must be >= 5.)" >&2
    exit 1
fi
if (( LINKER_LEN < 5 )); then
    echo "ERROR: --linker_len must be >= 5.)" >&2
    exit 1
fi
if (( BC_LEN < 5 )); then
    echo "ERROR: --bc_len must be >= 5." >&2
    exit 1
fi

names=()
ltr_seqs=()
linker_seqs=()
index1_seqs=()
index2_seqs=()
counter=0
while IFS=$'\t' read -r -a cols || [ -n "${cols[*]}" ]; do
    ((counter++))
    name_clean="${cols[$((NAME_COL-1))]// /}"
    ltr_clean="${cols[$((LTR_COL-1))]// /}"
    linker_clean="${cols[$((LINKER_COL-1))]// /}"

    name_clean="${name_clean#"${name_clean%%[![:space:]]*}"}"
    name_clean="${name_clean%"${name_clean##*[![:space:]]}"}"
    ltr_clean="${ltr_clean#"${ltr_clean%%[![:space:]]*}"}"
    ltr_clean="${ltr_clean%"${ltr_clean##*[![:space:]]}"}"
    linker_clean="${linker_clean#"${linker_clean%%[![:space:]]*}"}"
    linker_clean="${linker_clean%"${linker_clean##*[![:space:]]}"}"

    names+=("$name_clean")
    ltr_seqs+=("$ltr_clean")
    linker_seqs+=("$linker_clean")

    if (( INDEX_COL1 > 0 )); then
        index1_clean="${cols[$((INDEX_COL1-1))]// /}"
        index1_clean="${index1_clean#"${index1_clean%%[![:space:]]*}"}"
        index1_clean="${index1_clean%"${index1_clean##*[![:space:]]}"}"
        index1_seqs+=("$index1_clean")
    fi

    if (( INDEX_COL2 > 0 )); then
        index2_clean="${cols[$((INDEX_COL2-1))]// /}"
        index2_clean="${index2_clean#"${index2_clean%%[![:space:]]*}"}"
        index2_clean="${index2_clean%"${index2_clean##*[![:space:]]}"}"
        index2_seqs+=("$index2_clean")
    fi
done < "$TSV_FILE"

if [[ ${#names[@]} -ne $counter ]]; then
    echo "Length of names (${#names[@]}) does not match counter ($counter)." >&2
    exit 1
fi

echo "Sample number:"
echo "${counter}"

echo "Sample names:"
echo "${names[*]}" | tr ' ' ','

echo "Output files:"
Rltr_files=()
Rlink_files=()
for n in "${names[@]}"; do
    Rltr_files+=("${PATH_PREFIX}${n}_Rltr.${FILE_EXT}")
        if [[ "$SE" != true ]]; then
            Rlink_files+=("${PATH_PREFIX}${n}_Rlink.${FILE_EXT}")
    fi
done
echo "${Rltr_files[*]}" | tr ' ' ','
if [[ "$SE" != true ]]; then
    echo "${Rlink_files[*]}" | tr ' ' ','
fi

extract_segments() {
    local seq="$1" len_seq="$2" len_bc="$3"
    local l=${#seq}
    local len_check=$(( len_seq + len_bc ))
    if (( l < len_check )); then
        echo "ERROR: Sequence '$seq' is < ${len_check}" >&2
        return 1
    fi
    local seg="${seq: -$len_seq}"
    echo "$seg"
}

BC_LEN_CHECK=$( [[ "$NO_BC" == true ]] && echo 0 || echo "$BC_LEN" )

echo "ltr:"
ltr_arr=()
for s in "${ltr_seqs[@]}"; do
    ltr_arr+=( "$(extract_segments "$s" "$LTR_LEN" "$BC_LEN_CHECK")" )
done
echo "${ltr_arr[*]}" | tr ' ' ','

echo "linker:"
linker_arr=()
for s in "${linker_seqs[@]}"; do
    linker_arr+=( "$(extract_segments "$s" "$LINKER_LEN" "$BC_LEN_CHECK")" )
done
echo "${linker_arr[*]}" | tr ' ' ','

if [[ "$NO_BC" != true ]]; then
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

    if (( INDEX_COL1 > 0 )); then
        if [ -z "$MPLEX" ]; then
            index1_file="index1_bc.fa"
        else
            index1_file="${MPLEX}_index1_bc.fa"
        fi
        echo "Writing ${index1_file}..."
        : > "$index1_file"
        for i in "${!names[@]}"; do
            echo ">${names[$i]}" >> "$index1_file"
            echo "${index1_seqs[$i]}" >> "$index1_file"
        done
    fi

    if (( INDEX_COL2 > 0 )); then
        if [ -z "$MPLEX" ]; then
            index2_file="index2_bc.fa"
        else
            index2_file="${MPLEX}_index2_bc.fa"
        fi
        echo "Writing ${index2_file}..."
        : > "$index2_file"
        for i in "${!names[@]}"; do
            echo ">${names[$i]}" >> "$index2_file"
            echo "${index2_seqs[$i]}" >> "$index2_file"
        done
    fi
fi

echo "Done."