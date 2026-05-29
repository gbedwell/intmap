#!/bin/bash

# This script replaces primer IDs in a tab-separated input file with their
# corresponding primer sequences, as defined in a tab-separated primer reference file.
# By default, the reference file is assumed to contain primer IDs in column 1 and 
# primer sequences in column 2. Similarly, the input file is assumed to contain
# sample names in column 1 and primer IDs for integrant-end and linker-end primers
# in columns 2 and 3, respectively. The output file mirrors this structure, but
# replaces primer IDs with the actual primer sequence.

INPUT=""
OUTPUT=""
PRIMER_REF=""
PRIMER_ID_COL=1
PRIMER_SEQ_COL=2
SAMPLE_COL=1
P1_COL=2
P2_COL=3

while [[ $# -gt 0 ]]; do
    case $1 in
        -input|-i)
            INPUT="$2"
            shift 2
            ;;
        -output|-o)
            OUTPUT="$2"
            shift 2
            ;;
        -primer_ref|-r)
            PRIMER_REF="$2"
            shift 2
            ;;
        -primer_id_col|--primer_id_col)
            PRIMER_ID_COL="$2"
            shift 2
            ;;
        -primer_seq_col|--primer_seq_col)
            PRIMER_SEQ_COL="$2"
            shift 2
            ;;
        -sample_col|--sample_col)
            SAMPLE_COL="$2"
            shift 2
            ;;
        -p1_col|--p1_col)
            P1_COL="$2"
            shift 2
            ;;
        -p2_col|--p2_col)
            P2_COL="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo
            echo "Required arguments:"
            echo "  -input, -i STR              Input TSV file"
            echo "  -output, -o STR             Output TSV file"
            echo "Optional arguments:"
            echo "  -primer_ref, -r STR         Primer reference file"
            echo "  -primer_id_col N            Primer ID column in reference file (default: 1)"
            echo "  -primer_seq_col N           Primer sequence column in reference file (default: 2)"
            echo "  -sample_col N               Sample name column in input file (default: 1)"
            echo "  -p1_col N                   First primer ID column in input file (default: 2)"
            echo "  -p2_col N                   Second primer ID column in input file (default: 3)"
            echo "  -h, --help                  Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$INPUT" ]]; then
    echo "ERROR: Input file is required (-input/-i)." >&2
    exit 1
fi

if [[ -z "$OUTPUT" ]]; then
    echo "ERROR: Output file is required (-output/-o)." >&2
    exit 1
fi

if [[ -z "$PRIMER_REF" ]]; then
    echo "ERROR: Primer reference file is required (-primer_ref/-r)." >&2
    exit 1
fi

if [[ ! -f "$PRIMER_REF" ]]; then
    echo "ERROR: Primer reference file '$PRIMER_REF' not found." >&2
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input file '$INPUT' not found." >&2
    exit 1
fi

awk -v p1="$P1_COL" -v p2="$P2_COL" \
    -v pid="$PRIMER_ID_COL" -v pseq="$PRIMER_SEQ_COL" \
    'BEGIN { FS=OFS="\t"; in_ref=1 }
     FNR==1 && NR>1 { in_ref=0 }
     in_ref {
         if (FNR==1 && $pseq !~ /^[ACGTNacgtn]+$/) next
         map[$pid] = toupper($pseq)
         next
     }
     FNR==1 && !($p1 in map) {
         sub(/_id/, "_seq", $p1)
         sub(/_id/, "_seq", $p2)
         print
         next
     }
     {
         if (!($p1 in map)) print "Warning: ID " $p1 " not found in primer reference file" > "/dev/stderr"
         else $p1 = map[$p1]
         if (!($p2 in map)) print "Warning: ID " $p2 " not found in primer reference file" > "/dev/stderr"
         else $p2 = map[$p2]
         print
     }' "$PRIMER_REF" "$INPUT" > "$OUTPUT"

echo "Done. Output written to $OUTPUT"
