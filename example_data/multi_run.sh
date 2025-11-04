if [ ! -f "multi_args.json" ]; then
    # If multi_args.json does not exist, run the first command
    intmap_multi \
      -s multi_setup.txt \
      -n multi_args > multi_output.txt
else
    # If multi_args.json already exists, run the second command
    intmap_multi \
      -a multi_args.json > multi_output.txt
fi