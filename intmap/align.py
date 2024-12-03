#!/usr/local/bin/python

import os
import subprocess

# Make sure that the given genome index exists and is a bowtie2 index.
# Adapted from UMAP run_bowtie.py
def check_genome_idx(idx_dir, idx_name):
    bt2_idx_files = [
        file for file in
        os.listdir(idx_dir) if idx_name in file and file.endswith('.bt2')]
    if len(bt2_idx_files) < 6:
        raise ValueError("Given bowtie2 index does not exist.")
    
# Make sure that the given bowtie2 and samtools exectuables exist.
def check_executable(path, ex_name):
    if os.path.isfile(path) and os.access(path, os.X_OK):
        if os.path.basename(path) == ex_name:
            return path
        else:
            raise ValueError(f'Given executable does not appear to be {ex_name}.')
    
    if os.path.isdir(path):
        path = os.path.join(path, ex_name)
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
        else:
            raise ValueError('The executable does not exist at the given location.')

# Align cropped reads
def align_global(bt2_path, sam_path, bt2_idx_dir, bt2_idx_name, 
                file1, file2, name, min_frag_len, max_frag_len, nthr,
                v_idx_dir, v_idx_name):
    
    bt2_idx = f'{bt2_idx_dir}/{bt2_idx_name}'
    
    v_idx = f'{v_idx_dir}/{v_idx_name}'
        
    nonv_R1 = name + '_non_viral_R1.fq'
    nonv_R2 = name + '_non_viral_R2.fq'
    nonv_null = name + '_non_viral_null.fq'
    nonv_single = name + '_non_viral_single.fq'
    final_out = name + '_aln.bam'
    index_out = name + '_aln.bai'
    
    v_aln_cmd = (f'{bt2_path} --very-sensitive -x {v_idx} -1 {file1} -2 {file2} -p {nthr} -I {min_frag_len} -X {max_frag_len} '
                    '--sam-append-comment --no-mixed --no-discordant | ')
    rm_v_cmd1 = f'{sam_path} view -@ {nthr} -h -f 4 -b - | '
    rm_v_cmd2 = f'{sam_path} collate -O -@ {nthr} - | '
    rm_v_cmd3 = f'{sam_path} fastq -@ {nthr} -T RX,OX -1 {nonv_R1} -2 {nonv_R2} -0 {nonv_null} -n -s {nonv_single} -'
    
    remove_viral_cmd = v_aln_cmd + rm_v_cmd1 + rm_v_cmd2 + rm_v_cmd3

    bt2_cmd1 = (f'{bt2_path} --very-sensitive -x {bt2_idx} -1 {nonv_R1} -2 {nonv_R2} -p {nthr} -I {min_frag_len} -X {max_frag_len} '
                '--sam-append-comment --no-unal --no-mixed --no-discordant | ')
    
    bt2_cmd2 = f'{sam_path} view -@ {nthr} -h -f 2 -F 2308 -b | '
    bt2_cmd3 = f'{sam_path} sort -@ {nthr} -o {final_out} -'
    
    global_cmd = bt2_cmd1 + bt2_cmd2 + bt2_cmd3
    
    idx_cmd = f'{sam_path} index -b -@ {nthr} {final_out} -o {index_out}'
        
    print('Removing viral reads from cropped output...')
    subprocess.call(remove_viral_cmd, shell=True)
    
    print('\n')
    print('Performing global (end-to-end) alignment...')
    subprocess.call(global_cmd, shell=True)

    print('\n')
    print('Indexing output bam file...')
    subprocess.call(idx_cmd, shell=True)
    
    for f in [nonv_R1, nonv_R2, nonv_null, nonv_single]:
        os.remove(f)