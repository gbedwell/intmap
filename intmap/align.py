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
            raise ValueError('Given executable does not appear to be {}.'.format(ex_name))
    
    if os.path.isdir(path):
        path = os.path.join(path, ex_name)
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
        else:
            raise ValueError('The executable does not exist at the given location.')

# Align cropped reads
def align(bt2_path, sam_path, bt2_idx_dir, bt2_idx_name, 
        file1, file2, name, min_frag_len, max_frag_len, nthr,
        v_idx_dir, v_idx_name):
    
    bt2_idx = '{}/{}'.format(bt2_idx_dir, bt2_idx_name)
    
    v_idx = '{}/{}'.format(v_idx_dir, v_idx_name)
        
    bt2_out = name + '_full_aln.bam'
    ee_out = name + '_global_aln.bam'
    nonv_R1 = name + '_non_viral_R1.fq'
    nonv_R2 = name + '_non_viral_R2.fq'
    nonv_null = name + '_non_viral_null.fq'
    nonv_single = name + '_non_viral_single.fq'
    unmap_R1 = name + '_unmapped_R1.fq'
    unmap_R2 = name + '_unmapped_R2.fq'
    unmap_null = name + '_unmapped_null.fq'
    unmap_single = name + '_unmapped_single.fq'
    loc_out = name + '_local_aln.bam'
    final_out = name + '_aln.bam'
    index_out = name + '_aln.bai'
    
    v_aln_cmd = ('{} --very-sensitive -x {} -1 {} -2 {} -p {} -I {} -X {} '
            '--sam-append-comment --no-mixed --no-discordant | ').format(
            bt2_path, v_idx, file1, file2, nthr, min_frag_len, max_frag_len)
    rm_v_cmd1 = '{} view -@ {} -h -f 4 -b - | '.format(sam_path, nthr)
    rm_v_cmd2 = '{} collate -O -@ {} - | '.format(sam_path, nthr)
    rm_v_cmd3 = '{} fastq -@ {} -T RX,OX -1 {} -2 {} -0 {} -n -s {} -'.format(sam_path, nthr, nonv_R1, 
                                                                                nonv_R2, nonv_null, nonv_single)
    
    remove_viral_cmd = v_aln_cmd + rm_v_cmd1 + rm_v_cmd2 + rm_v_cmd3

    bt2_cmd = ('{} --very-sensitive -x {} -1 {} -2 {} -p {} -I {} -X {} '
            '--sam-append-comment --no-mixed --no-discordant | ').format(
            bt2_path, bt2_idx, nonv_R1, nonv_R2, nthr, min_frag_len, max_frag_len)
    
    bt2_to_bam = '{} view -@ {} -h -b -o {}'.format(sam_path, nthr, bt2_out)
    
    bt2_cmd = bt2_cmd + bt2_to_bam
    
    ee_cmd1 = '{} view -@ {} -h -f 2 -F 2308 -b {} | '.format(sam_path, nthr, bt2_out)
    ee_cmd2 = '{} sort -@ {} -o {} -'.format(sam_path, nthr, ee_out)
        
    global_cmd = ee_cmd1 + ee_cmd2
    
    unmap_cmd1 = '{} view -@ {} -h -f 4 -b {} | '.format(sam_path, nthr, bt2_out)
    unmap_cmd2 = '{} collate -O -@ {} - | '.format(sam_path, nthr)
    unmap_cmd3 = '{} fastq -@ {} -T RX,OX -1 {} -2 {} -0 {} -n -s {} -'.format(sam_path, nthr, unmap_R1, 
                                                                                unmap_R2, unmap_null, unmap_single)
    
    unmap_cmd = unmap_cmd1 + unmap_cmd2 + unmap_cmd3
    
    bt2_loc = ('{} --local --very-sensitive-local -x {} -1 {} -2 {} -p {} -I {} -X {} '
            '--sam-append-comment --no-mixed --no-discordant --no-unal | ').format(
            bt2_path, bt2_idx, unmap_R1, unmap_R2, nthr, min_frag_len, max_frag_len)
    
    loc_cmd1 = '{} view -@ {} -h -f 2 -F 2308 -b - | '.format(sam_path, nthr)
    loc_cmd2 = '{} sort -@ {} -o {} -'.format(sam_path, nthr, loc_out)
    
    local_cmd = bt2_loc + loc_cmd1 + loc_cmd2
        
    print('Removing viral reads from cropped output...')
    subprocess.call(remove_viral_cmd, shell=True)
    
    print('\n')
    print('Performing first-pass global (end-to-end) alignment...')
    subprocess.call(bt2_cmd, shell=True)
    
    print('\n')
    print('Cleaning up globally aligned reads...')
    subprocess.call(global_cmd, shell=True)
    
    print('\n')
    print('Extracting unmapped reads from global alignment...')
    subprocess.call(unmap_cmd, shell=True)
    
    print('\n')
    print('Performing second-pass local alignment on unaligned reads...')
    subprocess.call(local_cmd, shell=True)
    
    merge_cmd1 = '{} merge -n -@ {} - {} {} | '.format(sam_path, nthr, ee_out, loc_out)
    merge_cmd2 = '{} sort -@ {} -o {} -'.format(sam_path, nthr, final_out)
    
    merge_cmd = merge_cmd1 + merge_cmd2
    
    idx_cmd = '{} index -b -@ {} {} -o {}'.format(sam_path, nthr, final_out, index_out)

    print('\n')
    print('Merging bam files and indexing the output...')
    subprocess.call(merge_cmd, shell=True)
    subprocess.call(idx_cmd, shell=True)
    
    for f in [bt2_out, ee_out, unmap_R1, unmap_R2, unmap_null, unmap_single, loc_out,
            nonv_R1, nonv_R2, nonv_null, nonv_single]:
        os.remove(f)