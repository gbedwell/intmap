#!/usr/local/bin/python

import os
import subprocess

# Make sure that the given genome index exists and is a bowtie2 index.
# Adapted from UMAP run_bowtie.py
def check_genome_idx(idx_dir, idx_name):
    extensions = ['bt2', 'amb', 'ann', 'pac', 'sa', 'bwt']
    bt2_idx_files = [
        file for file in 
        os.listdir(idx_dir) if idx_name in file and any(file.endswith(ext) for ext in extensions)
        ]
    if len(bt2_idx_files) < 6: # changed to 5 for BWA compatibility, Bowtie2 has 6 index files
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
                file1, file2, out_nm, min_frag_len, max_frag_len, nthr,
                v_idx_dir, v_idx_name):
    
    bt2_idx = f'{bt2_idx_dir}/{bt2_idx_name}'
    
    v_idx = f'{v_idx_dir}/{v_idx_name}'
        
    nonv_R1 = out_nm + '_non_viral_R1.fq'
    nonv_R2 = out_nm + '_non_viral_R2.fq'
    nonv_null = out_nm + '_non_viral_null.fq'
    nonv_single = out_nm + '_non_viral_single.fq'
    final_out = out_nm + '_aln.bam'
    index_out = out_nm + '_aln.bai'
    
    seed_len = 22
    
    v_aln_cmd = (f'{bt2_path} --end-to-end -x {v_idx} -1 {file1} -2 {file2} -p {nthr} -I {min_frag_len} -X {max_frag_len} '
                    f'-L {seed_len} -N 0 -D 20 -R 3 -i S,1,0.50 --sam-append-comment --no-mixed --no-discordant | ')
    # v_aln_cmd = (f'{bwa_path} mem -t {nthr} -k {seed_len} -C  -L 1000,1000 {v_idx} {file1} {file2} | ')
    rm_v_cmd1 = f'{sam_path} view -@ {nthr} -h -f 4 -b - | '
    rm_v_cmd2 = f'{sam_path} collate -O -@ {nthr} - | '
    rm_v_cmd3 = f'{sam_path} fastq -@ {nthr} -T RX,OX -1 {nonv_R1} -2 {nonv_R2} -0 {nonv_null} -n -s {nonv_single} -'
    
    remove_viral_cmd = v_aln_cmd + rm_v_cmd1 + rm_v_cmd2 + rm_v_cmd3

    bt2_cmd1 = (f'{bt2_path} --end-to-end -x {bt2_idx} -1 {nonv_R1} -2 {nonv_R2} -p {nthr} -I {min_frag_len} -X {max_frag_len} '
                f'-L {seed_len} -N 0 -D 20 -R 3 -i S,1,0.50 --sam-append-comment --no-mixed --no-discordant | ')
    # bwa_cmd1 = (f'{bwa_path} mem -t {nthr} -k {seed_len} -C  -L 1000,1000 {bwa_idx} {file1} {file2} | ')
    bt2_cmd2 = f'{sam_path} view -@ {nthr} -h -f 2 -F 2316 -b | '
    bt2_cmd3 = f'{sam_path} sort -@ {nthr} -o {final_out} -'
    
    global_cmd = bt2_cmd1 + bt2_cmd2 + bt2_cmd3
    
    idx_cmd = f'{sam_path} index -b -@ {nthr} {final_out} -o {index_out}'
        
    print('Removing viral reads from cropped output...')
    subprocess.call(remove_viral_cmd, shell = True, stderr = subprocess.STDOUT)
    
    print('\n')
    print('Performing global (end-to-end) alignment...')
    subprocess.call(global_cmd, shell = True, stderr = subprocess.STDOUT)

    print('\n')
    print('Indexing output bam file...')
    subprocess.call(idx_cmd, shell = True, stderr = subprocess.STDOUT)
    
    for f in [nonv_R1, nonv_R2, nonv_null, nonv_single]:
        os.remove(f)