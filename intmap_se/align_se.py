import os
import subprocess

# Adapted from UMAP run_bowtie.py
def check_genome_idx(idx_dir, idx_name, long_read):
    # For BWA: ['amb', 'ann', 'pac', 'sa', 'bwt']
    extensions = ['bt2', 'mmi']
    idx_files = [
        file for file in 
        os.listdir(idx_dir) if idx_name in file and any(file.endswith(ext) for ext in extensions)
        ]
    
    if not long_read:
        if len(idx_files) < 6:
            raise ValueError("Given bowtie2 index does not exist.")
    else:
        mmi_files = [file for file in idx_files if file.endswith('mmi')]
        if not mmi_files:
            raise ValueError("Given minimap2 index does not exist.")
    
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

def align_global(bt2_path, sam_path, bt2_idx_dir, bt2_idx_name, 
                file1, out_nm, min_frag_len, max_frag_len, nthr,
                v_idx_dir, v_idx_name, aln_seed_len, aln_seed_mismatch):
    
    bt2_idx = f'{bt2_idx_dir}/{bt2_idx_name}'
    
    v_idx = f'{v_idx_dir}/{v_idx_name}'
        
    nonv_R1 = out_nm + '_non_viral_R1.fq'
    final_out = out_nm + '_aln.bam'
    index_out = out_nm + '_aln.bai'

    v_aln_cmd = (f'{bt2_path} --end-to-end -x {v_idx} -U {file1} -p {nthr} ' 
                    f'-I {min_frag_len} -X {max_frag_len} -L {aln_seed_len} -N {aln_seed_mismatch} '
                    f'-D 20 -R 3 -i S,1,0.50 --sam-append-comment | ')
    # v_aln_cmd = (f'{bwa_path} mem -t {nthr} -k {seed_len} -C  -L 1000,1000 {v_idx} {file1} {file2} | ')
    rm_v_cmd1 = f'{sam_path} view -@ {nthr} -h -f 4 -b - | '
    rm_v_cmd2 = f'{sam_path} collate -O -@ {nthr} - | '
    rm_v_cmd3 = f'{sam_path} fastq -@ {nthr} -T RX,OX -0 {nonv_R1} -n -'
    
    remove_viral_cmd = v_aln_cmd + rm_v_cmd1 + rm_v_cmd2 + rm_v_cmd3

    bt2_cmd1 = (f'{bt2_path} --end-to-end -x {bt2_idx} -U {nonv_R1} -p {nthr} '
                f'-I {min_frag_len} -X {max_frag_len} -L {aln_seed_len} -N {aln_seed_mismatch} '
                f'-D 20 -R 3 -i S,1,0.50 --sam-append-comment | ')
    # bwa_cmd1 = (f'{bwa_path} mem -t {nthr} -k {seed_len} -C  -L 1000,1000 {bwa_idx} {file1} {file2} | ')
    bt2_cmd2 = f'{sam_path} view -@ {nthr} -h -b | '
    bt2_cmd3 = f'{sam_path} sort -@ {nthr} -o {final_out} -'
    
    global_cmd = bt2_cmd1 + bt2_cmd2 + bt2_cmd3
    
    idx_cmd = f'{sam_path} index -b -@ {nthr} {final_out} -o {index_out}'
        
    print('Removing viral reads from cropped output...', flush = True)
    subprocess.call(remove_viral_cmd, shell = True, stderr = subprocess.STDOUT)
    
    print('\n', flush = True)
    print('Performing global (end-to-end) alignment...', flush = True)
    subprocess.call(global_cmd, shell = True, stderr = subprocess.STDOUT)

    print('\n', flush = True)
    print('Indexing output bam file...', flush = True)
    subprocess.call(idx_cmd, shell = True, stderr = subprocess.STDOUT)
    
    for f in [nonv_R1]:
        os.remove(f)
        
        
def align_lr(mm_path, sam_path, bt2_idx_dir, bt2_idx_name, 
            file1, out_nm, nthr, v_idx_dir, v_idx_name, lr_type):
    
    bt2_idx = f'{bt2_idx_dir}/{bt2_idx_name}'
    
    v_idx = f'{v_idx_dir}/{v_idx_name}'
        
    nonv_R1 = out_nm + '_non_viral_R1.fq'
    final_out = out_nm + '_aln.bam'
    index_out = out_nm + '_aln.bai'

    v_aln_cmd = (f'{mm_path} -Layx {lr_type} -t {nthr} --MD --secondary no {v_idx}.mmi {file1} | ')
    rm_v_cmd1 = f'{sam_path} view -@ {nthr} -h -f 4 -b - | '
    rm_v_cmd2 = f'{sam_path} collate -O -@ {nthr} - | '
    rm_v_cmd3 = f'{sam_path} fastq -@ {nthr} -T RX,OX -0 {nonv_R1} -n -'
    
    remove_viral_cmd = v_aln_cmd + rm_v_cmd1 + rm_v_cmd2 + rm_v_cmd3

    bt2_cmd1 = (f'{mm_path} -Layx {lr_type} -t {nthr} --MD --sam-hit-only {bt2_idx}.mmi {nonv_R1} | ')
    bt2_cmd2 = f'{sam_path} view -@ {nthr} -h -F 4 -b | '
    bt2_cmd3 = f'{sam_path} sort -@ {nthr} -o {final_out} -'
    
    global_cmd = bt2_cmd1 + bt2_cmd2 + bt2_cmd3
    
    idx_cmd = f'{sam_path} index -b -@ {nthr} {final_out} -o {index_out}'
        
    print('Removing viral reads from cropped output...', flush = True)
    subprocess.call(remove_viral_cmd, shell = True, stderr = subprocess.STDOUT)
    
    print('\n', flush = True)
    print('Performing global (end-to-end) alignment...', flush = True)
    subprocess.call(global_cmd, shell = True, stderr = subprocess.STDOUT)

    print('\n', flush = True)
    print('Indexing output bam file...', flush = True)
    subprocess.call(idx_cmd, shell = True, stderr = subprocess.STDOUT)
    
    for f in [nonv_R1]:
        os.remove(f)