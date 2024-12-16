import subprocess
import argparse
import os
import shutil

parser = argparse.ArgumentParser(description='Download usable fastq files from SRA.')
parser.add_argument('-accessions_file', '-a',
                    type = str, 
                    help = '''Path to the file holding the accession numbers of interest. 
                    Expects 1 accession number per line.''')
parser.add_argument('-output_dir', '-o',
                    type = str, 
                    help = '''The output path. Used for both the final output and the temporary directories.''')
parser.add_argument('-nthr',
                    type = int, 
                    help = '''The number of cores to use.''')

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

with open(args.accessions_file, 'r') as accessions:
    sra = [line.strip() for line in accessions]

for acc in sra:
    prefetch = 'prefetch ' + acc + ' -O ' + args.output_dir
    subprocess.call(prefetch, shell=True)
    
    # Run prefetch twice on purpose.
    # sra-tools is spotty and this is the safest way to ensure
    # that all of the requisite information is downloaded
    subprocess.call(prefetch, shell=True)
    
    f_dump = ("fasterq-dump {}/{} -e {} -O {} -t {} "
            "--split-3 --skip-technical --seq-defline '@$sn:$sg'").format(
                args.output_dir, acc, args.nthr, args.output_dir, args.output_dir)
    
    subprocess.call(f_dump, shell=True)
    
    fastq_1 = os.path.join(args.output_dir, f'{acc}_1.fastq')
    fastq_2 = os.path.join(args.output_dir, f'{acc}_2.fastq')
    
    subprocess.call(f'gzip {fastq_1}', shell=True)
    subprocess.call(f'gzip {fastq_2}', shell=True)
    
    acc_dir = os.path.join(args.output_dir, acc)
    if os.path.exists(acc_dir):
        shutil.rmtree(acc_dir)
