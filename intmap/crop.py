import subprocess
import pandas as pd
import regex
from rapidfuzz.distance import Levenshtein
from joblib import Parallel, delayed
import os
import glob
import math
from intmap.utils import *

def check_crop_input(ltr, linker, ltr_primer, ltr_error_rate, 
                     linker_error_rate, contamination):
    
    if ltr_primer is not None:
        ltr_prime_check = ltr_primer.replace(' ', '')
        ltr_prime_char = regex.sub('[ATGC]', '', ltr_prime_check.upper())
        if ltr_prime_char != '':
            raise Exception('Declared LTR primer may only contain A,T,G, and C nucleotides.')

        if len(ltr_prime_check) < 10:
            raise Exception('Declared LTR primer must be at least 10 nucleotides long.')
        
    ltr_check = ltr.replace(' ', '')
    ltr_char = regex.sub('[ATGC]', '', ltr_check.upper())
    if ltr_char != '':
        raise Exception('Declared LTR sequence may only contain A,T,G, and C nucleotides.')

    if len(ltr_check) < 10:
        raise Exception('Declared LTR sequence must be at least 10 nucleotides long.')
    
    if ltr_error_rate > 0.4:
        raise Exception('LTR sequence error rate cannot be > 0.4.')
    
    if linker:
        link_check = linker.replace(' ', '')
        link_char = regex.sub('[ATGC]', '', link_check.upper())
        if link_char != '':
            raise Exception('Declared linker sequence may only contain A,T,G and C nucleotides.')

        if len(link_check) < 10:
            raise Exception('Declared linker sequence must be at least 10 nucleotides long.')
        
        if linker_error_rate > 0.4:
            raise Exception('Linker sequence error rate cannot be > 0.4.')
            
    if contamination:
        for cc in contamination:
            cont = cc.replace(' ', '')
            cont_char = regex.sub('[ATGC]', '', cont.upper())
            if cont_char != '':
                raise Exception('Declared contamination(s) may only contain A,T,G, and C nucleotides.')
            
            if len(cont) < 10:
                raise Exception('Declared contamination(s) must be at least 10 nucleotides long.')
      
def parse_matches(ltr_matches, linker_matches,
                  ltr_umi_len, ltr_umi_offset, ltr_umi_pattern,
                  linker_umi_len, linker_umi_offset, linker_umi_pattern,
                  single_end, long_read):
  
    reads_to_remove = set()
  
    if ltr_umi_pattern:
        ltr_regex_pattern = regex.sub(r'N+', lambda m: f"[ATGC]{{{len(m.group())}}}", ltr_umi_pattern)
    if linker_umi_pattern:
        linker_regex_pattern = regex.sub(r'N+', lambda m: f"[ATGC]{{{len(m.group())}}}", linker_umi_pattern)
        
    if ltr_umi_len > 0:
        ltr_matches[4] = ltr_matches[4].fillna('N').str.pad(
            width=(ltr_umi_len + ltr_umi_offset), side='left', fillchar='N'
            ).str[-(ltr_umi_len + ltr_umi_offset):]
        if ltr_umi_offset > 0:
            ltr_matches[4] = ltr_matches[4].str[:-ltr_umi_offset]
        if ltr_umi_pattern:
            ltr_matches[4] = ltr_matches[4].where(ltr_matches[4].str.contains(ltr_regex_pattern, regex=True), 'N')
        reads_to_remove.update(ltr_matches[ltr_matches[4].str.contains('N', na=False)][0])
    else:
        ltr_matches[4] = ltr_matches[4].apply(lambda x: 'N')
    
    if not single_end or long_read:
        if linker_umi_len > 0:
            linker_matches[4] = linker_matches[4].fillna('N').str.pad(
                width=(linker_umi_len + linker_umi_offset), side='left', fillchar='N'
                ).str[-(linker_umi_len + linker_umi_offset):]
            if linker_umi_offset > 0:
                linker_matches[4] = linker_matches[4].str[:-linker_umi_offset]
            if linker_umi_pattern:
                linker_matches[4] = linker_matches[4].where(linker_matches[4].str.contains(linker_regex_pattern, regex=True), 'N')
            reads_to_remove.update(linker_matches[linker_matches[4].str.contains('N', na=False)][0])
        else:
            linker_matches[4] = linker_matches[4].apply(lambda x: 'N')

    if len(reads_to_remove) > 0:
        ltr_matches = ltr_matches[~ltr_matches[0].isin(reads_to_remove)]
        linker_matches = linker_matches[~linker_matches[0].isin(reads_to_remove)]
    
    ltr_matches = ltr_matches.sort_values(0).reset_index(drop=True)
    linker_matches = linker_matches.sort_values(0).reset_index(drop=True)

    if not single_end or long_read:
        if not ltr_matches[0].equals(linker_matches[0]):
            raise Exception('Mismatch in info-file read names.')
          
    return ltr_matches, linker_matches
      
def make_new_headers(ltr_data, linker_data, out_nm, single_end):
    combined_data = pd.DataFrame({
        'read_name': ltr_data.iloc[:, 0],
        'ltr_umi': ltr_data.iloc[:, 1],
        'linker_umi': linker_data.iloc[:, 1] if not linker_data.empty else 'N',
        'ltr_found': ltr_data.iloc[:, 2],
        'linker_found': linker_data.iloc[:, 2] if not linker_data.empty else 'N'
    })

    combined_data['header1'] = (
        '@' +
        combined_data['read_name'] +
        f'\tCO:Z:1:{out_nm}' +
        '\tRX:Z:' +
        combined_data['ltr_umi'].astype(str) +
        '-' +
        combined_data['linker_umi'].astype(str) +
        '\tOX:Z:' +
        combined_data['ltr_found'].astype(str) +
        '-' +
        combined_data['linker_found'].astype(str)
    )

    if not single_end:
        combined_data['header2'] = (
            '@' +
            combined_data['read_name'] +
            f'\tCO:Z:2:{out_nm}' +
            '\tRX:Z:' +
            combined_data['ltr_umi'].astype(str) +
            '-' +
            combined_data['linker_umi'].astype(str) +
            '\tOX:Z:' +
            combined_data['ltr_found'].astype(str) +
            '-' +
            combined_data['linker_found'].astype(str)
        )

        headers = {
            row.read_name: {'header1': row.header1, 'header2': row.header2}
            for row in combined_data.to_records(index=False)
        }
    else:
        headers = {
            row.read_name: {'header1': row.header1}
            for row in combined_data.to_records(index=False)
        }

    return headers
                
def write_new_fastq(in_fq, out_fq, headers, type, nthr = 1, kept_occ = None):
    with open_file(in_fq, 'rt', nthr = nthr) as f_in, open_file(out_fq, 'wt', nthr = nthr) as f_out:
        buffer = []
        written = set()
        occ_count = {}
        i = 0
        for line in f_in:
            if line.startswith('@') and i % 4 == 0:
                read_name = line.strip().split()[0][1:]
                occ = occ_count.get(read_name, 0)
                occ_count[read_name] = occ + 1
                target_occ = kept_occ.get(read_name, 0) if kept_occ else 0
                if read_name in headers and read_name not in written and occ == target_occ:
                    buffer.append(f"{headers[read_name][type]}\n")
                    written.add(read_name)
                else:
                    i += 4
                    next(f_in, None)
                    next(f_in, None)
                    next(f_in, None)
                    continue
            else:
                buffer.append(line)

            i += 1

            if len(buffer) >= 1200:
                f_out.writelines(buffer)
                buffer = []

        if buffer:
            f_out.writelines(buffer)
            
def ttr_subseqs(sequence, min_ttr_len, error_rate, cap_ttr_error):
    seq_len = len(sequence)
    if cap_ttr_error:
        reduction = max(1, math.floor(min_ttr_len * error_rate))
    chunks = []
    while seq_len >= min_ttr_len:
        chunks.append(sequence[0:seq_len])
        if not cap_ttr_error:
            reduction = max(1, math.floor(seq_len * error_rate))
        seq_len -= reduction
    return chunks

def save_ttr_subseqs(seq1, seq2, min_ttr_len, error_rate, processed_directory, write, rc, leftmost, cap_ttr_error):
    all_sequences = []
    chunks = ttr_subseqs(seq1, min_ttr_len, error_rate, cap_ttr_error)
    all_sequences.extend(chunks)

    if seq2:
        chunks = ttr_subseqs(seq2, min_ttr_len, error_rate, cap_ttr_error)
        all_sequences.extend(chunks)
        
    if rc:
        all_sequences = [revcomp(seq) for seq in all_sequences]
    
    all_sequences = sorted(set(all_sequences), key = len, reverse = True)
        
    if write:
        outfile = os.path.join(processed_directory, 'ttr.fa')
        
        with open(f'{outfile}', 'w') as fasta_file:
            for i, seq in enumerate(all_sequences):
                fasta_file.write(f">Sequence_{i + 1}\n")
                fasta_file.write(f"{seq}\n")
    else:
        if not rc:
            adapter_strs = [f"\"{seq};min_overlap={min_ttr_len}{';rightmost' if not leftmost else ''}\"" for seq in all_sequences]
            adpt_pos = '-g'
        else:
            adapter_strs = [f'\"{seq};min_overlap=5\"' for seq in all_sequences]
            adpt_pos = '-a'
        return ' '.join(f"{adpt_pos} {adpt}" for adpt in adapter_strs)

def crop(r1, r2, ltr, ltr2, linker, ltr_error_rate, linker_error_rate, min_frag_len,
         ltr_umi_len, ltr_umi_offset, ltr_umi_pattern, linker_umi_len, linker_umi_offset, 
         linker_umi_pattern, processed_directory, out_nm, nthr, single_end, long_read, ttr, 
         min_ttr_len, cut_path, contam, leftmost, mixed, seqtk_path, cap_ttr_error):
        
    if r1 == r2:
        raise Exception('r1 and r2 file names are identical.')

    if not ttr:
        ltr_seq = f"\"-g {ltr};min_overlap={len(ltr)}{';rightmost' if not leftmost else ''}\""
    else:
        if ltr and ltr2:
            ttr_left = list(set([ltr[:min_ttr_len], ltr2[:min_ttr_len]]))
            if len(ttr_left) == 1:
                ltr_init = (f'\"-g {ttr_left[0]};min_overlap={min_ttr_len}\" ')
            else:
                ltr_init = (f'\"-g {ttr_left[0]};min_overlap={min_ttr_len}\" '
                            f'-g \"{ttr_left[1]};min_overlap={min_ttr_len}\"')
        else:
            ltr_init = f'"-g {ltr[:min_ttr_len]};min_overlap={min_ttr_len}"'

        ltr_trim0 = (
            f'{cut_path} '
            f'-e {ltr_error_rate} '
            f'-j {nthr} '
            f'--discard-untrimmed '
            f'-m {min_frag_len} '
            f'{ltr_init} '
            f'--action=none '
            f"--quiet -o {os.path.join(processed_directory, 'ttr_init_crop_Rltr.fq.gz')} "
            + (f"-p {os.path.join(processed_directory, 'ttr_init_crop_Rlink.fq.gz')} {r1} {r2}" if r2 else f"{r1}")
        )
    
        subprocess.call(ltr_trim0, shell = True)
        
        r1 = os.path.join(processed_directory, 'ttr_init_crop_Rltr.fq.gz')
        if r2:
            r2 = os.path.join(processed_directory, 'ttr_init_crop_Rlink.fq.gz')

        ltr_seq = save_ttr_subseqs(seq1 = ltr, seq2 = ltr2, min_ttr_len = min_ttr_len, 
                                   error_rate = ltr_error_rate, processed_directory = processed_directory,
                                   write = False, rc = False, leftmost = leftmost, cap_ttr_error = cap_ttr_error)
        
    if mixed and long_read:
        ltr_seq_rc = f"\"-a {revcomp(ltr)};min_overlap={len(ltr)}{';rightmost' if leftmost else ''}\""

        mixed_trim = (
            f'{cut_path} '
            f'-e {ltr_error_rate} '
            f'-j {nthr} '
            f'--discard-untrimmed '
            f'--action=none '
            f'-m {min_frag_len} '
            f'{ltr_seq} '
            f"--quiet -o {os.path.join(processed_directory, 'mixed_for.fq.gz')} "
            f"{r1}"
        )
        
        mixed_rc_trim = (
            f'{cut_path} '
            f'-e {ltr_error_rate} '
            f'-j {nthr} '
            f'--discard-untrimmed '
            f'--action=none '
            f'-m {min_frag_len} '
            f'{ltr_seq_rc} '
            f"--quiet -o {os.path.join(processed_directory, 'mixed_rev.fq.gz')} "
            f"{r1}"
        )
        
        mixed_ltr_rc = (
            f"{seqtk_path} seq -r {os.path.join(processed_directory, 'mixed_rev.fq.gz')} | gzip > " 
            f"{os.path.join(processed_directory, 'mixed_rc.fq.gz')}"
        )

        mixed_ltr_cat = (
            f"gzip -dc {os.path.join(processed_directory, 'mixed_for.fq.gz')} " 
            f"{os.path.join(processed_directory, 'mixed_rc.fq.gz')} | gzip > "
            f"{os.path.join(processed_directory, 'mixed_lr.fq.gz')}"
        )

        trim_call = f'{mixed_trim} && {mixed_rc_trim} && {mixed_ltr_rc} && {mixed_ltr_cat}'
        
        r1 = os.path.join(processed_directory, 'mixed_lr.fq.gz')
    else:
        trim_call = None

    ltr_trim1 = (
        f'{cut_path} '
        f'-e {ltr_error_rate} '
        f'-j {nthr} '
        f"--info-file {os.path.join(processed_directory, 'ltr_info.txt')} "
        f'--discard-untrimmed '
        f'-m {min_frag_len} '
        f'{ltr_seq} '
        f"--quiet -o {os.path.join(processed_directory, 'ltr1_crop_Rltr.fq.gz')} "
        + (f"-p {os.path.join(processed_directory, 'ltr1_crop_Rlink.fq.gz')} {r1} {r2}" if r2 else f"{r1}")
    )
    
    if not trim_call:
        trim_call = f'{ltr_trim1}'
    else:
        trim_call = f'{trim_call} && {ltr_trim1}'
    
    if mixed and r2:
        ltr_trim_mixed = (
            f'{cut_path} '
            f'-e {ltr_error_rate} '
            f'-j {nthr} '
            f"--info-file {os.path.join(processed_directory, 'ltr_mixed_info.txt')} "
            f'--discard-untrimmed '
            f'-m {min_frag_len} '
            f'{ltr_seq} '
            f"--quiet -o {os.path.join(processed_directory, 'ltr1_mixed_crop_Rltr.fq.gz')} "
            f"-p {os.path.join(processed_directory, 'ltr1_mixed_crop_Rlink.fq.gz')} {r2} {r1}"
        )
        
        mixed_cat1 = (
            f'gzip -dc {os.path.join(processed_directory, "ltr1_crop_Rltr.fq.gz")} '
            f'{os.path.join(processed_directory, "ltr1_mixed_crop_Rltr.fq.gz")} '
            f'| gzip > {os.path.join(processed_directory, "ltr1_mixed_Rltr_cat.fq.gz")} '
            f'&& mv {os.path.join(processed_directory, "ltr1_mixed_Rltr_cat.fq.gz")} '
            f'{os.path.join(processed_directory, "ltr1_crop_Rltr.fq.gz")}'
        )
        
        mixed_cat2 = (
            f"gzip -dc {os.path.join(processed_directory, 'ltr1_crop_Rlink.fq.gz')} "
            f"{os.path.join(processed_directory, 'ltr1_mixed_crop_Rlink.fq.gz')} " 
            f"| gzip > {os.path.join(processed_directory, 'ltr1_mixed_Rlink_cat.fq.gz')} "
            f"&& mv {os.path.join(processed_directory, 'ltr1_mixed_Rlink_cat.fq.gz')} "
            f"{os.path.join(processed_directory, 'ltr1_crop_Rlink.fq.gz')}"
        )
        
        mixed_cat3 = (
            f"cat {os.path.join(processed_directory, 'ltr_info.txt')} "
            f"{os.path.join(processed_directory, 'ltr_mixed_info.txt')} " 
            f"> {os.path.join(processed_directory, 'ltr_cat_info.txt')} "
            f"&& mv {os.path.join(processed_directory, 'ltr_cat_info.txt')} "
            f"{os.path.join(processed_directory, 'ltr_info.txt')}"
        )
        
        trim_call = (
            f'{trim_call} && {ltr_trim_mixed} && {mixed_cat1} && {mixed_cat2} && {mixed_cat3}'
        )

    if linker:
        ltr1_r1 = f"{os.path.join(processed_directory, 'ltr1_crop_Rltr.fq.gz')}"
        if r2:
            ltr1_r2 = f"{os.path.join(processed_directory, 'ltr1_crop_Rlink.fq.gz')}"

        ltr_trim2 = (
            f'{cut_path} '
            f'-e {ltr_error_rate} '
            f'-j {nthr} '
            f"--info-file {os.path.join(processed_directory, 'linker_rc_info.txt')} "
            f'{"--discard-untrimmed" if long_read else ""} '
            f'-m {min_frag_len} '
            f'-a \"{revcomp(linker)};'
            f'min_overlap={5}\" '
            f"--quiet -o {os.path.join(processed_directory, 'ltr2_crop_Rltr.fq.gz')} "
            + (f"-p {os.path.join(processed_directory, 'ltr2_crop_Rlink.fq.gz')} {ltr1_r1} {ltr1_r2}" if r2 else f"{ltr1_r1}")
        )

        trim_call = f'{trim_call} && {ltr_trim2}'
        
    subprocess.call(trim_call, shell = True)
    
    if contam:
        if not isinstance(contam, list):
            contam = [contam]
        contam_r1 = os.path.join(processed_directory, f'ltr{2 if linker else 1}_crop_Rltr.fq.gz')
        contam_out_r1 = os.path.join(processed_directory, 'contam_crop_Rltr.fq.gz')
        if r2:
            contam_r2 = os.path.join(processed_directory, f'ltr{2 if linker else 1}_crop_Rlink.fq.gz')
            contam_out_r2 = os.path.join(processed_directory, 'contam_crop_Rlink.fq.gz')
                
        contam_strs = [f'\"{seq};min_overlap=10\"' for seq in contam]
        contam_strs = ' '.join(f"-g {cc}" for cc in contam_strs)
        
        if r2:
            contam_mv = f'&& mv {contam_out_r1} {contam_r1} && mv {contam_out_r2} {contam_r2}'
        else:
            contam_mv = f'&& mv {contam_out_r1} {contam_r1}'

        contam_trim = (
            f'{cut_path} '
            f'-e 0.3 '
            f'-j {nthr} '
            f"--info-file {os.path.join(processed_directory, 'contam_info.txt')} "
            f'--discard-trimmed '
            f'{contam_strs} '
            f'--quiet -o {contam_out_r1} '
            + (f"-p {contam_out_r2} {contam_r1} {contam_r2} " if r2 else f"{contam_r1} ")
            + contam_mv
        )
        
        subprocess.call(contam_trim, shell = True)
    
    ltr_file = os.path.join(processed_directory, 'ltr_info.txt')
    linker_rc_file = os.path.join(processed_directory, 'linker_rc_info.txt')
    
    try:
        ltr_data = pd.read_csv(
            ltr_file,
            sep = '\t',
            header = None,
            names = [0,1,2,3,4,5],
            dtype = {0: str, 1: int, 2: str, 3: str, 4: str, 5: str},
            usecols = [0,1,2,3,4,5],
            on_bad_lines = 'skip',
            low_memory = True
        )
    except pd.errors.ParserError:
        try:
            ltr_data = pd.read_csv(
                ltr_file,
                sep = '\t',
                header = None,
                names = [0,1,2,3,4,5],
                dtype = {0: str, 1: int, 2: str, 3: str, 4: str, 5: str},
                usecols = [0,1,2,3,4,5],
                on_bad_lines = 'skip',
                low_memory = False
            )
        except pd.errors.ParserError:
            raise Exception('No LTR matches found.')
    ltr_data = ltr_data[ltr_data.iloc[:, 1] != -1]
    ltr_data[0] = ltr_data[0].str.split(' ').str[0]
    # Added to eliminate a rare but annoying "mixed-type" warning in columns 2, 3
    # Also hard-coded for 2, 3 as str
    ltr_data[2] = pd.to_numeric(ltr_data[2], errors = 'coerce')
    ltr_data[3] = pd.to_numeric(ltr_data[3], errors = 'coerce')
    
    kept_occ = {}
    if mixed:
        dup_mask = ltr_data[0].duplicated(keep=False)
        if dup_mask.any():
            ltr_data['occurence'] = ltr_data.groupby(0).cumcount()
            dup_data = ltr_data[dup_mask].copy()
            dup_data['edit_dist'] = dup_data[5].apply(
                lambda s: Levenshtein.normalized_distance(s, ltr)
            )

            def choose_dup(grp):
                min_dist = grp['edit_dist'].min()
                best = grp[grp['edit_dist'] == min_dist]
                if len(best) == 1:
                    return best
                elif not leftmost:
                    return best.loc[[best[2].idxmax()]]
                else:
                    return best.loc[[best[2].idxmin()]]

            chosen_dups = (
                dup_data.groupby(0, sort = False, group_keys = False)
                .apply(choose_dup)
                .drop(columns = 'edit_dist')
            )
            del dup_data
            kept_occ = chosen_dups.set_index(0)['occurence'].to_dict()
            ltr_data = pd.concat(
                [ltr_data[~dup_mask], chosen_dups], ignore_index = True
            )
            ltr_data = ltr_data.drop(columns=['occurence'])

    try:
        linker_rc_data = pd.read_csv(
            linker_rc_file,
            sep = '\t',
            header = None,
            names = [0,1,2,3,4,5,6],
            dtype = {0: str, 1: int, 2: str, 3: str, 4: str, 5: str, 6: str},
            usecols = [0,1,2,3,4,5,6],
            on_bad_lines = 'skip'
        )
    except pd.errors.ParserError:
        linker_rc_data = pd.read_csv(
            linker_rc_file,
            sep = '\t',
            header = None,
            names = [0,1,2,3],
            dtype = {0: str, 1: int},
            on_bad_lines = 'skip'
        )
        linker_rc_data[[4, 5, 6]] = None
        linker_rc_data[0] = linker_rc_data[0].str.split(' ').str[0]
    except (pd.errors.EmptyDataError, FileNotFoundError):
        linker_rc_data = pd.DataFrame(columns=[0,1,2,3,4,5,6])
    linker_rc_data[2] = pd.to_numeric(linker_rc_data[2], errors = 'coerce')
    linker_rc_data[3] = pd.to_numeric(linker_rc_data[3], errors = 'coerce')
    
    short_reads = set(linker_rc_data.loc[
        (linker_rc_data[1] != -1) & (linker_rc_data[4].str.len() < min_frag_len), 0])
    
    if contam:
        contam_data = pd.read_csv(
            os.path.join(processed_directory, 'contam_info.txt'), 
            sep = '\t',
            header = None,
            names = [0,1,2,3,4,5],
            dtype = {0: str, 1: int, 4: str, 5: str}, 
            usecols = [0,1,2,3,4,5],
            on_bad_lines = 'skip'
          )

        contam_data = contam_data.loc[contam_data[1] != -1]
        short_reads.update(contam_data[0])
        del contam_data

    ltr_data = ltr_data[~ltr_data[0].isin(short_reads)]
    ltr_data = ltr_data.drop(columns=[1,2,3])
    if single_end:
        linker_data = linker_rc_data.drop(columns=[1,2,3]).fillna('N')
        if long_read:
            linker_data[4] = linker_data[6].apply(
                lambda x: revcomp(x) if x != 'N' else 'N'
            )
        linker_data = linker_data.drop(columns=[6])
        linker_data = linker_data[~linker_rc_data[0].isin(short_reads)]
        del linker_rc_data
    else:
        del linker_rc_data
        if linker:
            linker_trim1 = (
                f'{cut_path} '
                f'-e {linker_error_rate} '
                f'-j {nthr} '
                f"--info-file {os.path.join(processed_directory, 'linker_info.txt')} "
                f'--discard-untrimmed '
                f'-m {min_frag_len} '
                f"-g \"{linker};min_overlap={len(linker)}{';rightmost' if not leftmost else ''}\" "
                f"--quiet  -o {os.path.join(processed_directory, 'linker1_crop_Rlink.fq.gz')} "
                f"-p {os.path.join(processed_directory, 'linker1_crop_Rltr.fq.gz')} "
                f"{os.path.join(processed_directory, 'ltr2_crop_Rlink.fq.gz')} {os.path.join(processed_directory, 'ltr2_crop_Rltr.fq.gz')}"
            )
            
            linker1_r1 = f"{os.path.join(processed_directory, 'linker1_crop_Rltr.fq.gz')}"
            linker1_r2 = f"{os.path.join(processed_directory, 'linker1_crop_Rlink.fq.gz')}"
            
            if not ttr:
                ltr_rc_seq = f'-a \"{revcomp(ltr)};min_overlap=5\"'
            else:
                ltr_rc_seq = save_ttr_subseqs(seq1 = ltr, seq2 = ltr2, min_ttr_len = min_ttr_len, 
                                              error_rate = 0, processed_directory = processed_directory,
                                              write = False, rc = True, leftmost = leftmost, cap_ttr_error = cap_ttr_error)

            linker_trim2 = (
                f'{cut_path} '
                f'-e {linker_error_rate} '
                f'-j {nthr} '
                f"--info-file {os.path.join(processed_directory, 'ltr_rc_info.txt')} "
                f'-m {min_frag_len} '
                f'{ltr_rc_seq} '
                f"--quiet -o {os.path.join(processed_directory, 'linker2_crop_Rlink.fq.gz')} "
                f"-p {os.path.join(processed_directory, 'linker2_crop_Rltr.fq.gz')} {linker1_r2} {linker1_r1}"
            )

            linker_trim = f'{linker_trim1} && {linker_trim2}'
            subprocess.call(linker_trim, shell = True)
            
            linker_file = os.path.join(processed_directory, 'linker_info.txt')
            ltr_rc_file = os.path.join(processed_directory, 'ltr_rc_info.txt')
            
            try:
                linker_data = pd.read_csv(
                    linker_file,
                    sep = '\t',
                    header = None,
                    names = [0,1,2,3,4,5],
                    dtype = {0: str, 1: int, 2: str, 3: str, 4: str, 5: str},
                    usecols = [0,1,2,3,4,5],
                    on_bad_lines = 'skip',
                    low_memory = True
                )
            except pd.errors.ParserError:
                try:
                    linker_data = pd.read_csv(
                        linker_file,
                        sep = '\t',
                        header = None,
                        names = [0,1,2,3,4,5],
                        dtype = {0: str, 1: int, 2: str, 3: str, 4: str, 5: str},
                        usecols = [0,1,2,3,4,5],
                        on_bad_lines = 'skip',
                        low_memory = False
                    )
                except pd.errors.ParserError:
                    raise Exception('No linker matches found.')
            linker_data[0] = linker_data[0].str.split(' ').str[0]
            linker_data[2] = pd.to_numeric(linker_data[2], errors = 'coerce')
            linker_data[3] = pd.to_numeric(linker_data[3], errors = 'coerce')

            if kept_occ:
                linker_data['occurence'] = linker_data.groupby(0).cumcount()
                dup_reads = linker_data[0].isin(kept_occ)
                if dup_reads.any():
                    dup_linker = linker_data[dup_reads].copy()
                    dup_linker = dup_linker[
                        dup_linker.apply(lambda r: r['occurence'] == kept_occ[r[0]], axis=1)
                    ]
                    linker_data = pd.concat(
                        [linker_data[~dup_reads], dup_linker], ignore_index=True
                    )
                linker_data = linker_data.drop(columns=['occurence'])

            try:
                ltr_rc_data = pd.read_csv(
                    ltr_rc_file, 
                    sep = '\t',
                    header = None, 
                    names = [0,1,2,3,4,5],
                    dtype = {0: str, 1: int, 4: str, 5: str}, 
                    usecols = [0,1,2,3,4,5],
                    on_bad_lines = 'skip'
                )
            except pd.errors.ParserError:
                ltr_rc_data = pd.read_csv(
                    ltr_rc_file, 
                    sep = '\t',
                    header = None,
                    names = [0,1,2,3],
                    dtype = {0: str, 1: int},
                    on_bad_lines = 'skip'
                )
                ltr_rc_data[[4, 5]] = None
            
            short_reads.update(ltr_rc_data.loc[
                (ltr_rc_data[1] != -1) & (ltr_rc_data[4].str.len() < min_frag_len), 0])
            del ltr_rc_data
            
            linker_data = linker_data[~linker_data[0].isin(short_reads)]
            linker_data = linker_data[linker_data.iloc[:, 1] != -1]
            linker_data = linker_data.drop(columns=[1,2,3])
            ltr_data = ltr_data[ltr_data.iloc[:, 0].isin(linker_data.iloc[:, 0])]
        else:
            no_linker_r2 = os.path.join(processed_directory, 'ltr1_crop_Rlink.fq.gz')
            no_linker_r1 = os.path.join(processed_directory, 'ltr1_crop_Rltr.fq.gz')

            if not ttr:
                ltr_rc_seq = f'-a \"{revcomp(ltr)};min_overlap=5\"'
            else:
                ltr_rc_seq = save_ttr_subseqs(seq1 = ltr, seq2 = ltr2, min_ttr_len = min_ttr_len,
                                              error_rate = 0, processed_directory = processed_directory,
                                              write = False, rc = True, leftmost = leftmost, cap_ttr_error = cap_ttr_error)

            ltr_rc_trim = (
                f'{cut_path} '
                f'-e {ltr_error_rate} '
                f'-j {nthr} '
                f"--info-file {os.path.join(processed_directory, 'ltr_rc_info.txt')} "
                f'-m {min_frag_len} '
                f'{ltr_rc_seq} '
                f"--quiet -o {os.path.join(processed_directory, 'linker2_crop_Rlink.fq.gz')} "
                f"-p {os.path.join(processed_directory, 'linker2_crop_Rltr.fq.gz')} {no_linker_r2} {no_linker_r1}"
            )

            subprocess.call(ltr_rc_trim, shell=True)

            ltr_rc_file = os.path.join(processed_directory, 'ltr_rc_info.txt')

            try:
                ltr_rc_data = pd.read_csv(
                    ltr_rc_file,
                    sep='\t',
                    header=None,
                    names=[0,1,2,3,4,5],
                    dtype={0: str, 1: int, 2: str, 3: str, 4: str, 5: str},
                    usecols=[0,1,2,3,4,5],
                    on_bad_lines='skip'
                )
            except pd.errors.ParserError:
                ltr_rc_data = pd.read_csv(
                    ltr_rc_file,
                    sep='\t',
                    header=None,
                    names=[0,1,2,3],
                    dtype={0: str, 1: int},
                    on_bad_lines='skip'
                )
                ltr_rc_data[[4, 5]] = None
            except (pd.errors.EmptyDataError, FileNotFoundError):
                ltr_rc_data = pd.DataFrame(columns=[0,1,2,3,4,5])

            short_reads.update(ltr_rc_data.loc[
                (ltr_rc_data[1] != -1) & (ltr_rc_data[4].str.len() < min_frag_len), 0])
            del ltr_rc_data

            linker_data = ltr_data[[0]].copy()
            linker_data[4] = 'N'
            linker_data[5] = 'N'

    ltr_data, linker_data = parse_matches(
        ltr_matches = ltr_data, 
        linker_matches = linker_data,
        ltr_umi_len = ltr_umi_len, 
        ltr_umi_offset = ltr_umi_offset, 
        ltr_umi_pattern = ltr_umi_pattern,
        linker_umi_len = linker_umi_len, 
        linker_umi_offset = linker_umi_offset, 
        linker_umi_pattern = linker_umi_pattern,
        single_end = single_end,
        long_read = long_read
    )

    headers = make_new_headers(
        ltr_data = ltr_data,
        linker_data = linker_data,
        out_nm = out_nm,
        single_end = single_end
    )
    
    del ltr_data
    del linker_data

    if not single_end:
        if linker:
            in1 = os.path.join(processed_directory, 'linker2_crop_Rltr.fq.gz')
            in2 = os.path.join(processed_directory, 'linker2_crop_Rlink.fq.gz')
        else:
            in1 = os.path.join(processed_directory, 'linker2_crop_Rltr.fq.gz')
            in2 = os.path.join(processed_directory, 'linker2_crop_Rlink.fq.gz')
        out1 = os.path.join(processed_directory, f'{out_nm}_Rltr_cropped.fq.gz')
        out2 = os.path.join(processed_directory, f'{out_nm}_Rlink_cropped.fq.gz')

        Parallel(n_jobs = 2 if nthr > 1 else 1)(
            delayed(write_new_fastq)(
                in_fq = file_in,
                out_fq = file_out,
                headers = headers,
                type = type,
                nthr = math.floor(nthr / 2),
                kept_occ = kept_occ
            ) for (file_in, file_out, type) in [(in1, out1, 'header1'), (in2, out2, 'header2')]
          )
    else:
        in1 = os.path.join(processed_directory, f'ltr{2 if linker else 1}_crop_Rltr.fq.gz')
        out1 = os.path.join(processed_directory, f'{out_nm}_Rltr_cropped.fq.gz')
        
        write_new_fastq(
            in_fq = in1,
            out_fq = out1,
            headers = headers,
            type = 'header1',
            nthr = nthr,
            kept_occ = kept_occ
        )
    
    for pattern in [
        os.path.join(processed_directory, '*_crop_R*.fq.gz'),
        os.path.join(processed_directory, '*info.txt'),
        os.path.join(processed_directory, 'ttr.fa'),
        os.path.join(processed_directory, 'contam.fa'),
        os.path.join(processed_directory, 'mixed_*.fq.gz'),
        os.path.join(processed_directory, '*_cat.fq.gz')
    ]:
        for tmp_file in glob.glob(pattern):
            os.remove(tmp_file)