"""
script to annotate vaf of panel of normal bams
2021.03.01 cjyoon
2021.03.18 cjyoon pon only calculate the specific alt base for calculation. 
"""

import os
import subprocess
import shlex
import re
import argparse
import pandas as pd
import numpy as np
import multiprocessing as mp
import pysam


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    parser.add_argument('-p', '--pon_list_tsv', required=True, help='text file that lists of PON bam files to be used (parent bamf files)')
    parser.add_argument('-r', '--reference', required=True, help='Reference fasta file')
    parser.add_argument('-n', '--multicore', required=False, type=int, default=1, help='Number of multiprocessor to parallelize computation')
    parser.add_argument('-c', '--column_name', required=False, default='pon_vaf', help='New columan name to be added for PON vafs')
    parser.add_argument('-e', '--exclusion_list', required=False, nargs='+', default=[], help='Bam files to exclude from PON calculation')
    parser.add_argument('--polyn', required=False, default=0, type=int, help='add polyN region info. run only when the value is 1')

    args = vars(parser.parse_args())

    return args['input'], args['output_dir'], args['pon_list_tsv'], args['reference'], args['multicore'], args['column_name'], args['exclusion_list'], args['polyn']


def annotate_polyn(query_position, reference):
    chrom, pos = query_position.split(':')
    pos = int(pos)
    ref = pysam.FastaFile(reference)
    up, down = ref.fetch(chrom,pos-2,pos-1), ref.fetch(chrom,pos,pos+1)
    n_up, n_down = 1, 1
    while True:
        if ref.fetch(chrom,pos-n_up-2,pos-n_up-1) == up: n_up+=1
        else: break
    while True:
        if ref.fetch(chrom,pos+n_down,pos+n_down+1) == down: n_down+=1
        else: break
    return f'{up}:{n_up}:{down}:{n_down}'

def vaf(alt_count, depth):
    """handles division by 0 error"""
    if depth!=0:
        return float(alt_count/depth)
    else:
        return 0


def calculate_pon_vaf(query_position, ref_alt, normal_bam_list, reference):
    '''calculates the VAF in the normal Bams in the list of a given query_position'''
    bamlistString = ' '.join(normal_bam_list)
    chrom, pos = query_position.split(':')
    ref, alt = ref_alt.split('_')
    
    alt_ = ''
    if len(ref) == len(alt):
        alt_  = alt
    elif len(ref) > len(alt):
        alt_ = 'del'
    elif len(ref) < len(alt):
        alt_ = 'ins'
        
    
        
    query_position = f'{chrom}:{pos}-{pos}'
    cmd = f'samtools mpileup -Q 0 -q 0 -f {reference} -r {query_position} {bamlistString}'
    mpileup = subprocess.check_output(shlex.split(cmd), stderr=subprocess.DEVNULL)
    split_mpileup = mpileup.decode("utf-8").split('\t'); #print(split_mpileup)
    totalDepth = 0 
    mismatches = 0
    totalCharacters = 0
    error_lists = [0 for i in range(0, len(normal_bam_list))]
    for i in range(1, int(len(split_mpileup)/3)):
        
        # initialize mismatch dict for each bam
        mismatch_dict = dict({'A': 0, 'C': 0, 'G': 0, 'T': 0, 'ins': 0, 'del': 0})
        
        bases_idx = 3*i + 1
        depths_idx = 3*i
        depth = int(split_mpileup[depths_idx])
        totalDepth += depth
        mpiledup = split_mpileup[bases_idx].upper()
        insertions = re.findall(r'\+[0-9]+[ACGTNacgtn]+', mpiledup)
        deletions = re.findall(r'-[0-9]+[ACGTNacgtn]+', mpiledup)
        
        mismatch_dict['ins'] = len(insertions)
        mismatch_dict['del'] = len(deletions)
        
        mpileupsnv = re.sub(r'\+[0-9]+[ACGTNacgtn]+|-[0-9]+[ACGTNacgtn]+', '', mpiledup)
        
        
        mismatch_dict['A'] = mpileupsnv.count('A')
        mismatch_dict['T'] = mpileupsnv.count('T')
        mismatch_dict['G'] = mpileupsnv.count('G')
        mismatch_dict['C'] = mpileupsnv.count('C')
        mismatchCount = mismatch_dict[alt_]
        mismatches += mismatchCount
        error_lists[i-1] = vaf(mismatchCount, depth)

    return round(float(mismatches/totalDepth), 3), totalDepth, mismatches, error_lists


def annotate_pon_vaf(df, PANEL_OF_NORMAL_BAMS, REFERENCE, new_column_name='pon_vaf', polyn=0):
    """Created a standalond function for starmap use"""
    df['pon_calc'] = df.apply(lambda x: calculate_pon_vaf(x['POS'], x['REF_ALT'], PANEL_OF_NORMAL_BAMS, REFERENCE), axis=1)
    df[[f'{new_column_name}_pon_vaf', f'{new_column_name}_pon_depth', f'{new_column_name}_alt_count', f'{new_column_name}_error_list']] = pd.DataFrame(df['pon_calc'].tolist(), index=df.index)
    if(polyn==1):
        df["polyn"]   = df.apply(lambda x: annotate_polyn(x['POS'], REFERENCE), axis=1)
    return df


def parse_pon_list(pon_list_tsv, exclusion_list):
    """parses a text file that has one line of full bam path for each panel of normal bams"""
    exclusion_list = [os.path.basename(i) for i in exclusion_list]
    pon_list = []
    with open(pon_list_tsv, 'r') as f:
        for line in f:
            bam = line.strip()
            if not os.path.basename(bam) in exclusion_list:
                pon_list.append(bam)

    print(pon_list)
    return pon_list


def main():
    input_tsv, output_dir, pon_list_tsv, REFERENCE, ncore, new_column_name, exclusion_list, polyn = argument_parser()
    output_file_path = os.path.join(output_dir, re.sub(r'.tsv$|.txt$', '.pon.tsv', os.path.basename(input_tsv)))

    # parse panel of normal tsv file
    PANEL_OF_NORMAL_BAMS = parse_pon_list(pon_list_tsv, exclusion_list)

    # read twin comparison tsv file into datafrae
    df = pd.read_csv(input_tsv, delimiter='\t')

    # split df into multicore chunks
    df_split = np.array_split(df, ncore)

    # prepare multicore starmap arg list
    arg_list = [] 
    for df in df_split:
        arg_list.append((df, PANEL_OF_NORMAL_BAMS, REFERENCE, new_column_name, polyn))

    

    # annotate PON on the data fame
    # df_pon = df.apply(lambda x: calculate_pon_vaf(x['POS'], PANEL_OF_NORMAL_BAMS, REFERENCE), axis=1)
    # multicore usage
    with mp.Pool(ncore) as pool:
        result = pool.starmap(annotate_pon_vaf, arg_list)

    # combine the multicore runs
    df_combined = pd.concat(result)

    # sort the result by position
    df_combined.sort_values('POS')

    # write the final table
    df_combined.to_csv(output_file_path, sep='\t', index=False)


if __name__=='__main__':
    main()
