"""
Script to indentify potential post-zygotic somatic mutations in the child

Christopher J. Yoon (cjyoon@kaist.ac.kr) 

2019.06.27 cjyoon v0.1e increased depth threshold after finding out the pipeline was missing many true de novo mutations using too stringent depth filter. Few de novo variants were discovered in the Novaseq runs were simply filtered due to depth cut off. 
2019.10.15 cjyoon v0.1f changing output file prefix
2020.02.12 cjyoon v0.1f.1 added support for deepsequencing where depth threshold is ignored
2020.05.20 cjyoon v0.1f.2 added sample names from the three bam files that were used as inputs into vcf header
2021.02.25 cjyoon v1.0 calculating NM tag between variant reads and reference reads.
2021.05.15 cjyoon v1.1 read_pair_contains_variant fixed that caused error in proper phasing to het snp. Added phased readcount column in final tsv file
"""

import cyvcf2
import subprocess
import shlex
import pysam
import os
import sys
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import norm
import re
import argparse
import ast
import scipy.stats
import yaml
import datetime
import multiprocessing as mp

__version__ = "v1.1"

# config PATHS used in this script
path_config_path = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), 'path_config.yaml')

with open(path_config_path, 'r') as stream:
    try:
        path_config = (yaml.load(stream))
    except yaml.YAMLError as exc:
        print(exc)
        print('path_config cannot be parsed. Make sure it is in proper YAML format')
        sys.exit(1)

samtools = path_config['SAMTOOLS']
bgzip = path_config['BGZIP']
tabix = path_config['TABIX']
bcftools = path_config['BCFTOOLS']
cache_version_mapper = path_config['cache_version_mapper']
gnomad_genome_mapper = path_config['gnomad_genome_mapper']
VEP_PATH = path_config['VEP_PATH']
VEP_CACHE_DIR = path_config['VEP_CACHE_DIR']
VEP_ENVIRON = path_config['VEP_ENVIRON']
PERL5LIB = path_config['PER5LIB']


def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name


def natural_sort(l):
    def convert(text): return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key): return [convert(c)
                                   for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def get_chromosome_list(vcf):
    """get all chromosome list that appears in a vcf that is not an alternative contig"""
    chromosomes = set()
    for variant in cyvcf2.VCF(vcf):
        if re.search(r'^chr[0-9XY]+$|^[0-9XY]+$', variant.CHROM):
            chromosomes.add(variant.CHROM)

    return natural_sort(list(chromosomes))


def argument_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('--varscan_vcf', type=str, required=True,
                        help='Family level Varscan VCF file ({family_id}_varscan.vtdcn.atgc.dbsnpa.vcf.gz) that has been post-processed with varscan_parallel.smk')
    parser.add_argument('--father_bam', type=str, required=True,
                        help='BAM file of the father')
    parser.add_argument('--mother_bam', type=str, required=True,
                        help='BAM file of the mother')
    parser.add_argument('--child_bam', type=str, required=True,
                        help='BAM file of the child')
    parser.add_argument('-r', '--reference_fasta', required=True,
                        help='Reference fasta file')
    parser.add_argument('-a', '--assembly', type=str, choices=['GRCh37', 'GRCh38'], required=True,
                        help='Reference assembly')
    parser.add_argument('--father_id', type=str, required=False, default='',
                        help='ID for the father, DEFAULT=read from BAM header @RG SM: from father BAM file')
    parser.add_argument('--mother_id', type=str, required=False, default='',
                        help='ID for the mother, DEFAULT=read from BAM header @RG SM: from mother BAM file')
    parser.add_argument('--child_id', type=str, required=False, default='',
                        help='ID for the child (can be either the proband or the sibling, DEFAULT=read from BAM header @RG SM: from child BAM file')
    parser.add_argument('-o', '--output_dir', required=False,
                        default=os.getcwd(), help='Output directory, DEFAULT=current working directory')
    parser.add_argument('-t', '--min_read_count', default=2, type=int,
                        help='Read count cutoff to call a variant, DEFAULT=2')
    parser.add_argument('-v', '--min_vaf', default=0.01, type=float,
                        help='Read count cutoff to call a variant, DEFAULT=0.01')
    parser.add_argument('-g', '--germline_analysis', type=int, default=0, choices=[0, 1, 2],
                        help='0=only run postzygotic variant calling, 1=run both postzygotic and germline variant calling, 2=only run germline variant calling, DEFAULT=0')
    parser.add_argument('-n', '--ncore', type=int, default=1,
                        help='Number of cores to utilize for parallelization, DEFAULT=1')
    parser.add_argument('-m', '--mismatch_threshold', type=int, default=5,
                        help='Number of mismatches allowed in the variant read to be considered, DEFAULT=5')
    parser.add_argument('-p', '--prefix', required=False, default='', help='File name prefix for the result files')
    parser.add_argument('-d', '--deepseq', required=False, default=0, type=int, help='Ignoring depth threshold when using on deep sequencing results')

    args = vars(parser.parse_args())

    if args['min_vaf'] > 1:
        print('--min_vaf should be less than 1, exiting...')
        sys.exit(1)

    return args['varscan_vcf'], args['father_id'], args['mother_id'], args['child_id'], args['father_bam'], args['mother_bam'], args['child_bam'], args['reference_fasta'], args['output_dir'], args['min_read_count'], args['min_vaf'], args['assembly'], args['germline_analysis'], args['ncore'], args['mismatch_threshold'], args['prefix'], args['deepseq']


def get_id(individual_id, bamfile):
    if individual_id == '':
        return sampleNameBam(bamfile)
    else:
        return individual_id


def bgzip_tabix(vcf_file):
    """compress and index vcf file"""
    global bgzip
    global tabix

    indexed_vcf = vcf_file + '.gz'
    cmd = f'{bgzip} -f {vcf_file}'
    compress = subprocess.Popen(shlex.split(cmd))
    # print(cmd)
    compress.wait()

    cmd = f'{tabix} -f -p vcf {indexed_vcf}'
    index_run = subprocess.Popen(shlex.split(cmd))
    # print(cmd)
    index_run.wait()

    return indexed_vcf


def identify_trio_index(father_id, mother_id, child_id, vcf):
    """get the index positions of father mother child trio in the 
    given vcf file. 
    Infer this from the VCF header

    """
    vcf_handle = cyvcf2.VCF(vcf)

    try:
        # start indexing from 9 since first 9 is required vcf columns
        father_index = [i for i, j in enumerate(
            vcf_handle.samples) if j == father_id][0]
        mother_index = [i for i, j in enumerate(
            vcf_handle.samples) if j == mother_id][0]
        child_index = [i for i, j in enumerate(
            vcf_handle.samples) if j == child_id][0]
        return father_index, mother_index, child_index

    except IndexError:
        print(f'Given vcf file does not contain the family members')
        sys.exit(1)


def vaf(alt_count, depth):
    if depth > 0:
        return float(alt_count/depth)
    else:
        return 0


def count_int(count):
    """handle '.' in alt/ref counts"""
    if count == None or count < 0:
        return 0
    else:
        return count


def varscan_info_annotate(varscan_vcf, output_vcf, father_id, mother_id, child_id, command_string=None, timestamp=None, ncore=23):
    arg_list = []
    chrom_list = get_chromosome_list(varscan_vcf)
    for chrom in chrom_list:
        arg_list.append((varscan_vcf, output_vcf, father_id,
                         mother_id, child_id, command_string, timestamp, chrom))

    with mp.Pool(ncore) as pool:
        split_files = pool.starmap(varscan_info_annotate_chrom, arg_list)

    concat_cleanup_index(split_files, output_vcf)
    return output_vcf


def varscan_info_annotate_chrom(varscan_vcf, output_vcf, father_id, mother_id, child_id, command_string, timestamp, chrom):
    """Annotates Varscan VCF with allele counts for father, child, and mother
    and filter for those that can be potential parental mosaics. 

    Returns the filtered varscan_vcf_file
    """
    # GET GATK Filter list
    output_vcf = re.sub(r'.vcf.gz$', f'.{chrom}.vcf', output_vcf)

    father_index, mother_index, child_index = identify_trio_index(
        father_id, mother_id, child_id, varscan_vcf)
    # print(father_index, mother_index, child_index)
    varscan_putative_counts = 0
    varscan_putative_list = []
    varscan_vcf_handle = cyvcf2.VCF(varscan_vcf)
    if command_string != None:
        varscan_vcf_handle.add_to_header(f'##pztrioCommand={command_string}')

    if timestamp != None:
        varscan_vcf_handle.add_to_header(f'##pztrioTimestamp={timestamp}')

    varscan_vcf_handle.add_to_header(f'##pztrio={__version__}')

    varscan_vcf_handle.add_to_header(f'##child_id={child_id}')
    varscan_vcf_handle.add_to_header(f'##father_id={father_id}')
    varscan_vcf_handle.add_to_header(f'##mother_id={mother_id}')


    varscan_vcf_handle.add_info_to_header(
        {'ID': 'father_ref', 'Description': 'Ref count in father', 'Type': 'Integer', "Number": '1'})
    varscan_vcf_handle.add_info_to_header(
        {'ID': 'mother_ref', 'Description': 'Ref count in father', 'Type': 'Integer', "Number": '1'})
    varscan_vcf_handle.add_info_to_header(
        {'ID': 'child_ref', 'Description': 'Ref count in father', 'Type': 'Integer', "Number": '1'})

    varscan_vcf_handle.add_info_to_header(
        {'ID': 'father_alt', 'Description': 'Ref count in father', 'Type': 'Integer', "Number": '1'})
    varscan_vcf_handle.add_info_to_header(
        {'ID': 'mother_alt', 'Description': 'Ref count in father', 'Type': 'Integer', "Number": '1'})
    varscan_vcf_handle.add_info_to_header(
        {'ID': 'child_alt', 'Description': 'Ref count in father', 'Type': 'Integer', "Number": '1'})

    varscan_vcf_handle.add_info_to_header(
        {'ID': 'father_vaf', 'Description': 'VAF in father', 'Type': 'Float', 'Number': '1'})
    varscan_vcf_handle.add_info_to_header(
        {'ID': 'mother_vaf', 'Description': 'VAF in mother', 'Type': 'Float', 'Number': '1'})
    varscan_vcf_handle.add_info_to_header(
        {'ID': 'child_vaf', 'Description': 'VAF in child', 'Type': 'Float', 'Number': '1'})

    varscan_vcf_handle.add_info_to_header(
        {'ID': 'father_depth', 'Description': 'Depth in father', 'Type': 'Integer', 'Number': '1'})
    varscan_vcf_handle.add_info_to_header(
        {'ID': 'mother_depth', 'Description': 'Depth in mother ', 'Type': 'Integer', 'Number': '1'})
    varscan_vcf_handle.add_info_to_header(
        {'ID': 'child_depth', 'Description': 'Depth in child', 'Type': 'Integer', 'Number': '1'})

    varscan_filtered_output = cyvcf2.Writer(output_vcf, varscan_vcf_handle)
    for variant in varscan_vcf_handle(chrom):
        father_alt = count_int(variant.format('AD')[father_index][0])
        mother_alt = count_int(variant.format('AD')[mother_index][0])
        child_alt = count_int(variant.format('AD')[child_index][0])

        father_ref = count_int(variant.format('RD')[father_index][0])
        mother_ref = count_int(variant.format('RD')[mother_index][0])
        child_ref = count_int(variant.format('RD')[child_index][0])

        father_depth = int(father_ref) + int(father_alt)
        mother_depth = int(mother_ref) + int(mother_alt)
        child_depth = int(child_ref) + int(child_alt)

        father_vaf = vaf(father_alt, father_depth)
        mother_vaf = vaf(mother_alt, mother_depth)
        child_vaf = vaf(child_alt, child_depth)

        variant.INFO['father_ref'] = int(father_ref)
        variant.INFO['mother_ref'] = int(mother_ref)
        variant.INFO['child_ref'] = int(child_ref)
        variant.INFO['father_alt'] = int(father_alt)
        variant.INFO['mother_alt'] = int(mother_alt)
        variant.INFO['child_alt'] = int(child_alt)
        variant.INFO['father_vaf'] = float(father_vaf)
        variant.INFO['mother_vaf'] = float(mother_vaf)
        variant.INFO['child_vaf'] = float(child_vaf)
        variant.INFO['father_depth'] = father_depth
        variant.INFO['mother_depth'] = mother_depth
        variant.INFO['child_depth'] = child_depth

        varscan_filtered_output.write_record(variant)

    varscan_vcf_handle.close()
    varscan_filtered_output.close()

    return bgzip_tabix(output_vcf)


def is_dbsnp(variant_id):
    """looks at the ID column of VCF to see if is a member of dbsnp"""
    if re.search(r'^rs', variant_id):
        return 'T'
    else:
        return 'F'


def quality_annotate(input_vcf, output_vcf, bamfile, ncore=23):
    """annotates with variant site quality by reading in the child bam file"""

    chrom_list = get_chromosome_list(input_vcf)
    print(chrom_list)
    # infer read length
    readlength = infer_bam_readlength(bamfile)

    arg_list = []
    for chrom in chrom_list:
        arg_list.append((input_vcf, output_vcf, bamfile, chrom, readlength))

    with mp.Pool(ncore) as pool:
        quality_chrom = pool.starmap(quality_annotate_chrom, arg_list)

    concat_cleanup_index(quality_chrom, output_vcf)

    return output_vcf


def clip_length(cigarstring):
    """return length of clipped portion of cigar"""
    clip_lengths = ([int(re.sub(r'[S|H]', '', i))
                     for i in (re.findall(r'[0-9]+[H|S]', cigarstring))])
    return sum(clip_lengths)


def quality_annotate_chrom(input_vcf, output_vcf, bamfile, chrom, readlength):
    """per chrom split run of quality annotation"""
    output_vcf = re.sub('.vcf.gz', f'.{chrom}.vcf', output_vcf)
    vcf_handle = cyvcf2.VCF(input_vcf)
    vcf_handle.add_info_to_header(
        {'ID': 'MQMEAN', 'Description': 'MQ mean', 'Type': 'Float', "Number": 1})
    vcf_handle.add_info_to_header(
        {'ID': 'MQMEDIAN', 'Description': 'MQ median', 'Type': 'Float', "Number": 1})
    vcf_handle.add_info_to_header(
        {'ID': 'MQ0', 'Description': 'MQ0 count', 'Type': 'Integer', "Number": 1})
    vcf_handle.add_info_to_header(
        {'ID': 'CLIP', 'Description': 'CLIPPED read count', 'Type': 'Integer', "Number": 1})
    vcf_handle.add_info_to_header(
        {'ID': 'READCOUNT', 'Description': 'Read count', 'Type': 'Integer', "Number": 1})
    vcf_handle.add_info_to_header(
        {'ID': 'CLIP_FRACTION_REGION', 'Description': 'Clipped base fraction +- read length (150bp) window', 'Type': 'Float', "Number": 1})

    output_handle = cyvcf2.Writer(output_vcf, vcf_handle)
    bam = pysam.AlignmentFile(bamfile)

    for variant in vcf_handle(chrom):
        mq = []
        mq0_count = 0
        read_count = 0
        clip_read_count = 0
        dbsnp_member = is_dbsnp(str(variant.ID))
        for read in bam.fetch(variant.CHROM, variant.POS - 1, variant.POS):
            mq.append(read.mapping_quality)
            read_count += 1
            if read.mapping_quality == 0:
                mq0_count += 1

            if read.cigarstring != None:
                if re.search(r'S|H', read.cigarstring):
                    clip_read_count += 1
        try:
            mq0_fraction = round(float(mq0_count/read_count), 3)
        except ZeroDivisionError:
            mq0_fraction = -1
#            print('mq0_fraction_division_error' + str(variant))

        try:
            clip_read_fraction = round(float(clip_read_count/read_count), 3)
        except ZeroDivisionError:
            clip_read_fraction = -1
#            print('clipread fraction division_error' + str(variant))

        # get clip fraction +- readlength surrounding this variant region
        clipped_base = 0
        total_read_region = 0
        for read in bam.fetch(variant.CHROM, variant.POS-readlength, variant.POS+readlength):
            if read.cigarstring != None:
                clipped_base += (clip_length(read.cigarstring))
                total_read_region += 1
        total_bases = total_read_region * readlength
        clip_fraction_region = round(float(clipped_base/total_bases), 3)

        if read_count > 0:
            variant.INFO['MQMEAN'] = np.mean(mq)
            variant.INFO['MQMEDIAN'] = np.median(mq)
            variant.INFO['MQ0'] = mq0_count
            variant.INFO['CLIP'] = clip_read_count
            variant.INFO['READCOUNT'] = read_count
            variant.INFO['CLIP_FRACTION_REGION'] = clip_fraction_region

            output_handle.write_record(variant)

    output_handle.close()
    vcf_handle.close()
    return bgzip_tabix(output_vcf)


def quality_filter(input_vcf, output_vcf):
    """variant sites annotated with 'quality_annotate' is filtered with various info annotations
    changed clip fraction threshold to 20% 2019.03.13
    down to 15% 2019.03.14 9pm 
    up to 20% with regional clip filter 2019.03.21 8pm
    """
    output_vcf = re.sub('.vcf.gz', '.vcf', output_vcf)

    vcf_handle = cyvcf2.VCF(input_vcf)
    output_handle = cyvcf2.Writer(output_vcf, vcf_handle)
    total_count = 0
    after_filtering = 0
    for variant in vcf_handle:
        total_count += 1
        if variant.INFO['MQMEDIAN'] > 40 and variant.INFO['MQ0'] < 10:
            clip_fraction = float(
                variant.INFO['CLIP']/variant.INFO['READCOUNT'])
            if clip_fraction < 0.20 and variant.INFO['CLIP_FRACTION_REGION'] < 0.05:
                output_handle.write_record(variant)
                after_filtering += 1
    fraction = round(float(after_filtering/total_count), 3)
    print(f'# varaints input file: {total_count}')
    print(f'# variants after quality filter: {after_filtering}')
    print(f'Fraction remaining: {fraction}')
    output_handle.close()
    vcf_handle.close()

    return bgzip_tabix(output_vcf)


def somatic_count_vaf_filter(info_annotated_vcf, output_vcf, father_depth_cutoff, mother_depth_cutoff, child_depth_cutoff, min_read_count, min_vaf):
    """takes in INFO annotated VCF and filters variant based on their depth, and alt counts of the trio
    only looks at somatic (postzygotic or de novo)
    """

    output_vcf = re.sub(r'.gz$', '', output_vcf)
    vcf_handle = cyvcf2.VCF(info_annotated_vcf)
    filtered_output = cyvcf2.Writer(output_vcf, vcf_handle)

    for variant in vcf_handle:
        father_depth = variant.INFO['father_depth']
        mother_depth = variant.INFO['mother_depth']
        child_depth = variant.INFO['child_depth']
        if father_depth > 10 and mother_depth > 10 and child_depth > 10:
            if father_depth < father_depth_cutoff and mother_depth < mother_depth_cutoff and child_depth < child_depth_cutoff:
                if variant.INFO['child_alt'] >= min_read_count and variant.INFO['child_vaf'] > min_vaf and variant.INFO['father_vaf'] < 0.001 and variant.INFO['mother_vaf'] < 0.001:
                    filtered_output.write_record(variant)

    vcf_handle.close()
    filtered_output.close()

    return bgzip_tabix(output_vcf)


def germline_count_vaf_filter(info_annotated_vcf, output_vcf, father_depth_cutoff, mother_depth_cutoff, child_depth_cutoff, min_read_count, min_vaf):
    """takes in INFO annotated VCF and filters variant based on their depth, and alt counts of the trio
    only looks at germline variant (has to be inherited from at least one parent)
    """

    output_vcf = re.sub(r'.gz$', '', output_vcf)

    vcf_handle = cyvcf2.VCF(info_annotated_vcf)
    vcf_handle.add_info_to_header(
        {'ID': 'GT', 'Description': 'HET if 0.3<VAF<0.9, HOM if VAF > 0.9', 'Type': 'Character', "Number": '1'})
    filtered_output = cyvcf2.Writer(output_vcf, vcf_handle)

    for variant in vcf_handle:
        father_depth = variant.INFO['father_depth']
        mother_depth = variant.INFO['mother_depth']
        child_depth = variant.INFO['child_depth']
        if father_depth > 10 and mother_depth > 10 and child_depth > 10:
            if father_depth < father_depth_cutoff and mother_depth < mother_depth_cutoff and child_depth < child_depth_cutoff:
                if variant.INFO['child_alt'] >= min_read_count and variant.INFO['child_vaf'] > min_vaf and (variant.INFO['father_vaf'] > min_vaf or variant.INFO['mother_vaf'] > min_vaf):
                    if variant.INFO['child_vaf'] < 0.9:
                        variant.INFO['GT'] = 'HET'
                    else:
                        variant.INFO['GT'] = 'HOM'

                    filtered_output.write_record(variant)

    vcf_handle.close()
    filtered_output.close()

    return bgzip_tabix(output_vcf)


def annotate_variant_kmer(input_vcf, output_vcf, bamfile, ncore=23):
    """annotate output vcf with the most commmon kmer sequence around the variant in input_vcf"""
    arg_list = []
    chrom_list = get_chromosome_list(input_vcf)
    for chrom in chrom_list:
        arg_list.append((input_vcf, output_vcf, bamfile, chrom))

    with mp.Pool(ncore) as pool:
        split_files = pool.starmap(annotate_variant_kmer_chrom, arg_list)

    print(split_files)
    concat_cleanup_index(split_files, output_vcf)

    return (output_vcf)


def cleanup_files(file_list):
    """remove files given as list"""
    file_string = ' '.join(file_list)
    cleanup_cmd = f'rm -rf {file_string}'
    cleanup_run = subprocess.Popen(shlex.split(cleanup_cmd))
    cleanup_run.wait()
    return 0


def concat_cleanup_index(vcf_files, output_vcf):
    """merge list of vcf files, cleanup input vcf files, and index the merged vcf file"""
    global bcftools
    global tabix
    index_files = [i + '.tbi' for i in vcf_files]
    file_string = ' '.join(natural_sort(vcf_files))
    index_file_string = ' '.join(index_files)
    concat_cmd = f'{bcftools} concat {file_string} -o {output_vcf} -O z -a '
    concat = subprocess.Popen(shlex.split(concat_cmd))
    concat.wait()
    # clean up intermediate split vcf files
    cleanup_files(vcf_files + index_files)

    # tabix
    tabix_cmd = f'{tabix} -p vcf {output_vcf}'
    tabix_run = subprocess.Popen(shlex.split(tabix_cmd))
    tabix_run.wait()
    return output_vcf

def median_nm(array):
    if len(array)==0:
        return np.nan
    else:
        return int(np.median(array))
    
def delta_nm(ref_nm_list, var_nm_list):
    """calculates the difference in median nm values
    if no reference read -> then just report variant NM values"""
    if len(var_nm_list) > 0:
        if len(ref_nm_list) > 0:
            median_delta_nm = median_nm(var_nm_list) - median_nm(ref_nm_list)
        else:
            median_delta_nm  = median_nm(var_nm_list)
    else:
        median_delta_nm = -100
        
    return median_delta_nm
            
    
def annotate_variant_kmer_chrom(input_vcf, output_vcf, bamfile, chrom):
# input_vcf = '/home/users/cjyoon/notebook/twin_quad/test/ST3023.varscan.info.vaf.qual.qfilter.kmer.vcf.gz'
# output_vcf = '/home/users/cjyoon/notebook/twin_quad/test/test.vcf.gz'
# bamfile = child_bam
# chrom = 'chr22'
    """annotate output vcf with the most commmon kmer sequence around the variant in input_vcf"""
    output_vcf = re.sub(r'.vcf.gz$', f'.{chrom}.vcf', output_vcf)
    vcf_handle = cyvcf2.VCF(input_vcf)
    vcf_handle.add_info_to_header(
        {'ID': 'variant_kmer', 'Description': 'Kmer around variant read', 'Type': 'Character', "Number": '1'})
    vcf_handle.add_info_to_header(
        {'ID': 'median_NM', 'Description': 'Median value of NM tag on variant read', 'Type': 'Integer', "Number": '1'})
    vcf_handle.add_info_to_header(
        {'ID': 'median_deltaNM', 'Description': 'Median value of the difference in NM tag on variant reads vs reference reads. Also adjust variant length, so if difference is only explained by the variant, then 0 is the expected value.', 'Type': 'Integer', "Number": '1'})
    output_vcf_handle = cyvcf2.Writer(output_vcf, vcf_handle)
    for variant in vcf_handle(chrom):
        try:
            varcount = 0
            reads = get_reads_in_position(
                chrom=variant.CHROM, variant_pos=variant.POS, bamfile=bamfile)
            var_reads = []
            ref_reads = []
            variant_length = max(len(variant.REF), len(variant.ALT[0])) -1 
            for read in reads:
                if read_contains_variant(read, variant_chrom=variant.CHROM, variant_pos=variant.POS, variant_ref=variant.REF, variant_alt=variant.ALT[0]) == True:
                    var_reads.append(read)
                else:
                    ref_reads.append(read)

        #     print(len(ref_reads))
        #     print(len(var_reads))
        #     print(var_reads)
            kmer_list = []
            nmtag_counts_var_reads = []
            nmtag_counts_ref_reads = []
            for read in var_reads:
                kmer = kmer_around_variant(read, variant_chrom=variant.CHROM, variant_pos=variant.POS,
                                           variant_ref=variant.REF, variant_alt=variant.ALT[0], kmer_length=7)
                kmer_list.append(kmer)
                nmtag_counts_var_reads.append(read.get_tag('NM'))

            for read in ref_reads:
                try:
                    nmtag_counts_ref_reads.append(read.get_tag('NM'))
                except:
                    """unmapped read do not have NM values"""
                    pass
    #         print(f'NM tag var reads {nmtag_counts_var_reads}')
    #         print(np.median(nmtag_counts_var_reads), np.mean(nmtag_counts_var_reads))
    #         print(f'NM tag ref reads {nmtag_counts_ref_reads}')
    #         print(np.median(nmtag_counts_ref_reads), np.mean(nmtag_counts_ref_reads))



            if len(nmtag_counts_var_reads) > 0:
                variant.INFO['median_NM'] = int(median_nm(nmtag_counts_var_reads))
                variant.INFO['median_deltaNM'] = delta_nm(nmtag_counts_ref_reads, nmtag_counts_var_reads) - variant_length
            else:
                variant.INFO['median_NM'] = -100
                variant.INFO['median_deltaNM'] = -100

            # print(Counter(kmer_list).most_common(1))
            try:
                most_common_kmer = Counter(kmer_list).most_common(1)[0]
                variant.INFO['variant_kmer'] = str(
                    most_common_kmer[0]) + '_' + str(most_common_kmer[1])
            except IndexError:
                variant.INFO['variant_kmer'] = '_0'

            output_vcf_handle.write_record(variant)
        except ValueError:
            print(variant)
            print(len(ref_reads))
            print(len(var_reads))
            print(var_reads)
            print(f'NM tag var reads {nmtag_counts_var_reads}')
            print(np.median(nmtag_counts_var_reads), np.mean(nmtag_counts_var_reads))
            print(f'NM tag ref reads {nmtag_counts_ref_reads}')
            print(np.median(nmtag_counts_ref_reads), np.mean(nmtag_counts_ref_reads))

    print(f'done {chrom}')
    output_vcf_handle.close()
    vcf_handle.close()

    return bgzip_tabix(output_vcf)



def count_variant_kmer_chrom(input_vcf, output_vcf, father_bam, mother_bam, child_bam, chromosome):
    output_vcf = re.sub(r'.vcf.gz$', f'.{chromosome}.vcf', output_vcf)
    vcf_handle = cyvcf2.VCF(input_vcf)
    for member in ['father', 'mother', 'child']:
        vcf_handle.add_info_to_header(
            {'ID': f'{member}_variant_kmer_count', 'Description': f'Kmer around variant read in {member}', 'Type': 'Integer', "Number": '1'})
    output_vcf_handle = cyvcf2.Writer(output_vcf, vcf_handle)

    for variant in vcf_handle(chromosome):
        variant_kmer = variant.INFO['variant_kmer'].split('_')[0]
        variant_kmer_count = int(variant.INFO['variant_kmer'].split('_')[1])
        if variant_kmer_count > 2:  # if most common variant_kmer count is 1, probably a very bad site, so don't bother
            for member, bam in zip(['father', 'mother', 'child'], [father_bam, mother_bam, child_bam]):
                reads = get_all_reads(
                    chrom=variant.CHROM, variant_pos=variant.POS, bamfile=bam, mq_threshold=-1)
                var_read_count = 0
                for read in reads:
                    if re.search(variant_kmer, read.query_sequence):
                        var_read_count += 1
                variant.INFO[f'{member}_variant_kmer_count'] = var_read_count
            output_vcf_handle.write_record(variant)

    vcf_handle.close()
    output_vcf_handle.close()

    return bgzip_tabix(output_vcf)


def parse_kmer_count_text(kmer_text):
    kmer_count_dict = dict()
    with open(kmer_text, 'r') as f:
        for line in f:
            variant, count = line.strip().split('\t')
            count = int(count)
            kmer_count_dict.update({variant: count})

    return kmer_count_dict


def count_variant_kmer(input_vcf, output_vcf, father_bam, mother_bam, child_bam, ncore):
    """annotate with occurence of variant kmer in all individuals in trio
    2019.03.12 cjyoon parallelized this step for faster processing"""
    global bcftools
    chrom_list = get_chromosome_list(input_vcf)
    arg_list = []
    unsorted_output = output_vcf + '.tmp.vcf.gz'
    for chrom in chrom_list:
        arg_list.append((input_vcf, unsorted_output,
                         father_bam, mother_bam, child_bam, chrom))

    # parallelize each kmer counting in individuals as there are independent
    with mp.Pool(ncore) as pool:
        counted_result = pool.starmap(count_variant_kmer_chrom, arg_list)

    concat_cleanup_index(counted_result, unsorted_output)

    sort_cmd = f'{bcftools} sort -O z -o {output_vcf} {unsorted_output}'
    sort = subprocess.Popen(shlex.split(sort_cmd))
    sort.wait()
    # remove unsorted file
    cleanup_files([unsorted_output, unsorted_output + '.tbi'])
    return output_vcf


def kmer_around_variant(read, variant_chrom, variant_pos, variant_ref, variant_alt, kmer_length):
    """for a given read  see if this read contains SNV variant or reference"""
    pos = variant_pos - 1
    is_variant_read = 0
    kmer = None

    if (read.get_reference_positions()[0] - pos) * (read.get_reference_positions()[-1] - pos) <= 0:
        position_index_in_read = read.get_reference_positions(
            full_length=True).index(pos)
        base = (read.query_sequence[position_index_in_read]).upper()
        base_quality = read.query_qualities[position_index_in_read]
        is_variant_read = True
        kmer = (read.query_sequence[position_index_in_read -
                                     kmer_length: position_index_in_read + kmer_length]).upper()

    else:
        # print('non informative read included somehow')
        raise ValueError

    return kmer


def final_kmer_filter(input_vcf, output_vcf, mismatch_threshold):
    """looks at the INFO annotated kmer counts of the trio and 
    only write those with 0 count in parents and more than x kmers in the child"""
    count = 0
    output_vcf = re.sub(r'.vcf.gz$', '.vcf', output_vcf)
    vcf_handle = cyvcf2.VCF(input_vcf)
    output_handle = cyvcf2.Writer(output_vcf, vcf_handle)

    for variant in vcf_handle:
        if int(variant.INFO['father_variant_kmer_count']) == 0 and int(variant.INFO['mother_variant_kmer_count']) == 0:
            if variant.INFO['median_deltaNM'] <= mismatch_threshold:
                # print(variant)
                output_handle.write_record(variant)

    vcf_handle.close()
    output_handle.close()

    return bgzip_tabix(output_vcf)


def parse_region_of_interest(region):
    """standardize many different formats of region formats and returns as tuple"""
    region = re.sub(r',', '', region)
    if region.startswith('chr'):
        region = re.sub(r'^chr', '', region)

    chrom, pos = region.split(':')
    if len(pos.split('-')) > 1:
        start, end = pos.split('-')
    else:
        start = int(pos)
        end = start
    return chrom, start, end


def position_overlap(pos1, pos2):
    """checks if two positions overlap"""
    chr1, startend1 = pos1.split(':')
    [start1, end1] = [int(i) for i in startend1.split('-')]
    if start1 == end1:
        end1 = start1 + 1
    chr2, startend2 = pos2.split(':')
    [start2, end2] = [int(i) for i in startend2.split('-')]
    if start2 == end2:
        end2 = start2 + 1

    if chr1 == chr2:
        min_pos = min(start1, start2)
        max_pos = max(end1, end2)
        total_length = max_pos - min_pos
        pos1_length = end1 - start1
        pos2_length = end2 - start2
        if pos1_length + pos2_length >= total_length:
            return True
        else:
            return False
    else:
        return False


def outputfile(vcffile, output_dir):
    if vcffile.endswith('vcf.gz'):
        outputfile = os.path.join(output_dir, os.path.basename(
            re.sub(string=vcffile, pattern=r'vcf.gz$', repl='vep.vcf')))
    elif vcffile.endswith('vcf'):
        outputfile = os.path.join(output_dir, os.path.basename(
            re.sub(string=vcffile, pattern=r'vcf$', repl='vep.vcf')))
    else:
        print('Input must be a file that has a suffix of vcf or vcf.gz')
        raise ValueError

    return outputfile


def run_vep(input_vcf, temp_annotated_vcf, assembly_version, cache_version, gnomad_genome_vcf, display=True):
    """ if you want single annotation then use --per_gene option"""
    global VEP_PATH
    global VEP_CACHE_DIR
    global VEP_ENVIORN
    global PERL5LIB
    # DEFINE PERL5LIB PATH for other people's use
    os.environ['PERL5LIB'] = PERL5LIB
    os.environ['PATH'] = VEP_ENVIRON + ':' + os.environ['PATH']

    # 2018.11.13 cjyoon
    CACHE_VER = f'{cache_version}'
    ASSEMBLY_VER = assembly_version
    cmd = f'{VEP_PATH} --dir {VEP_CACHE_DIR} --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1  --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length  --offline --no_progress --no_stats  --polyphen b  --regulatory  -i {input_vcf} -o {temp_annotated_vcf} --force_overwrite --nearest symbol  --cache_version {CACHE_VER} --assembly {ASSEMBLY_VER} -custom {gnomad_genome_vcf},gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_POPMAX,POPMAX'

    if display == True:
        print(cmd)

    vep_cmd = subprocess.Popen(shlex.split(cmd))
    vep_cmd.wait()
    return 0


def find_canonical_index(vcf):
    """find which index CANONICAL annotation occurs by looking up the header info"""
    global bcftools
    vep_csq = re.sub(r'">', '', [i for i in subprocess.check_output(shlex.split(f'{bcftools} view -h {vcf}')).decode('utf-8').split('\n') if re.search(r'ID=CSQ', i)][0].split('Format: ')[1]).split('|')
    canonical_index = vep_csq.index('CANONICAL')
    return canonical_index


def find_canonical_annotation(vep_annotation_string, canonical_index):
    """VEP annotates with many alternative transcripts as well as canonical transcript
    this function finds the canonical transcript within vep_annotation_string.
    If there is no canonical transcript, which is usually the case fore intergenic,
    will just report the first annotation.
    """
    annotations = vep_annotation_string.split(',')
    return_status = 0
    for annotation in annotations:
        CANONICAL = annotation.split('|')[canonical_index]  # CANONICAL
        if CANONICAL == 'YES':
            return_status = 1
            return annotation

    if return_status == 0:
        return vep_annotation_string.split(',')[0]


def get_canonical_annotation(temp_annotated_vcf, annotated_vcf):
    """get only the canonical annotation into the CSQ INFO field"""
    canonical_index = find_canonical_index(temp_annotated_vcf)
    vcf_handle = cyvcf2.VCF(temp_annotated_vcf)
    vcf_writer = cyvcf2.Writer(annotated_vcf, vcf_handle)
    print('writing ' + annotated_vcf)
    for variant in vcf_handle:
        canonical_annotation = find_canonical_annotation(
            variant.INFO['CSQ'], canonical_index)
        variant.INFO['CSQ'] = canonical_annotation
        vcf_writer.write_record(variant)

    vcf_writer.close()

    return 0


def vcf_extract(input_vcf, region):
    global bcftools
    output_vcf = re.sub(r'.vcf.gz$|.vcf$', f'.{region}.vcf.gz', input_vcf)
    extract_cmd = f'{bcftools} view -O z -o {output_vcf} {input_vcf} {region}'
    extract_run = subprocess.Popen(shlex.split(extract_cmd))
    extract_run.wait()
    return output_vcf


def split_vcf_chrom(input_vcf, ncore=1):
    """given a vcf file, splits into per chromosome vcf"""
    chrom_list = get_chromosome_list(input_vcf)
    arg_list = []
    for chrom in chrom_list:
        arg_list.append((input_vcf, chrom))

    with mp.Pool(ncore) as pool:
        split_vcfs = pool.starmap(vcf_extract, arg_list)

    return split_vcfs


def vep_annotate(input_vcf, output_dir, genome_assembly, ncore):
    chrom_list = get_chromosome_list(input_vcf)
    split_vcf_list = split_vcf_chrom(input_vcf, ncore)
    output_vcf = re.sub(r'.vcf.gz$|.vcf$', '.vep.vcf.gz', input_vcf)
    arg_list = []
    for split_vcf in split_vcf_list:
        arg_list.append((split_vcf, output_dir, genome_assembly))

    with mp.Pool(ncore) as pool:
        vep_annotated_vcf_list = pool.starmap(vep_annotate_chrom, arg_list)

    concat_cleanup_index(vep_annotated_vcf_list, output_vcf)
    cleanup_files(split_vcf_list)

    return output_vcf


def vep_annotate_chrom(input_vcf, output_dir, genome_assembly):
    cache_version = cache_version_mapper[genome_assembly]
    gnomad_genome_vcf = gnomad_genome_mapper[genome_assembly]

    # prepare annotated output file path
    annotated_vcf = outputfile(input_vcf, output_dir)
    temp_annotated_vcf = annotated_vcf + '.temp.vcf'

    # run VEP
    run_vep(input_vcf, temp_annotated_vcf, genome_assembly,
            cache_version, gnomad_genome_vcf)

    # re-write with only canonical variants
    get_canonical_annotation(temp_annotated_vcf, annotated_vcf)
    # cleanup temp file
    cleanup_files([temp_annotated_vcf])
    # bgzip and tabix output vcf
    return bgzip_tabix(annotated_vcf)


def gnomad_filter(vcf, output_vcf, threshold):
    """filter a given vcf file with a gnomad population filter"""
    global bcftools
    output_vcf = re.sub(r'.gz$', '', output_vcf)

    vep_csq = re.sub(r'">', '', [i for i in subprocess.check_output(shlex.split(f'{bcftools} view -h {vcf}')).decode('utf-8').split('\n') if re.search(r'ID=CSQ', i)][0].split('Format: ')[1]).split('|')
    eas_af_index = vep_csq.index('gnomADg_AF_EAS')
    af_max_index = vep_csq.index('gnomADg_AF_POPMAX')

    vcf_handle = cyvcf2.VCF(vcf)
    output_vcf_handle = cyvcf2.Writer(output_vcf, vcf_handle)
    for variant in vcf_handle:
        csq = variant.INFO['CSQ'].split('|')

        # sometimes multiple freq annotated, in this case take the first one
        eas_af = (csq[eas_af_index].split('&')[0])
        af_max = (csq[af_max_index].split('&')[0])

        # if allele fractio not known, then to retain the variant assign 0 AF
        eas_af = float(eas_af) if eas_af not in ['', '.'] else 0
        af_max = float(af_max) if af_max not in ['', '.'] else 0

        if af_max < threshold:
            output_vcf_handle.write_record(variant)

    vcf_handle.close()
    output_vcf_handle.close()

    return bgzip_tabix(output_vcf)


def find_mate(bam_path, read):
    """finds the mate of read"""
    try:
        bamfile = pysam.AlignmentFile(bam_path)
        if not read.is_unmapped:
            if read.is_paired:
                for mate in bamfile.fetch(read.reference_name, read.next_reference_start, read.next_reference_start + 1):
                    if mate.query_name == read.query_name and mate.reference_name == read.reference_name:
                        return mate
        else:
            return None
    except:
        raise Exception


def get_informative_read_pairs(chrom, variant_pos, het_snp, bamfile):
    """retrieves informative read pair that contains both the variant_position and het_snp nearby"""
    informative_pair_count = 0
    informative_pairs = set()
    min_pos = min(variant_pos, het_snp)
    max_pos = max(variant_pos, het_snp)

    bam = pysam.AlignmentFile(bamfile)
    for read in bam.fetch(chrom, variant_pos-1, variant_pos+1):
        #         distance_to_snp = variant_pos-read.pos
        if read.is_paired and read.is_proper_pair:
            paired_read = find_mate(bamfile, read)
            # added mapping quality threshold 2019.02.07
            if paired_read != None and read.mapping_quality > 20 and paired_read.mapping_quality > 20:
                var_overlap = False
                snp_overlap = False

                for aread in [read, paired_read]:
                    # -1 to adjust for 0-based pysam
                    for pos_type, pos in zip(['var', 'snp'], [variant_pos - 1, het_snp - 1]):
                        if (aread.get_reference_positions()[0] - pos) * (aread.get_reference_positions()[-1] - pos) <= 0:
                            if pos_type == 'var':
                                var_overlap = True

                            if pos_type == 'snp':
                                snp_overlap = True

                    # if both variant and nearby snp is covered by this read pair
                if var_overlap and snp_overlap:
                    informative_pair_count += 1
#                             print(f'{aread} overlaps')

                    # to prevent writing the same pair twice, always record read1 first then read 2
                    if read.is_read1:
                        informative_pairs.add((read, paired_read))
                    else:
                        informative_pairs.add((paired_read, read))
                else:
                    pass  # print('not informative')

    return informative_pairs


def find_all_mate(bam_path, read):
    """finds the mate of read, EVEN THOSE THAT ARE ON DIFFERNT CHROMOSOMES
    Used for Kmer filtering
    """
    try:
        bamfile = pysam.AlignmentFile(bam_path)
        if not read.is_unmapped:
            if read.is_paired:
                for mate in bamfile.fetch(read.next_reference_name, read.next_reference_start, read.next_reference_start + 1):
                    return mate
        else:
            return None
    except:
        raise Exception


def get_reads_in_position(chrom, variant_pos, bamfile, mq_threshold=0):
    """Get all reads that align to that specific position
    DOES NOT retrieve mate pairs"""
    bam = pysam.AlignmentFile(bamfile)
    reads = set()
    for read in bam.fetch(chrom, variant_pos-1, variant_pos+1):
        if read.mapping_quality >= mq_threshold:
            reads.add(read)
    return list(reads)


def get_all_reads(chrom, variant_pos, bamfile, mq_threshold=-1):
    """retrieves all reads at a position
    different from get_informative_read_pairs in that it will retrieve ALL reads that map to position, even those that are not 
    properly mapped or mapping quality is poor. and also will NOT return as pairs
    Used for Kmer filtering
    """
    read_count = 0
    reads = set()

    bam = pysam.AlignmentFile(bamfile)
    for read in bam.fetch(chrom, variant_pos-1, variant_pos+1):
        #         distance_to_snp = variant_pos-read.pos
        if read.mapping_quality > mq_threshold:
                reads.add(read)
    return reads


def get_haplotype_linked_to_snp(alt_allele, haplotype):
    """haplotype is in Counter object
    find which Het SNP postzygotic mutation is linked to
    and haplotype dictionary of allele and counts 
    """
    pz_haplo = ''
    linked_snp = ''
    for haplo, count in haplotype.items():
        if haplo[0] == alt_allele:
            pz_haplo = haplo
            linked_snp = pz_haplo[-1]

    postzygotic_haplotype = dict()
    wt_haplotype = dict()
    supporting_read_counts = 0
    if linked_snp == '':
        postzygotic_haplotype = dict()
        wt_haplotype = dict()
    else:

        for haplo, count in haplotype.items():
            if haplo[1] == linked_snp:
                if haplo[0] == alt_allele:
                    postzygotic_haplotype.update({haplo: count})
                    supporting_read_counts += count
                else:
                    wt_haplotype.update({haplo: count})
                    supporting_read_counts += count

    return postzygotic_haplotype, wt_haplotype, supporting_read_counts


def get_haplotype_from_informative_pairs(informative_pairs, variant_pos, snp_pos):
    """from informative read pair containing variant_pos and snp_pos in paired read, 
    identify the base haplotypes

    also lookes at base quality threshold
    """
    haplo_list = []
    for r1, r2 in informative_pairs:
        haplo = ''
        # default is True, if any of the mutation/phasable snp is below threshold 20, do not append to list
        basequality_good = True
        for pos in [i-1 for i in [variant_pos, snp_pos]]:  # i-1 to adjust to 0-based pysam
            if (r1.get_reference_positions()[0] - pos) * (r1.get_reference_positions()[-1] - pos) <= 0:
                try:
                    position_index_in_read = r1.get_reference_positions(
                        full_length=True).index(pos)

                    base = (r1.query_sequence[position_index_in_read]).upper()
                    base_quality_r1 = r1.query_qualities[position_index_in_read]
                    if base_quality_r1 < 20:
                        basequality_good = False
                except:
                    base = 'D'  # if not present -> deletion

                haplo += (base)

            elif (r2.get_reference_positions()[0] - pos) * (r2.get_reference_positions()[-1] - pos) <= 0:
                try:
                    position_index_in_read_2 = r2.get_reference_positions(
                        full_length=True).index(pos)
                    base = (
                        r2.query_sequence[position_index_in_read_2]).upper()
                    base_quality_r2 = r2.query_qualities[position_index_in_read_2]
                    if base_quality_r2 < 20:
                        basequality_good = False

                except:
                    base = 'D'  # if not present -> deletion

                haplo += (base)
            else:
                # indicates nei
                raise ValueError

        if basequality_good == True:
            haplo_list.append(haplo)
        else:
            # base quality of either the potential postzygotic mutation or phasable SNP
            # do not meet the threshold 20
            pass

    return Counter(haplo_list)


def basecount_to_vaf(a_count, c_count, g_count, t_count):
    total = a_count + c_count + g_count + t_count
    vaf = float(
        sorted([a_count, c_count, g_count, t_count], reverse=True)[1])/total
    return vaf


def basecount_to_depth(a_count, c_count, g_count, t_count):
    total = a_count + c_count + g_count + t_count
    return total


def get_reference(variant_chrom, variant_pos, reference):
    ref = pysam.FastaFile(reference)
    base = ref.fetch(variant_chrom, variant_pos-1, variant_pos)
    return base


def get_alt_base(a_count, c_count, g_count, t_count, ref_base):
    """find the most frequent alternate base that is not ref base from count data"""
    basecounts = dict({'A': a_count, 'C': c_count,
                       'G': g_count, 'T': t_count, 'N': 0.5})
    try:
        basecounts.pop(ref_base)
    except:
        print(f'{ref_base} not ATGC')

    return (max(basecounts.keys(), key=(lambda key: basecounts[key])))


def is_het_snp_nearby(variant_chrom, variant_pos, bamfile, reference, depth_cutoff=100):
    """searches within 1000bp to see if there are any phasable SNPs near the variant position in the given bam file"""
    bam = pysam.AlignmentFile(bamfile, reference_filename=reference)
    het_snp = []
    search_window = 1000  # bp
    # identify het if  0.3 ~ 0.7 VAF nearby and depth > 15
    start = variant_pos - search_window
    end = variant_pos + search_window

    a, c, g, t = pysam.AlignmentFile(bamfile, reference_filename=reference).count_coverage(
        variant_chrom, start, end, quality_threshold=20)
    counts_df = pd.DataFrame({'A': a, 'C': c, 'G': g, 'T': t})
    counts_df['depth'] = counts_df.apply(
        lambda x: basecount_to_depth(x['A'], x['C'], x['G'], x['T']), axis=1)
    counts_df['vaf'] = counts_df.apply(
        lambda x: basecount_to_vaf(x['A'], x['C'], x['G'], x['T']), axis=1)
    counts_df['pos'] = counts_df.index + start + 1

    counts_df['ref'] = counts_df.apply(lambda x: get_reference(
        variant_chrom, x['pos'], reference), axis=1)
    counts_df['alt'] = counts_df.apply(lambda x: get_alt_base(
        x['A'], x['C'], x['G'], x['T'], x['ref']), axis=1)

    het_candidates = counts_df.loc[(counts_df['vaf'] >= 0.25) & (counts_df['vaf'] <= 0.75) & (
        counts_df['depth'] > 15) & (counts_df['depth'] < depth_cutoff) & (counts_df['pos'] != variant_pos)]

    het_candidates_string = []
    # convert to 1-based coordinate for the vcf file
    for pos, ref, alt in zip(het_candidates['pos'], het_candidates['ref'], het_candidates['alt']):
        het_candidates_string.append(f'{pos}_{ref}_{alt}')

    # for downstream purposes assign -1 for those without phasable snp
    if len(het_candidates_string) == 0:
        het_candidates_string = [-1]
    return het_candidates_string


def depth_distribution(vcf, father_id, mother_id, child_id, outputfile):
    """Infer the distribution of depth from a given VCF file
    Could do both het/homo, but here only getting mean/std for het since the results are comparable to each other
    2019.03.14 excluding those dbsnps in Y chromosome since this will skew males to lower depth and greater standard deviation
    """
    father_index, mother_index, child_index = identify_trio_index(
        father_id, mother_id, child_id, vcf)

    depth_dist = {key: [] for key in [father_index, mother_index, child_index]}
    depth_stats = []

    with open(outputfile, 'w') as f:
        for individual, individual_id, individual_index in zip(['father', 'mother', 'child'], [father_id, mother_id, child_id], [father_index, mother_index, child_index]):
            print(individual, individual_id, individual_index)
            for variant in cyvcf2.VCF(vcf):
                if re.search(r'rs', str(variant.ID)) and not re.search(r'Y', str(variant.CHROM)):
                    depth = int(variant.format('DP')[individual_index])
                    if depth < 100:
                        if (variant.genotypes[individual_index][0:2]) in [[0, 1], [1, 0]]:
                            depth_dist[individual_index].append(depth)

            if len(depth_dist[individual_index]) == 0:
                print(f'{vcf} is probably not annotated with common dbSNP ID. Check if your reference FASTA chromosome representation matches that of the dbSNP data.')
                print('Exiting...')
                sys.exit(1)

            mean = np.mean(depth_dist[individual_index])
            sd = np.std(depth_dist[individual_index])
            cutoff = mean + 1.96 * sd
            depth_stats.append((mean, sd, cutoff))
            print(f'{individual}\t{individual_id}\t{mean}\t{sd}\t{cutoff}\n')
            f.write(f'{individual}\t{individual_id}\t{mean}\t{sd}\t{cutoff}\n')

    return depth_stats


def parse_depth_stat(depth_stat):
    """parses depth stats file
    v0.1e increased threshold to twice the mean depth
    """
    cutoff_dict = dict()
    mean_depth_dict = dict()
    with open(depth_stat, 'r') as f:
        for line in f:
            member, member_id, mean_depth, sd_depth, cutoff_value = line.strip().split('\t')
            print(member)
            cutoff_value = float(cutoff_value)
            mean_depth = float(mean_depth)
            cutoff_dict.update({member:  cutoff_value})
            mean_depth_dict.update({member: mean_depth})
    return (mean_depth_dict['father']*2, mean_depth_dict['mother']*2, mean_depth_dict['child']*2), (mean_depth_dict['father'], mean_depth_dict['mother'], mean_depth_dict['child'])


def phase_parental_snp(input_vcf, output_vcf, father_bam, mother_bam, reference):
    """takes a look at the input vcf and searches for phasable snp in either father or the mother"""

    output_vcf = re.sub(r'.gz$', '', output_vcf)
    vcf_handle = cyvcf2.VCF(input_vcf)
    vcf_handle.add_info_to_header(
        {'ID': 'mosaic_parent', 'Description': 'Parent who can be a mosaic with respect to the variant', 'Type': 'Character', "Number": '1'})
    vcf_handle.add_info_to_header(
        {'ID': 'phasable_snp', 'Description': 'Position of potentially phasable snps within 1000bp', 'Type': 'Character', "Number": '1'})
    output_handle = cyvcf2.Writer(output_vcf, vcf_handle)
    for variant in vcf_handle:
        if variant.INFO['father_vaf'] == 0 and variant.INFO['mother_vaf'] > 0:
            bamfile = mother_bam
            variant.INFO['mosaic_parent'] = 'mother'
        elif variant.INFO['mother_vaf'] == 0 and variant.INFO['father_vaf'] > 0:
            bamfile = father_bam
            variant.INFO['mosaic_parent'] = 'father'
        else:
            print('VCF file is not appropriately filtered. Variants should have at least one of the parent with VAF=0')
            raise ValueError

        variant.INFO['phasable_snp'] = ','.join([str(i) for i in is_het_snp_nearby(
            variant.CHROM, variant.POS, bamfile, reference)])
        output_handle.write_record(variant)

    output_handle.close()
    vcf_handle.close()

    return bgzip_tabix(output_vcf)


def phase_child_snp_chrom(input_vcf, output_vcf, child_bam, reference, child_depth_cutoff, chrom):
    """takes a look at the input vcf and searches for phasable snp in the child"""
    output_vcf = re.sub(r'.vcf.gz$', f'.{chrom}.vcf', output_vcf)
    vcf_handle = cyvcf2.VCF(input_vcf)
    vcf_handle.add_info_to_header(
        {'ID': 'mosaic', 'Description': 'Parent who can be a mosaic with respect to the variant', 'Type': 'Character', "Number": '1'})
    vcf_handle.add_info_to_header(
        {'ID': 'phasable_snp', 'Description': 'Position of potentially phasable snps within 1000bp', 'Type': 'Character', "Number": '1'})
    output_handle = cyvcf2.Writer(output_vcf, vcf_handle)
    for variant in vcf_handle(chrom):
        variant.INFO['mosaic'] = 'child'

        variant.INFO['phasable_snp'] = ','.join([str(i) for i in is_het_snp_nearby(
            variant.CHROM, variant.POS, child_bam, reference, child_depth_cutoff)])
        output_handle.write_record(variant)

    output_handle.close()
    vcf_handle.close()

    return bgzip_tabix(output_vcf)


def phase_child_snp(input_vcf, output_vcf, child_bam, reference, child_depth_cutoff, ncore=23):
    """takes a look at the input vcf and searches for phasable snp in the child"""
    chrom_list = get_chromosome_list(input_vcf)
    arg_list = []
    for chrom in chrom_list:
        arg_list.append((input_vcf, output_vcf, child_bam, reference, child_depth_cutoff, chrom))

    with mp.Pool(ncore) as pool:
        phase_chrom = pool.starmap(phase_child_snp_chrom, arg_list)

    concat_cleanup_index(phase_chrom, output_vcf)

    return output_vcf


def read_contains_snv_variant(read, variant_chrom, variant_pos, variant_ref, variant_alt):
    """for a given read pair see if this read contains SNV variant or reference"""
    pos = variant_pos - 1
    is_variant_read = 0

    try:
        if (read.get_reference_positions()[0] - pos) * (read.get_reference_positions()[-1] - pos) <= 0:
            position_index_in_read = read.get_reference_positions(
                full_length=True).index(pos)
            base = (read.query_sequence[position_index_in_read]).upper()
            base_quality = read.query_qualities[position_index_in_read]
            if base == variant_alt and base_quality >= 20:
                is_variant_read = True
            elif base == variant_ref and base_quality >= 20:
                is_variant_read = False
            else:
                is_variant_read = None
    #             print(base, base_quality, base_quality >= 20
        else:
            # print('non informative read included somehow')
            raise ValueError
    except:
        is_variant_read = 0  # contains indels, do not use in calculation

    return is_variant_read


def read_contains_insertion_variant(read, variant_chrom, variant_pos, variant_ref, variant_alt):
    """for a given read pair see if this read contains Insertion variant or reference"""
    pos = variant_pos - 1
    variant_length = len(variant_alt)
    is_variant_read = 0
    try:
        if (read.get_reference_positions()[0] - pos) * (read.get_reference_positions()[-1] - pos) <= 0:
            insertion_start_index = read.get_reference_positions(
                full_length=True).index(pos) + 1
            insertion_end_index = read.get_reference_positions(
                full_length=True).index(pos) + variant_length
            insertion_base_quality = read.query_qualities[insertion_start_index -
                                                           1:insertion_end_index]
            # only consider variants with insertion base qualities  greater than 20
            if all(i > 10 for i in insertion_base_quality):
                if read.get_reference_positions(full_length=True)[insertion_start_index:insertion_end_index] == [None] * (variant_length-1) and read.query_sequence[insertion_start_index-1:insertion_end_index] == variant_alt:
                    #                 print('read contains variant')
                    is_variant_read = True
                else:
                    is_variant_read = False

        else:
            # print('non informative read included somehow')
            raise ValueError
    except:
        # print('exception')
        is_variant_read = 0

    return is_variant_read


def read_contains_deletion_variant(read, variant_chrom, variant_pos, variant_ref, variant_alt):
    """for a given read see if this read contains Deletion variant or reference"""
    pos = variant_pos - 1
    # since deletion variant length is calculated with variant_ref, assumes variant_alt is length 1
    variant_length = len(variant_ref)
    try:
        if (read.get_reference_positions()[0] - pos) * (read.get_reference_positions()[-1] - pos) <= 0:
            deletion_start_index = read.get_reference_positions(
                full_length=True).index(pos)
#             print(read.get_reference_positions(full_length=True))
#             print(deletion_start_index)
#             try:
            if read.get_reference_positions(full_length=True)[deletion_start_index + 1] == variant_pos + variant_length - 1:
                is_variant_read = True
            else:
                is_variant_read = False
#             except:
#                 is_variant_read = 0
        else:
            # print('non informative read included somehow')
            raise ValueError

    except:
        is_variant_read = 0
    return is_variant_read


def read_contains_variant(read, variant_chrom, variant_pos, variant_ref, variant_alt):
    """goes through read covering a position of a variant
    tries to identify those reads that support the provided variant and those that do not

    variant is a cyvcf2 variant class and variant can be SNV, small insertion, or small deletion

    assumes that variants have been decomposed into single alleleic variants
    and also left normalized (ideally using vt)
    2019.02.20 cjyoon

    """
    variant_length = len(variant_alt)
    is_variant_read = False

    # variant is an SNV
    if len(variant_ref) == 1 and len(variant_alt) == 1:
        # print('SNV')
        is_variant_read = read_contains_snv_variant(
            read, variant_chrom, variant_pos, variant_ref, variant_alt)

    # variant is a INSERTION
    elif len(variant_ref) < len(variant_alt):
        # print('INSERTION')
        is_variant_read = read_contains_insertion_variant(
            read, variant_chrom, variant_pos, variant_ref, variant_alt)

    # variant is a DELETION
    # for deletion there is no base quality to check for
    elif len(variant_ref) > len(variant_alt):
        # print('DELETION')
        is_variant_read = read_contains_deletion_variant(
            read, variant_chrom, variant_pos, variant_ref, variant_alt)

    else:
        print(f'{variant_chrom}_{variant_ref}>{variant_alt} is probably not left normalized single allelic record')
        raise ValueError
    return is_variant_read

def read_pair_contains_variant(areadpair, variant_chrom, variant_pos, variant_ref, variant_alt):
    """given a read pair checks if this read pair contains the variant of interest"""
    
    read1_status = read_contains_variant(areadpair[0], variant_chrom, variant_pos, variant_ref, variant_alt)
    read2_status = read_contains_variant(areadpair[1], variant_chrom, variant_pos, variant_ref, variant_alt)
    
    variant_status = read1_status or read2_status
    if None in [read1_status, read2_status]:
        return None
    else:
        return variant_status
    
    return

def get_haplotype_linked_to_snp(alt_allele, haplotype):
    """haplotype is in Counter object
    find which Het SNP postzygotic mutation is linked to
    and haplotype dictionary of allele and counts 
    """
    pz_haplo = ''
    linked_snp = ''
    for haplo, count in haplotype.items():
        variant, snp = haplo.split('_')
        if variant == alt_allele:
            pz_haplo = haplo
            linked_snp = pz_haplo.split('_')[1]

    postzygotic_haplotype = dict()
    wt_haplotype = dict()
    supporting_read_counts = 0
    if linked_snp == '':
        postzygotic_haplotype = dict()
        wt_haplotype = dict()
    else:

        for haplo, count in haplotype.items():
            if haplo.split('_')[1] == linked_snp:
                if haplo.split('_')[0] == alt_allele:
                    postzygotic_haplotype.update({haplo: count})
                    supporting_read_counts += count
                else:
                    wt_haplotype.update({haplo: count})
                    supporting_read_counts += count

    return postzygotic_haplotype, wt_haplotype, supporting_read_counts

def get_read_ids(read_pair_list):
    read_ids = set()
    for readpair in read_pair_list:
        read, pair = readpair
        read_ids.add(read.query_name)
    return read_ids


def phase(input_vcf, output_vcf, family_member, input_bam, output_dir):
    """write out phased result into vcf for a given individual and informative reads into a bam file"""
    global samtools
    output_vcf = re.sub(r'.gz$', '', output_vcf)
    output_bam = os.path.join(
        output_dir, os.path.basename(input_bam) + '.phase.bam')
    sorted_output_bam = os.path.join(
        output_dir, os.path.basename(input_bam) + '.phase.sorted.bam')
    bam_handle = pysam.AlignmentFile(input_bam)
    output_handle = pysam.AlignmentFile(output_bam, "wb", template=bam_handle)
    skip_count = 0
    phase_count = 0
    informative_pairs_all = set()

    vcf_handle = cyvcf2.VCF(input_vcf)
    vcf_handle.add_info_to_header({'ID': f'PHASE_{family_member}', 'Description': f'Phasing haplotype info of {family_member}', 'Type': 'Character', "Number": '1'})
    vcf_output_handle = cyvcf2.Writer(output_vcf, vcf_handle)
    # annotate vcf for future filtering
    for variant in vcf_handle:
        if variant.INFO['phasable_snp'] != '-1':
            snp_list = variant.INFO['phasable_snp'].split(',')
            haplo_dict = {str(i): {'post': {}, 'wt': {}} for i in snp_list}

            for snp in snp_list:
                snp_pos, snp_ref, snp_alt = snp.split('_')

                reads = get_informative_read_pairs(
                    chrom=variant.CHROM, variant_pos=variant.POS, het_snp=int(snp_pos), bamfile=input_bam)
                read_ids = get_read_ids(reads)

                informative_pairs_all = informative_pairs_all.union(read_ids)
                if len(reads) > 0:
                    haplotype_list = []
                    for areadpair in reads:
                        query_name = areadpair[0].query_name
                        is_snp_read = read_pair_contains_variant(areadpair, variant_chrom=variant.CHROM, variant_pos=int(
                            snp_pos), variant_ref=snp_ref, variant_alt=snp_alt)
                        read1_start_pos = (areadpair[0].reference_start)
                        read1_end_pos = (areadpair[0].reference_end)
                        read2_start_pos = (areadpair[0].reference_start)
                        read2_end_pos = (areadpair[0].reference_end)

                        is_variant_read = read_pair_contains_variant(
                            areadpair, variant_chrom=variant.CHROM, variant_pos=variant.POS, variant_ref=variant.REF, variant_alt=variant.ALT[0])

                        if is_variant_read == True and is_snp_read == True:
                            haplotype_list.append(f'{variant.ALT[0]}_{snp_alt}')
                        elif is_variant_read == True and is_snp_read == False:
                            haplotype_list.append(f'{variant.ALT[0]}_{snp_ref}')
                        elif is_variant_read == False and is_snp_read == True:
                            haplotype_list.append(f'{variant.REF}_{snp_alt}')
                        elif is_variant_read == False and is_snp_read == False:
                            haplotype_list.append(f'{variant.REF}_{snp_ref}')
                        else:
                            # bad read , do not consider in haplotypes
                            pass

                    haplotype = (Counter(haplotype_list))
                    postzygotic_haplotype, wt_haplotype, supporting_read_counts = (
                        get_haplotype_linked_to_snp(variant.ALT[0], haplotype=haplotype))
                    haplo_dict[str(snp)]['post'] = postzygotic_haplotype
                    haplo_dict[str(snp)]['wt'] = wt_haplotype

            # remove spaces for vcf INFO format
            variant.INFO[f'PHASE_{family_member}'] = str(haplo_dict).replace(' ', '')
            vcf_output_handle.write_record(variant)
        else:
            # not phasable SNV, but still write the final output
            variant.INFO[f'PHASE_{family_member}'] = 'NA'
            vcf_output_handle.write_record(variant)
        

    # Write a new bam file with informative reads only for visual inspection later.
    # print(informative_pairs_all)
    # print(len(informative_pairs_all))
    
    # recreate vcf_handle
    vcf_handle = cyvcf2.VCF(input_vcf)
    
    # go through the bam files and write the reads for the read in that were phase informative
    bam_handle = pysam.AlignmentFile(input_bam)
    for read in bam_handle:
        if read.query_name in informative_pairs_all:
            output_handle.write(read)

    bam_handle.close()
    output_handle.close()
    vcf_handle.close()
    vcf_output_handle.close()

    # sort and index output bam file
    cmd = f'{samtools} sort -O BAM -o {sorted_output_bam} {output_bam}'
    print(cmd)
    sort = subprocess.Popen(shlex.split(cmd))
    sort.wait()

    cmd = f'{samtools} index {sorted_output_bam}'
    indexing = subprocess.Popen(shlex.split(cmd))
    indexing.wait()
    cleanup_files([output_bam])

    return bgzip_tabix(output_vcf), sorted_output_bam


def parse_phase_info(phased_vcf, member):

    phase_dict = {}
    for variant in cyvcf2.VCF(phased_vcf):
        variant_string = f'{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}'
        phase_dict.update({variant_string: variant.INFO['PHASE_' + member]})

    return phase_dict


def combine_phased_trio(father_phased_vcf, mother_phased_vcf, child_phased_vcf, output_vcf):
    """combined split phased individuals into one single vcf"""
    output_vcf = re.sub(r'.vcf.gz$', '.vcf', output_vcf)
    vcf_handle = cyvcf2.VCF(child_phased_vcf)
    vcf_handle.add_info_to_header({'ID': f'PHASE_mother', 'Description': f'Phasing haplotype info of mother', 'Type': 'Character', "Number": '1'})
    vcf_handle.add_info_to_header({'ID': f'PHASE_father', 'Description': f'Phasing haplotype info of father', 'Type': 'Character', "Number": '1'})

    output_handle = cyvcf2.Writer(output_vcf, vcf_handle)

    mother_phase_dict = parse_phase_info(mother_phased_vcf, 'mother')
    father_phase_dict = parse_phase_info(father_phased_vcf, 'father')

    for variant in vcf_handle:
        variant_string = f'{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}'
        father_phase = father_phase_dict[variant_string]
        mother_phase = mother_phase_dict[variant_string]
        variant.INFO['PHASE_father'] = father_phase
        variant.INFO['PHASE_mother'] = mother_phase
        output_handle.write_record(variant)

    output_handle.close()
    vcf_handle.close()

    return bgzip_tabix(output_vcf)


def final_parental_postzygotic_phase_filter(input_vcf, output_vcf):
    """reads in input vcf annotated with phase information of trio
    filters for those that can be phased and output only those as final set"""
    postzygotic_counts = 0
    vcf_handle = cyvcf2.VCF(input_vcf)
    vcf_handle.add_info_to_header(
        {'ID': 'PHASE_SUMMARY', 'Description': 'Final phasing summary for parental mosaics', 'Type': 'Character', "Number": '1'})
    vcf_handle.add_info_to_header(
        {'ID': 'PHASE_SUPPORTING_READS', 'Description': 'Final phasing summary for parental mosaics', 'Type': 'Integer', "Number": '1'})
    output_vcf = re.sub(r'.vcf.gz$', '.vcf', output_vcf)

    output_vcf_handle = cyvcf2.Writer(output_vcf, vcf_handle)
    for variant in vcf_handle:
        if variant.INFO['mosaic_parent'] == 'mother':
            parent_haplotype = ast.literal_eval(variant.INFO['PHASE_mother'])
            wt_parent_haplotype = ast.literal_eval(
                variant.INFO['PHASE_father'])
        elif variant.INFO['mosaic_parent'] == 'father':
            parent_haplotype = ast.literal_eval(variant.INFO['PHASE_father'])
            wt_parent_haplotype = ast.literal_eval(
                variant.INFO['PHASE_mother'])

        else:
            print('INFO mosaic_parent has to be either father or mother')
            raise ValueError

        child_haplotype = ast.literal_eval(variant.INFO['PHASE_child'])
        is_postzygotic = False
        phasable_snp_with_read_support = []

        for phasable_snp in variant.INFO['phasable_snp'].split(','):
            if str(phasable_snp) in parent_haplotype.keys() and phasable_snp in child_haplotype.keys():
                if child_haplotype[str(phasable_snp)]['wt'] == dict() and parent_haplotype[str(phasable_snp)]['wt'] != dict() and parent_haplotype[str(phasable_snp)]['post'].keys() == child_haplotype[str(phasable_snp)]['post'].keys():
                    is_postzygotic = True
                    phasable_snp_with_read_support.append(phasable_snp)

        total_phase_support_reads = 0

        phase_summary_string = []
        if is_postzygotic == True:
            #             postzygotic_counts += 1
            for phasable_snp in phasable_snp_with_read_support:
                total_phase_support_reads += sum(
                    list(parent_haplotype[str(phasable_snp)]['wt'].values()))

                child_string = str(
                    child_haplotype[str(phasable_snp)]).replace(' ', '')
                mosaicparent_string = str(
                    parent_haplotype[str(phasable_snp)]).replace(' ', '')
                wtparent_string = str(
                    wt_parent_haplotype[str(phasable_snp)]).replace(' ', '')
                phase_summary_string.append(f'{phasable_snp}|child_{child_string}|mosaicparent_{mosaicparent_string}|wtparent_{wtparent_string}')

        if is_postzygotic == True:
            postzygotic_counts += 1
            variant.INFO['PHASE_SUMMARY'] = ('*'.join(phase_summary_string))
            variant.INFO['PHASE_SUPPORTING_READS'] = total_phase_support_reads
            output_vcf_handle.write_record(variant)

    print(f'Postzygotic mutation counts: {postzygotic_counts}')
    vcf_handle.close()
    output_vcf_handle.close()

    return bgzip_tabix(output_vcf)


def calculate_prob(input_vcf, output_vcf):
    """Calculate probability based on child, father, mother's read count"""
    output_vcf = re.sub(r'.vcf.gz$', '.vcf', output_vcf)
    vcf_handle = cyvcf2.VCF(input_vcf)
    vcf_handle.add_info_to_header(
        {'ID': 'pTRANSMISSION', 'Description': 'probability of observing more than child_depth reads given error rate estimated from the parents', 'Type': 'Float', "Number": '1'})
    vcf_handle.add_info_to_header(
        {'ID': 'pDENOVOHET', 'Description': 'probabilty of observing the read count given that the underlying site is a heterozygous site', 'Type': 'Float', "Number": '1'})

    output_vcf_handle = cyvcf2.Writer(output_vcf, vcf_handle)
    for variant in vcf_handle:
        if variant.INFO['father_depth'] > 0 and variant.INFO['mother_depth'] > 0:
            parent_var_count = max(
                variant.INFO['father_alt'] + variant.INFO['mother_alt'], 1)
            error_rate = parent_var_count / \
                (variant.INFO['father_depth'] + variant.INFO['mother_depth'])
            # probability of observing more than child_depth reads given error rate estimated from the parents
            transmission_prob = 1 - \
                scipy.stats.binom.cdf(
                    variant.INFO['child_alt'], variant.INFO['child_depth'], error_rate)
            # probabilty of observing the read count given that the underlying site is a heterozygous site
            # would expect 0.5x depth if true hetero, much lower if underlying site is a postzygotic site
            denovo_hetero_prob = scipy.stats.binom.cdf(
                variant.INFO['child_alt'], variant.INFO['child_depth'], 0.5)
            variant.INFO['pTRANSMISSION'] = transmission_prob
            variant.INFO['pDENOVOHET'] = denovo_hetero_prob
            output_vcf_handle.write_record(variant)

    output_vcf_handle.close()
    vcf_handle.close()

    return bgzip_tabix(output_vcf)

def pz_child_phase_summary(input_vcf, output_vcf):
    output_vcf = re.sub(r'.vcf.gz$', '.vcf', output_vcf)
    vcf_handle = cyvcf2.VCF(input_vcf)
    vcf_handle.add_info_to_header({'ID': f'PHASE_summary', 'Description': f'Summary of phasing info from family', 'Type': 'Character', "Number": '1'})
    vcf_handle.add_info_to_header({'ID': f'PHASE_summary2', 'Description': f'Summary of phasing info from family in read counts', 'Type': 'Character', "Number": '1'})
    output_handle = cyvcf2.Writer(output_vcf, vcf_handle)

    for variant in vcf_handle:
        if variant.INFO['phasable_snp'] != '-1':
            child_haplotype = ast.literal_eval(variant.INFO['PHASE_child'])
            father_haplotype = ast.literal_eval(variant.INFO['PHASE_father'])
            mother_haplotype = ast.literal_eval(variant.INFO['PHASE_mother'])
            phasable_snp_with_read_support = []
            phase_summary = []
            total_wt_phased = 0
            total_pz_phased = 0
            for phasable_snp in variant.INFO['phasable_snp'].split(','):
                print(phasable_snp)
                if str(phasable_snp) in child_haplotype.keys():
#                     print(child_haplotype[str(phasable_snp)])
#                     print(child_haplotype[str(phasable_snp)]['post'])

#                     print(child_haplotype[str(phasable_snp)]['post'])
#                     print(father_haplotype[str(phasable_snp)]['post'])
#                     print(mother_haplotype[str(phasable_snp)]['post'])

                    child_pz_status = (len(child_haplotype[str(phasable_snp)]['post']))
                    father_pz_status = (len(father_haplotype[str(phasable_snp)]['post']))
                    mother_pz_status = (len(mother_haplotype[str(phasable_snp)]['post']))

                    child_wt_linked_status = len(child_haplotype[str(phasable_snp)]['wt'])

#                     print(child_wt_linked_status, child_pz_status, father_pz_status, mother_pz_status)

#                     print(child_haplotype[str(phasable_snp)]['wt'])
#                     print(child_haplotype[str(phasable_snp)]['wt'] != dict())
#                     print(child_haplotype[str(phasable_snp)]['post'] != dict())
#                     print(father_haplotype[str(phasable_snp)]['post'] == dict() and mother_haplotype[str(phasable_snp)]['post'] == dict())
                    if child_pz_status==1 and child_wt_linked_status==0 and father_pz_status==0 and mother_pz_status==0:
                        is_postzygotic = True
                        phasable_snp_with_read_support.append(phasable_snp)

                        child_string = str(
                            child_haplotype[str(phasable_snp)]).replace(' ', '')
                        father_string = str(
                            father_haplotype[str(phasable_snp)]).replace(' ', '')
                        mother_string = str(
                            mother_haplotype[str(phasable_snp)]).replace(' ', '')

                        phase_summary.append(f'{phasable_snp}|child_{child_string}|father_{father_string}|mother_{mother_string}')
            phase_summary_string = '*'.join(phase_summary)

            if len(phase_summary) == 0:
                variant.INFO['PHASE_summary'] = 'NA'
                variant.INFO['PHASE_summary2'] = 'NA'
            else:
                variant.INFO['PHASE_summary'] = '*'.join(phase_summary)
                print(variant.INFO['PHASE_summary'])
                for info in phase_summary_string.split('*'):
                    snp, child_info, father_info, mother_info = info.split('|')
                    child_info_dict = (ast.literal_eval(re.sub(r'child_', '', child_info)))
                    total_pz_phased += (sum(child_info_dict['post'].values()))
                    total_wt_phased += (sum(child_info_dict['wt'].values()))

                variant.INFO['PHASE_summary2'] = f'{str(total_pz_phased)}/{str(total_wt_phased)}'

        else:
            variant.INFO['PHASE_summary'] = 'NA'
            variant.INFO['PHASE_summary2'] = 'NA'


        output_handle.write_record(variant)
    vcf_handle.close()
    output_handle.close()
    return bgzip_tabix(output_vcf)

def get_vep_annotation_index(vcf):
    """finds which index corresnponds to which VEP annotation"""
    vcf_handle = cyvcf2.VCF(vcf)
    vep_index_dict = {}
    for header in vcf_handle.header_iter():
        try:
            if (header['ID']) == 'CSQ':
                vep_annotation_list = header['Description'].split(
                    'Format: ')[1].strip('"').split('|')
                for vep_index, vep_index_item in enumerate(vep_annotation_list):
                    vep_index_dict.update({vep_index_item: vep_index})
                return vep_index_dict
        except KeyError:
            pass


def variant_type(ref, alt):
    """determines whether the variant is SNV, insertion, or deletion by looking at ref and alt length"""
    if len(ref) == 1 and len(alt) == 1:
        return 'SNV'
    elif len(ref) > len(alt):
        return "DEL"
    elif len(ref) < len(alt):
        return "INS"
    else:
        print('variant not properly left normalized')
        raise ValueError


def is_phased(phase_string):
    if phase_string == 'NA':
        return 'F'
    else:
        return 'T'


def get_variant_info(variant, info):
    """retrieves INFO annotated in the variant class, 
    returns NA if INFO is not annotated"""
    try:
        return variant.INFO[info]
    except KeyError:
        return 'NA'


def write_final_table(input_vcf, output_table, child_id):
    """create a final table format of vcf that can be used for visual insepction and annotation"""
    vcf_handle = cyvcf2.VCF(input_vcf)
    vep_annotation_index = get_vep_annotation_index(input_vcf)
    with open(output_table, 'w') as f:
        # write header
        f.write(f'CHROM\tPOS\tREF\tALT\tvar_type\tchild_id\tfather_ref\tmother_ref\tchild_ref\tfather_alt\tmother_alt\tchild_alt\t')
        f.write(f'father_vaf\tmother_vaf\tchild_vaf\tfather_depth\tmother_depth\tchild_depth\t')
        f.write(f'MQMEAN\tMQMEDIAN\tMQ0\tCLIP\tREADCOUNT\tCLIP_FRACTION_REGION\tmedian_NM\t')
        f.write(f'pDENOVOHET\tpTRANSMISSION\tconsequence\tsymbol\thgvsp\t')
        f.write(f'gnomADg_AF_POPMAX\tPHASE_summary\tPHASE_rc\tPHASE_status\n')
        for variant in vcf_handle:
            vep_annotation = variant.INFO['CSQ'].split('|')
            consequence = vep_annotation[vep_annotation_index['Consequence']]
            symbol = vep_annotation[vep_annotation_index['SYMBOL']]
            hgvsp = vep_annotation[vep_annotation_index['HGVSp']]
            gnomADg_AF_POPMAX = vep_annotation[vep_annotation_index['gnomADg_AF_POPMAX']]
            var_type = variant_type(variant.REF, variant.ALT[0])
            phase_status = is_phased(
                get_variant_info(variant, 'PHASE_summary'))

            for_table = [variant.CHROM, variant.POS, variant.REF, variant.ALT[0], var_type, child_id,
                         get_variant_info(variant, 'father_ref'), get_variant_info(
                             variant, 'mother_ref'), get_variant_info(variant, 'child_ref'),
                         get_variant_info(variant, 'father_alt'), get_variant_info(
                             variant, 'mother_alt'), get_variant_info(variant, 'child_alt'),
                         get_variant_info(variant, 'father_vaf'), get_variant_info(
                             variant, 'mother_vaf'), get_variant_info(variant, 'child_vaf'),
                         get_variant_info(variant, 'father_depth'), get_variant_info(
                             variant, 'mother_depth'), get_variant_info(variant, 'child_depth'),
                         get_variant_info(variant, 'MQMEAN'), get_variant_info(
                             variant, 'MQMEDIAN'), get_variant_info(variant, 'MQ0'),
                         get_variant_info(variant, 'CLIP'), get_variant_info(
                             variant, 'READCOUNT'), get_variant_info(variant, 'CLIP_FRACTION_REGION'), get_variant_info(variant, 'median_NM'),
                         get_variant_info(variant, 'pDENOVOHET'), get_variant_info(
                             variant, 'pTRANSMISSION'), consequence, symbol, hgvsp,
                         gnomADg_AF_POPMAX, get_variant_info(variant, 'PHASE_summary'), get_variant_info(variant, 'PHASE_summary2'), phase_status]

            for_table = [str(i) for i in for_table]
            f.write('\t'.join(for_table) + '\n')

    f.close()
    vcf_handle.close()

    return output_table


def infer_bam_readlength(bam_path):
    """reads first 1,000,000 reads of a bam file to get read length of the bam file"""
    bam = pysam.AlignmentFile(bam_path)
    read_sample_count = 1000000
    read_length_vector = np.zeros(read_sample_count)
    i = 0
    for read in bam.fetch():
        read_length_vector[i] = (read.infer_read_length())
        i += 1
        if i == read_sample_count:
            break
        else:
            continue
    return int(np.nanmedian(read_length_vector))


def main():
    # timestamp
    start_time = datetime.datetime.now()
    timestamp = start_time.strftime("%Y-%m-%d %H:%M %Z")

    # parse arguments from command line
    global reference_fasta
    varscan_vcf, father_id, mother_id, child_id, father_bam, mother_bam, child_bam, reference_fasta, output_dir, min_read_count, min_vaf, assembly, germline_analysis, ncore, mismatch_threshold, prefix, deepseq = argument_parser()
    command_string = ' '.join(sys.argv)
    print(command_string)

    # create output directory if it does not already exist
    os.system(f'mkdir -p {output_dir}')

    # if no ID is given for father/mother/child ids, then read from the BAM header
    father_id = get_id(father_id, father_bam)
    mother_id = get_id(mother_id, mother_bam)
    child_id = get_id(child_id, child_bam)

    ####################################################################
    ####################################################################
    # 1. Parse and initial filter for Varscan/Freebayes for initial Child Mosaics
    if prefix=='':
        varscan_info_annotate_vcf = os.path.join(output_dir, f'{child_id}.varscan.info.vcf.gz')
    else:
        varscan_info_annotate_vcf = os.path.join(output_dir, f'{prefix}.varscan.info.vcf.gz')

    varscan_initial_filtered_vcf = re.sub(
        r'.vcf.gz$', '.vaf.vcf.gz', varscan_info_annotate_vcf)
    varscan_qual_annotated_vcf = re.sub(
        r'.vcf.gz$', '.qual.vcf.gz', varscan_initial_filtered_vcf)
    varscan_qual_filtered_vcf = re.sub(
        r'.vcf.gz$', '.qfilter.vcf.gz', varscan_qual_annotated_vcf)
    varscan_kmer_annotated_vcf = re.sub(
        r'.vcf.gz$', '.kmer.vcf.gz', varscan_qual_filtered_vcf)
    varscan_kmer_counted_vcf = re.sub(
        r'.vcf.gz$', '.counted.vcf.gz', varscan_kmer_annotated_vcf)
    varscan_kmer_filtered_vcf = re.sub(
        r'.vcf.gz$', '.kfilter.vcf.gz', varscan_kmer_counted_vcf)

    depth_stats = os.path.join(output_dir, f'{child_id}.depth.stats')

    # 1.1 Get necessary depth information from dbSNP annotations
    if not os.path.isfile(depth_stats) or os.path.getsize(depth_stats) == 0:
        (father_depth_mean, father_depth_sd, father_depth_cutoff), (mother_depth_mean, mother_depth_sd, mother_depth_cutoff), (child_depth_mean,
                                                                                                                               child_depth_sd, child_depth_cutoff) = depth_distribution(varscan_vcf, father_id, mother_id, child_id, depth_stats)

    (father_depth_cutoff, mother_depth_cutoff, child_depth_cutoff), (father_depth_mean, mother_depth_mean, child_depth_mean) = parse_depth_stat(
        depth_stats)


    # If Bam file is coming from deep sequencing, then remove depth cutoff by setting it to infinity
    if deepseq == 1:
        father_depth_cutoff, mother_depth_cutoff, child_depth_cutoff = np.inf, np.inf, np.inf

    print(1.1, father_depth_cutoff, mother_depth_cutoff, child_depth_cutoff)

    # 1.2 annotate VCF with INFO section with necessary VAF, REF, ALT counts for downstream
    if not os.path.isfile(varscan_info_annotate_vcf):
        varscan_info_annotate_vcf = varscan_info_annotate(
            varscan_vcf, varscan_info_annotate_vcf, father_id, mother_id, child_id, command_string, timestamp, ncore)

    print(1.2, varscan_info_annotate_vcf)

    #############################################################################
    #############################################################################
    # run germline analysis only if -g is set to 1 or 2
    if germline_analysis > 0:
        varscan_germline_initial_filtered_vcf = re.sub(
            r'.vcf.gz$', '.germline.vaf.vcf.gz', varscan_info_annotate_vcf)
        varscan_germline_qual_annotated_vcf = re.sub(
            r'.vcf.gz$', '.qual.vcf.gz', varscan_germline_initial_filtered_vcf)
        varscan_germline_qual_filtered_vcf = re.sub(
            r'.vcf.gz$', '.qfilter.vcf.gz', varscan_germline_qual_annotated_vcf)
        # Germline
        # 1.3.1  Filter annotated VCF with DEPTH, and VAFs
        if not os.path.isfile(varscan_germline_initial_filtered_vcf):
            varscan_germline_initial_filtered_vcf = germline_count_vaf_filter(
                varscan_info_annotate_vcf, varscan_germline_initial_filtered_vcf, father_depth_cutoff, mother_depth_cutoff, child_depth_cutoff, min_read_count, min_vaf=0.3)
        print(1.31, varscan_germline_initial_filtered_vcf)
        # Germline
        # 1.4.1 Quality Annotate VCF
        if not os.path.isfile(varscan_germline_qual_annotated_vcf):
            varscan_germline_qual_annotated_vcf = quality_annotate(
                varscan_germline_initial_filtered_vcf, varscan_germline_qual_annotated_vcf, bamfile=child_bam, ncore=ncore)
        print(1.41, varscan_germline_qual_annotated_vcf)
        # Germline
        # 1.5.1 Quality Filter VCF
        if not os.path.isfile(varscan_germline_qual_filtered_vcf):
            varscan_germline_qual_filtered_vcf = quality_filter(
                varscan_germline_qual_annotated_vcf, varscan_germline_qual_filtered_vcf)
        print(1.51, varscan_germline_qual_filtered_vcf)
        # 3.1.1. VEP annotate
        vep_germline_annotated_vcf = re.sub(
            r'.vcf.gz$', '.vep.vcf.gz', varscan_germline_qual_filtered_vcf)
        if not os.path.isfile(vep_germline_annotated_vcf):
            vep_germline_annotated_vcf = vep_annotate(
                varscan_germline_qual_filtered_vcf, output_dir, assembly, ncore)
        print(3.11, vep_germline_annotated_vcf)
        # 3.2.1 Go through VCF and retain only variants with population frequency < 0.01
        germline_population_filtered_vcf = re.sub(
            r'.vcf.gz$', '.popfilter.vcf.gz', vep_germline_annotated_vcf)
        if not os.path.isfile(germline_population_filtered_vcf):
            germline_population_filtered_vcf = gnomad_filter(
                vep_germline_annotated_vcf, germline_population_filtered_vcf, 0.01)
        print(3.21, germline_population_filtered_vcf)

        # 3.3 write table
        germline_table = re.sub(r'.vcf.gz$', '.tsv',
                                vep_germline_annotated_vcf)
        if not os.path.isfile(germline_table):
            germline_table = write_final_table(
                vep_germline_annotated_vcf, germline_table, child_id)
        print(3.31, germline_table)

    #############################################################################
    #############################################################################

    if germline_analysis == 2:
        sys.exit()

    # 1.3 Filter annotated VCF with DEPTH, and VAFs
    if not os.path.isfile(varscan_initial_filtered_vcf):
        varscan_initial_filtered_vcf = somatic_count_vaf_filter(
            varscan_info_annotate_vcf, varscan_initial_filtered_vcf, father_depth_cutoff, mother_depth_cutoff, child_depth_cutoff, min_read_count, min_vaf)
    print(1.3, varscan_initial_filtered_vcf)

    # 1.4 Quality Annotate VCF
    if not os.path.isfile(varscan_qual_annotated_vcf):
        varscan_qual_annotated_vcf = quality_annotate(
            varscan_initial_filtered_vcf, varscan_qual_annotated_vcf, bamfile=child_bam, ncore=ncore)
    print(1.4, varscan_qual_annotated_vcf)

    # 1.5 Quality Filter VCF
    if not os.path.isfile(varscan_qual_filtered_vcf):
        varscan_qual_filtered_vcf = quality_filter(
            varscan_qual_annotated_vcf, varscan_qual_filtered_vcf)
    print(1.5, varscan_qual_filtered_vcf)

    # 1.6 variant Kmer annotate VCF
    if not os.path.isfile(varscan_kmer_annotated_vcf):
        varscan_kmer_annotated_vcf = annotate_variant_kmer(
            varscan_qual_filtered_vcf, varscan_kmer_annotated_vcf, child_bam, ncore)

    print(1.6, varscan_kmer_annotated_vcf)

    # 1.7 Kmer count all individuals VCF
    if not os.path.isfile(varscan_kmer_counted_vcf):
        varscan_kmer_counted_vcf = count_variant_kmer(
            varscan_kmer_annotated_vcf, varscan_kmer_counted_vcf, father_bam, mother_bam, child_bam, ncore)

    print(1.7, varscan_kmer_filtered_vcf)

    # 1.8 Kmer filter VCF
    if not os.path.isfile(varscan_kmer_filtered_vcf):
        varscan_kmer_filtered_vcf = final_kmer_filter(
            varscan_kmer_counted_vcf, varscan_kmer_filtered_vcf, mismatch_threshold)

    print(1.8, varscan_kmer_filtered_vcf)

    # 2.1 VEP annotation
    vep_annotated_vcf = re.sub(
        r'.vcf.gz$', '.vep.vcf.gz', varscan_kmer_filtered_vcf)
    print(vep_annotated_vcf)

    if not os.path.isfile(vep_annotated_vcf):
        vep_annotated_vcf = vep_annotate(
            varscan_kmer_filtered_vcf, output_dir, assembly, ncore)

    print(2.1, vep_annotated_vcf)

    # 2.2 Go through VCF and retain only variants with population frequency < 0.01
    population_filtered_vcf = re.sub(
        r'.vcf.gz$', '.popfilter.vcf.gz', vep_annotated_vcf)
    if not os.path.isfile(population_filtered_vcf):
        population_filtered_vcf = gnomad_filter(
            vep_annotated_vcf, population_filtered_vcf, 0.01)
    print(2.2, population_filtered_vcf)

    # 3.3 Calculate transmission probability and de novo heterozygous probability to identify postzygotic mutation
    pcalculated_vcf = re.sub(
        r'.vcf.gz$', '.pcal.vcf.gz', population_filtered_vcf)
    if not os.path.isfile(pcalculated_vcf):
        pcalculated_vcf = calculate_prob(
            population_filtered_vcf, pcalculated_vcf)
    print(2.3, pcalculated_vcf)

    ####################################################################
    ####################################################################
    # 3.1 Annotate with phasable snps
    snp_phased_vcf = re.sub(r'.vcf.gz$', '.phased.vcf.gz', pcalculated_vcf)

    if not os.path.isfile(snp_phased_vcf):
        snp_phased_vcf = phase_child_snp(
            pcalculated_vcf, snp_phased_vcf, child_bam, reference_fasta, child_depth_mean *2, ncore)
    print(3.1, snp_phased_vcf)

    # 3.2~3.4. Annotate with haplotype info and write out bam file
    # sequentially annotate same vcf with father, mother, then child haplotypes
    # each annotation will take approximately 1 hour
    father_haplotype_vcf = re.sub('.vcf.gz$', '.f.vcf.gz', snp_phased_vcf)
    mother_haplotype_vcf = re.sub('.vcf.gz$', '.m.vcf.gz', snp_phased_vcf)
    child_haplotype_vcf = re.sub('.vcf.gz$', '.c.vcf.gz', snp_phased_vcf)
    trio_haplotype_vcf = re.sub('.vcf.gz$', '.fmc.vcf.gz', snp_phased_vcf)

    if not (os.path.isfile(father_haplotype_vcf) and os.path.isfile(mother_haplotype_vcf) and os.path.isfile(child_haplotype_vcf)):
        with mp.Pool(min(3, ncore)) as pool:
            (father_haplotype_vcf, father_phased_bam), (mother_haplotype_vcf, mother_phased_bam), (child_haplotype_vcf, child_phased_bam) = pool.starmap(phase, [(snp_phased_vcf, father_haplotype_vcf, 'father', father_bam, output_dir),
                                                                                                                                                                 (snp_phased_vcf, mother_haplotype_vcf, 'mother', mother_bam, output_dir),
                                                                                                                                                                 (snp_phased_vcf, child_haplotype_vcf, 'child', child_bam, output_dir)])

    print(3.2, father_haplotype_vcf, mother_haplotype_vcf, child_haplotype_vcf)
    combine_phased_trio(father_haplotype_vcf, mother_haplotype_vcf,
                        child_haplotype_vcf, trio_haplotype_vcf)

    ####################################################################
    ####################################################################
    # 4. Now go through the final vcf file and select for those variants with sufficient phasable read supporting
    final_child_postzygotic_vcf = re.sub(
        r'.vcf.gz$', '.child_pz.vcf.gz', trio_haplotype_vcf)
    if not os.path.isfile(final_child_postzygotic_vcf):
        final_parental_postzygotic_vcf = pz_child_phase_summary(
            trio_haplotype_vcf, final_child_postzygotic_vcf)
    print(4, final_child_postzygotic_vcf)

    # 5 Create a final table from final vcf file
    final_table = re.sub(r'.vcf.gz$', '.tsv', final_child_postzygotic_vcf)
    if not os.path.isfile(final_table):
        final_table = write_final_table(
            final_child_postzygotic_vcf, final_table, child_id)
    print(5, final_table)

    ####################################################################
    # Get total time to run the code
    end_time = datetime.datetime.now()
    time_delta = str(end_time - start_time)
    print(f'Total time to run the code: {time_delta}')

    sys.exit(0)


if __name__ == '__main__':
    main()
