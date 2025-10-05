"""script to find common/discordant mutations in twin quad sample
Test with one quad
2019.03.13 cjyoon
2019.03.14 cjyoon

clip fraction changed to 20 and rerun
run for all 30x quads and also 60x 

2019.04.20 cjyoon rerun with the final updated output file using v0.1c
2019.04.22 cjyoon remove those with median_deltaNM >5 as FP
2019.06.27 cjyoon run same pipeline on Novaseq result for T266
2019.10.15 cjyoon incorpoating into pipeline
"""
import cyvcf2
import re
import os
import argparse

def argument_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--family_id', required=True, help='family id')
    parser.add_argument('-1', '--twin1_vcf', required=True, help='twin1 final vcf from pztrio')
    parser.add_argument('-2', '--twin2_vcf', required=True, help='twin2 final vcf from pztrio')
    parser.add_argument('-r1', '--twin1_rawvcf', required=True, help='twin1 initial raw varscan vcf from pztrio')
    parser.add_argument('-r2', '--twin2_rawvcf', required=True, help='twin2 initial raw varscan vcf from pztrio')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    args = vars(parser.parse_args())

    return args['twin1_vcf'], args['twin2_vcf'], args['twin1_rawvcf'], args['twin2_rawvcf'], args['output_dir'], args['family_id']


def get_variant_sites(vcf):
    sites = [] 
    for variant in cyvcf2.VCF(vcf):
        if variant.INFO['child_alt'] >= 5 and variant.INFO['median_deltaNM'] <=5:
            sites.append(f'{variant.CHROM}:{variant.POS}\t{variant.REF}_{variant.ALT[0]}')
    return sites

def is_indel(ref_alt):
    ref, alt = ref_alt.split('_')
    if len(ref) != len(alt):
        return True
    else:
        return False
    

def get_noninherited_vafs(vcf1, vcf2, raw_vcf1, raw_vcf2, family_id, output_prefix, output_dir):
    child1_noninherited = get_variant_sites(vcf1)
    child2_noninherited = get_variant_sites(vcf2)
    print(len(child1_noninherited), len(child2_noninherited))
    noninherited_sets = list(set(child1_noninherited).union(set(child2_noninherited)))
    print(len(noninherited_sets))
    len(set(child1_noninherited).intersection(set(child2_noninherited)))
    total = sorted(list(set(child1_noninherited).union(set(child2_noninherited))))
    in_child1_list = []
    in_child2_list = []
    for i in total:
        in_child1 = False
        in_child2 = False
        if i in child1_noninherited:
            in_child1 = True
        
        in_child1_list.append(in_child1)

        if i in child2_noninherited:
             in_child2 = True
        in_child2_list.append(in_child2)

    vaf1_list = [] 
    vaf2_list = []
    with open(os.path.join(output_dir, output_prefix + '.noninherited.txt'), 'w') as f:
        f.write(f'POS\tREF_ALT\tALT1\tDEPTH1\tVAF1\tALT2\tDEPTH2\tVAF2\tIS_INDEL\tFAMILY\n')
        for variant_string, in_child1_s, in_child2_s in zip(total, in_child1_list, in_child2_list):
            chrom, pos = variant_string.split('\t')[0].split(':')
            ref_alt = variant_string.split('\t')[1]
            is_indel_s = is_indel(ref_alt)
            for variant in cyvcf2.VCF(raw_vcf1)(f'{chrom}:{pos}-{pos}'):
                try:
                    depth1 = variant.INFO['child_depth']
                    alt1 = variant.INFO['child_alt']
                    vaf1 = round(float(alt1/depth1), 4)

                except ZeroDivisionError:
                    vaf1 = -1
                    alt1, depth1 = 1, 1
            for variant in cyvcf2.VCF(raw_vcf2)(f'{chrom}:{pos}-{pos}'):
                try:
                    alt2 = variant.INFO['child_alt']
                    depth2 = variant.INFO['child_depth']
                    vaf2 = round(float(alt2/depth2), 4)    
                except ZeroDivisionError:
                    vaf2 = -1

            if vaf1 >=0 and vaf2 >= 0:
                f.write(f'{variant_string}\t{alt1}\t{depth1}\t{vaf1}\t{alt2}\t{depth2}\t{vaf2}\t{is_indel_s}\t{family_id}\n')
                vaf1_list.append(vaf1)
                vaf2_list.append(vaf2)
                             
    return vaf1_list, vaf2_list
    
    



def main():
    final_vcf1, final_vcf2, raw_vcf1, raw_vcf2, output_dir, family_id = argument_parser()
    get_noninherited_vafs(final_vcf1, final_vcf2, raw_vcf1, raw_vcf2, family_id, family_id, output_dir)


if __name__=='__main__':
    main()



