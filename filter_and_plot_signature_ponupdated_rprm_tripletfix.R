# 2021.06.17 cjyoon. After meeting with ysju, needed to change the minimal filter condition >=3 reads instead of starting with VAF > 0.1 
# 2021.06.19 cjyoon. Other families were having erros on the MZtwin_annotator function. But as the function had too many arguments and long function, it was difficult to debug 
# therefore splitting the annotator function into separate filter function, and then combining results into one column like mztwin_annotator output
# 2021.06.22 redo with new output files
# 2021.06.28 also do it with trinucleotide context added output file.
# 2021.07.01 modify script so that filtering still works for buccal only samples which do not have columns associated with blood data
# 2021.07.02 having 3 minimal read threshold made too many calls in samples without blood. Maybe need to adjust the number of read counts. 
# draw signature based on different read counts?
# see how mutation counts change?
# 2021.07.15 cjyoon fixed parent, child pon value to updated pon list, 
# 2021.07.15  constitutional_twin prob1 < 0.05 
# 2021.07.15 also included ST501 
# 2021.07.22 cjyoon. After manual review of missing true variants, realized RP needs to be removed or changed for those with few reads. 
# also need to change ALT_COUNT threshold to add ALT3_buccal for triplet -> separate function


library(tidyverse)
library(ggplot2)
library(gridExtra)

# read in mutation signature plot script 
source('~/notebook/twin_quad/mutation_signature_plot.R')
result_dir = '/home/users/cjyoon/Projects/twin_quad/analysis/buccal_batch4_v1.1/twin_comparison/'
family_ids = gsub(pattern = '.noninherited.baminfo.pon.pon.pon.tri.tsv', replacement = '', list.files(path=result_dir, pattern='.noninherited.baminfo.pon.pon.pon.tri.tsv'))

delta_nm <- function(nm1, nm2){
  if(is.na(nm2) & is.na(nm1)){
    return(0)
  }else if(is.na(nm2) & !(is.na(nm1))){
    return(nm1)
  }else if(is.na(nm1) & !(is.na(nm2))){
    return(-nm2)
  }else{
    return(nm1-nm2)
  }
}
delta_nm_vec = Vectorize(delta_nm)

constitutional_twin <- function(prob1, prob2){
  if(prob1 >= 0.01 & prob2 >= 0.01){
    return('denovo')
  }else if(prob1 < 0.01 & prob2 >= 0.01){
    return('one')
  }else if(prob1 >= 0.01 & prob2 < 0.01){
    return('one')
  }else{
    return('none')
  }
}
constitutional_twin_vec = Vectorize(constitutional_twin)


constitutional_triplet <- function(prob1, prob2, prob3){
  all_probs = c(prob1, prob2, prob3)
  clonal_counts = length(which(all_probs > 0.01))
  if(clonal_counts==3) {
    return('denovo')
  }else if(clonal_counts==2){
    return('two')
  }else if(clonal_counts==1){
    return('one')
  }else{
    return('none')
  }
}
constitutional_triplet_vec = Vectorize(constitutional_triplet)



clip_fraction <- function(clipped_read, alt_read){
  if(alt_read==0){
    return(0)
  }
  else if(is.na(alt_read)){
    return(0)
  }else if(is.na(clipped_read)){
    return(0)
  }else{
    return(clipped_read/alt_read)
  }
}
clip_fraction_vec<-Vectorize(clip_fraction)

test_na_arguments <- function(...){
  if(length(which(is.na(...)))>0){
    return(TRUE) # NA is present
  }else{
    return(FALSE) # NA is not present
  }
}

parent_pon_annotator <- function(parent_pon_vaf){
  if(!(parent_pon_vaf <0.006)){
    return('parentPON;')
  }else{
    return('')
  }
}
parent_pon_annotator_vec = Vectorize(parent_pon_annotator)

child_pon_annotator <- function(child_pon_vaf){
  if(!(child_pon_vaf <0.007)){
    return('childPON;')
  }else{
    return('')
  }
}
child_pon_annotator_vec = Vectorize(child_pon_annotator)

sb_annotator<-function(SB1_buccal, SB2_buccal){
  if(is.na(SB1_buccal) | is.na(SB2_buccal)){
    return('')
  }else if(!(SB1_buccal > 0.005 & SB2_buccal > 0.005)==T){
    return('SB;')
  }
  else{
    return('')
  }
}
sb_annotator_vec = Vectorize(sb_annotator)

rp_annotator<-function(mean_RP_alt1_buccal, mean_RP_alt2_buccal){
  if((!(is.na(mean_RP_alt1_buccal) | (mean_RP_alt1_buccal<0.9 & mean_RP_alt1_buccal>0.1))==T|
      !(is.na(mean_RP_alt2_buccal) | (mean_RP_alt2_buccal<0.9 & mean_RP_alt2_buccal>0.1))==T)==T){
    return("RP;")
  }else{
    return ('')
  }
}
rp_annotator_vec = Vectorize(rp_annotator)
blood_support_annotator<-function(ALT1_blood, ALT2_blood){
  if(is.na(ALT1_blood)|is.na(ALT2_blood)){
    return('')
  }else{
    if(!(ALT1_blood + ALT2_blood >=1)==T){
      return("NotInBlood;")
    }else{
      return('')
    }
  }
}
blood_support_annotator_vec = Vectorize(blood_support_annotator)

clip_annotator<-function(CLIPPED_ALT1_buccal, CLIPPED_ALT2_buccal, ALT1_buccal, ALT2_buccal){
  buccal1_clip_fraction = clip_fraction_vec(CLIPPED_ALT1_buccal, ALT1_buccal)
  buccal2_clip_fraction = clip_fraction_vec(CLIPPED_ALT2_buccal, ALT2_buccal)
  if(((buccal1_clip_fraction > 0.3) | (buccal2_clip_fraction > 0.3))==T){
    return('clip;')
  }else{
    return('')
  }
}
clip_annotator_vec = Vectorize(clip_annotator)

repeat_annotator<-function(up_polyn, down_polyn){
  if (!(up_polyn<10 & down_polyn<10)==T){
    return("Repeat;")
  }else{
    return('')
  }
}
repeat_annotator_vec = Vectorize(repeat_annotator)

alt_annotator<-function(ALT1_buccal, ALT2_buccal){
  MIN_ALT_THRESHOLD = 3
  if(is.na(ALT1_buccal) | is.na(ALT2_buccal)){
    return('ALT;')
  }else if(!(ALT1_buccal+ALT2_buccal>=MIN_ALT_THRESHOLD)==T){
    return('ALT;')
  }else{
    return('')
  }
}
alt_annotator_vec = Vectorize(alt_annotator)

alt_annotator_triplet<-function(ALT1_buccal, ALT2_buccal, ALT3_buccal){
  MIN_ALT_THRESHOLD = 3
  if(is.na(ALT1_buccal) | is.na(ALT2_buccal) | is.na(ALT3_buccal)){
    return('ALT;')
  }else if(!(ALT1_buccal+ALT2_buccal + ALT3_buccal>=MIN_ALT_THRESHOLD)==T){
    return('ALT;')
  }else{
    return('')
  }
}
alt_annotator_triplet_vec = Vectorize(alt_annotator_triplet)


summarize_filter<-function(parent_pon_filter, child_pon_filter, sb_filter, blood_filter, clip_filter, repeat_filter, alt_filter){
  result = paste0(parent_pon_filter, child_pon_filter, sb_filter, blood_filter, clip_filter, repeat_filter, alt_filter, sep = '')
  if(result==''){
    result='PASS'
  }
  return(result)
}
summarize_filter_vec = Vectorize(summarize_filter)


plot_family <- function(family_id){
  
  print(family_id)
  #   for(MIN_ALT_THRESHOLD in 5:10){
  df_pon_path = paste0(result_dir, family_id, '.noninherited.baminfo.pon.pon.pon.tri.tsv')
  df_pon = read_delim(df_pon_path, delim='\t')
  df_pon = df_pon %>% separate(polyn, c('ref_upstream_base', 'upstream_polyN', 'ref_downstream_base', 'downstream_polyN'), sep=':')
  df_pon$upstream_polyN = as.integer(df_pon$upstream_polyN)
  df_pon$downstream_polyN = as.integer(df_pon$downstream_polyN)
  
  # check and see if this family has twin blood data or not
  
  if(!('ALT1_blood' %in% colnames(df_pon))){
    # if this family has no blood related columns, then add all those columns with NA value
    blood_columns = c('DEPTH1_blood','REF1_blood','ALT1_blood','CLIPPED_ALT1_blood','OTHER1_blood','VAF1_blood','RDF1_blood','RDR1_blood','ADF1_blood','ADR1_blood','SB1_blood','mean_NM_ref1_blood','mean_NM_alt1_blood','median_NM_ref1_blood','median_NM_alt1_blood','mean_RP_alt1_blood','ID_context_ref1_blood','ID_context_alt1_blood','ref_up_polyn1_blood','ref_down_polyn1_blood','DEPTH2_blood','REF2_blood','ALT2_blood','CLIPPED_ALT2_blood','OTHER2_blood','VAF2_blood','RDF2_blood','RDR2_blood','ADF2_blood','ADR2_blood','SB2_blood','mean_NM_ref2_blood','mean_NM_alt2_blood','median_NM_ref2_blood','median_NM_alt2_blood','mean_RP_alt2_blood','ID_context_ref2_blood','ID_context_alt2_blood','ref_up_polyn2_blood','ref_down_polyn2_blood')
    df_pon[blood_columns] = NA 
  }
  df_pon = df_pon %>% mutate(VAF1_hetprobcdf=pbinom(ALT1_buccal, DEPTH1_buccal, 0.5), VAF2_hetprobcdf=pbinom(ALT2_buccal, DEPTH2_buccal, 0.5)) %>% 
    mutate(dNM_buccal1 = delta_nm_vec(median_NM_alt1_buccal, median_NM_ref1_buccal)) %>% 
    mutate(dNM_buccal2 = delta_nm_vec(median_NM_alt2_buccal, median_NM_ref2_buccal)) %>% 
    mutate(parent_pon_filter = parent_pon_annotator_vec(parent_pon_vaf)) %>% 
    mutate(child_pon_filter = child_pon_annotator_vec(child_pon_vaf)) %>% 
    mutate(sb_filter = sb_annotator_vec(SB1_buccal, SB2_buccal)) %>% 
    # mutate(rp_filter = rp_annotator_vec(mean_RP_alt1_buccal, mean_RP_alt2_buccal)) %>% 
    mutate(blood_filter = blood_support_annotator_vec(ALT1_blood, ALT2_blood)) %>% 
    mutate(clip_filter = clip_annotator_vec(CLIPPED_ALT1_buccal, CLIPPED_ALT2_buccal, ALT1_buccal, ALT2_buccal)) %>% 
    mutate(repeat_filter = repeat_annotator_vec(up_polyn, down_polyn))
    
  # for triplet use different alt filter to account for the extra child 
  if(family_id=='ST501'){
    df_pon = df_pon %>% mutate(alt_filter = alt_annotator_triplet_vec(ALT1_buccal, ALT2_buccal, ALT3_buccal)) %>% 
      mutate(VAF3_hetprobcdf=pbinom(ALT3_buccal, DEPTH3_buccal, 0.5)) %>% 
      mutate(DENOVO=constitutional_triplet_vec(VAF1_hetprobcdf, VAF2_hetprobcdf, VAF3_hetprobcdf))
  }else{
    df_pon = df_pon %>% mutate(alt_filter = alt_annotator_vec(ALT1_buccal, ALT2_buccal)) %>% 
      mutate(DENOVO=constitutional_twin_vec(VAF1_hetprobcdf, VAF2_hetprobcdf))
  }  
  df_pon = df_pon %>% mutate(FILTER=summarize_filter_vec(parent_pon_filter, child_pon_filter, sb_filter, blood_filter, clip_filter, repeat_filter, alt_filter)) %>%
    select(-c(parent_pon_filter, child_pon_filter, sb_filter, blood_filter, clip_filter, repeat_filter, alt_filter))
  
  df_filtered = df_pon %>% filter(FILTER=='PASS') %>% mutate(family=family_id)
  count_summary = (df_filtered %>% group_by(DENOVO) %>% summarise(n=n()))
  count_summary$family = family_id
  
  pdf(paste0('~/Projects/twin_quad/analysis/metric_filter/', family_id, '_202100722_metricfilter.pdf'), width=7, height=7, useDingbats = F)
  p <- df_filtered %>%
    ggplot(aes(x=VAF1_buccal, y=VAF2_buccal)) + geom_point(aes(col=DENOVO), alpha=0.5, cex=5) + xlab('Twin1 VAF') + ylab('Twin2 VAF') + scale_color_manual(values=c('denovo'="black", 'none'="blue", 'one'="red", 'two'='orange')) +
    xlim(0, 1) + ylim(0, 1) + theme_classic(base_size=24)  + theme(legend.position = "none") +# ggtitle(paste0(family_id, ' PON + dNM filtered\n ALT1+ALT2>=', MIN_ALT_THRESHOLD)) +
    annotate("text", x = 0.5, y = 1, label=family_id, size=10)
  # grid.text(family_id, x=unit(1, 'npc') - unit(0.5, 'char'), y=unit(1, 'npc') -  unit(0.5, 'char'), just=c('right', 'top'), gp=gpar(fontsize=30))
  print(p, newpage=F)
  dev.off()
  
  # Write final filtered Dataframe
  output_df_path = paste0('~/Projects/twin_quad/analysis/metric_filter/', family_id, '_20210722_metricfilter.tsv')
  write_delim(df_filtered, output_df_path, delim='\t')
  
  # Also write the original with annotation for debugging and rescuing
  output_df_path_all_annotated = paste0('~/Projects/twin_quad/analysis/metric_filter/', family_id, '_20210722_metricfilter_annotated_all.tsv')
  write_delim(df_pon, output_df_path_all_annotated, delim='\t')
  
  
  return(list(counts = count_summary, variants = df_filtered))
}


# THESE BOTTOM CODES ARE COMMENTED OUT TO USE THIS SCRIPT AS SOURCE FOR SAMPLES THAT WERE SEQUENCED LATER. COMMENTED OUT ON OCT 17, 2021. cjyoon
# 
# 
# all_filtered_variants = data.frame()
# count_summary_all = data.frame()
# 
# for(family_id in family_ids){
#   result = plot_family(family_id)
#   count_summary_all = bind_rows(count_summary_all, result$counts)
#   all_filtered_variants = bind_rows(all_filtered_variants, result$variants)
#   
# }
# count_summary_all
# count_summary_all %>% spread(key=DENOVO, value=n) %>% filter(str_detect(family, 'ST[1|2|3|5]')) -> variant_count_stats
# variant_count_stats
# total_sites = sum(variant_count_stats$denovo, na.rm=T) + sum(variant_count_stats$one, na.rm=T) + sum(variant_count_stats$none, na.rm=T) + sum(variant_count_stats$two, na.rm=T) 
# 
# # Write all filtered variants into a single file
# output_df_all_path = paste0('~/Projects/twin_quad/analysis/metric_filter/all_snuh_20210722_metricfilter.tsv')
# write_delim(all_filtered_variants, output_df_all_path, delim='\t')
# 
# # extract signatures for the filtered variants
# # de novo
# df_96 = generate_empty_signature_df()
# 
# pdf('~/Projects/twin_quad/analysis/metric_filter/mutation_signature_denovo_20210722.pdf', useDingbats = F, width=14, height=6)
# all_filtered_variants %>% filter(str_detect(family, 'ST[1|2|3]')) %>% select(POS, REF_ALT, tricontext, DENOVO) %>% separate(tricontext, into=c('tri', 'sub'), sep='_') %>% filter(tri != 'indel') %>% group_by(DENOVO, sub, tri) %>% summarise(count=n()) %>% 
#   filter(DENOVO=='denovo') -> denovo_counts 
# denovo_counts %>% plot_mutation_signature
# denovo_counts = df_96 %>% left_join(denovo_counts) %>% mutate(count = replace_na(count, 0)) %>% select(sub, tri, count) 
# colnames(denovo_counts) = c('Substitution', 'Trinucleotide', 'Count')
# output_count_path = paste0('~/Projects/twin_quad/analysis/metric_filter/mutation_count_denovo_20210722.tsv')
# write_delim(denovo_counts, output_count_path, delim='\t')
# dev.off()
# 
# # one
# pdf('~/Projects/twin_quad/analysis/metric_filter/mutation_signature_one_20210722.pdf', useDingbats = F, width=14, height=6)
# all_filtered_variants %>% filter(str_detect(family, 'ST[1|2|3]')) %>% select(POS, REF_ALT, tricontext, DENOVO) %>% separate(tricontext, into=c('tri', 'sub'), sep='_') %>% filter(tri != 'indel') %>% group_by(DENOVO, sub, tri) %>% summarise(count=n()) %>% 
#   filter(DENOVO=='one') -> one_counts 
# one_counts %>% plot_mutation_signature
# 
# one_counts = df_96 %>% left_join(one_counts) %>% mutate(count = replace_na(count, 0)) %>% select(sub, tri, count) 
# colnames(one_counts) = c('Substitution', 'Trinucleotide', 'Count')
# output_count_path = paste0('~/Projects/twin_quad/analysis/metric_filter/mutation_count_one_20210722.tsv')
# write_delim(one_counts, output_count_path, delim='\t')
# dev.off()
# 
# 
# 
# # none
# pdf('~/Projects/twin_quad/analysis/metric_filter/mutation_signature_none_20210722.pdf', useDingbats = F, width=14, height=6)
# all_filtered_variants %>% filter(str_detect(family, 'ST[1|2|3]')) %>% select(POS, REF_ALT, tricontext, DENOVO) %>% separate(tricontext, into=c('tri', 'sub'), sep='_') %>% filter(tri != 'indel') %>% group_by(DENOVO, sub, tri) %>% summarise(count=n()) %>% 
#   filter(DENOVO=='none') -> none_counts
# none_counts %>% plot_mutation_signature
# none_counts = df_96 %>% left_join(none_counts) %>% mutate(count = replace_na(count, 0)) %>% select(sub, tri, count) 
# colnames(none_counts) = c('Substitution', 'Trinucleotide', 'Count')
# output_count_path = paste0('~/Projects/twin_quad/analysis/metric_filter/mutation_count_none_20210722.tsv')
# write_delim(none_counts, output_count_path, delim='\t')
# dev.off()
# 
# 
# 
# ######################
# # 2021.07.02 
# MIN_COUNT=7
# grid.newpage()
# all_filtered_variants %>% mutate(ALT_TOTAL = ALT1_buccal + ALT2_buccal) %>% 
#   filter(ALT_TOTAL >= MIN_COUNT) %>% filter(str_detect(family, 'ST[1|2|3]')) %>% select(POS, REF_ALT, tricontext, DENOVO) %>% separate(tricontext, into=c('tri', 'sub'), sep='_') %>% filter(tri != 'indel') %>% group_by(DENOVO, sub, tri) %>% summarise(count=n()) %>% 
#   filter(DENOVO=='denovo') %>% right_join(df_96) %>% plot_mutation_signature
