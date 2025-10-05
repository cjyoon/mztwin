# script for finalizing p and q for mz twin pairs
# uses two dimensional clustering if true-identical
# uses one dimensional clustering if sub-identical with known clonal mutation (from read phasing)

library(tidyverse)
library(readxl)
library(grid)
library(plotly)
library(ggrepel)


setwd("~/Documents/twin_2025/FINAL_2D_TABLE/")
family_ids = str_extract(list.files('~/Documents/twin_2025/FINAL_2D_TABLE//', pattern = '^ST[1-5][0-9]+.final.vaf.table.tsv$'), '^ST[1-5][0-9]+')
twin_type <- function(family_id){
  if(str_detect(family_id, '^ST1')){
    return('MCMA')
  }else if(str_detect(family_id, '^ST2')){
    return('MCDA')
  }else if(str_detect(family_id, '^ST3')){
    return('DCDA')
  }else if(str_detect(family_id, '^ST5')){
    return('DCTA Triplet')
  }else if(str_detect(family_id, '^ST9')){
    return("DZ")
  }else{
    return(NA)
  }
}
twin_type_vec = Vectorize(twin_type)

anchor_pq = function(df, p_clonal = FALSE) {
  # find anchor pq positions
  df_eem = df %>%
    filter(dnm == "eem") %>%
    filter(!is.na(vaf1) & !is.na(vaf2)) %>%
    select(vaf1, vaf2) %>%
    rowid_to_column("id")
  # function to assign clusters
  assign_clusters = function(df, threshold = 0.1) {
    clusters = list()
    assigned = rep(FALSE, nrow(df))
    
    for (i in seq_len(nrow(df))) {
      if (assigned[i]) next
      current = df[i, ]
      cluster_indices = which(
        !assigned &
          abs(df$vaf1 - current$vaf1) < threshold &
          abs(df$vaf2 - current$vaf2) < threshold
      )
      clusters[[length(clusters) + 1]] <- df[cluster_indices, ]
      assigned[cluster_indices] <- TRUE
    }
    
    return(clusters)
  }
  
  # create clusters and centroids
  clusters = assign_clusters(df_eem, threshold = 0.1)
  centroids = lapply(clusters, function(cluster) colMeans(cluster[, c("vaf1", "vaf2")]))
  
  # Try to find best pair
  best_pair = NULL
  best_midpoint = NULL
  best_dist = Inf
  
  for (i in seq_along(centroids)) {
    for (j in seq_along(centroids)) {
      if (i >= j) next
      midpoint = (centroids[[i]] + centroids[[j]]) / 2
      dist = sqrt(sum((midpoint - c(0.25, 0.25))^2))
      if (dist < best_dist) {
        best_dist = dist
        best_pair = list(centroids[[i]], centroids[[j]])
        best_midpoint = midpoint
      }
    }
  }
  
  # use pair or fallback to one cluster
  use_fallback = is.null(best_pair) ||
    abs(best_midpoint[1] - 0.25) > 0.1 ||
    abs(best_midpoint[2] - 0.25) > 0.1
  
  if (use_fallback) {
    # Find the cluster whose centroid is closest to (0.25, 0.25)
    fallback_centroid = centroids[[which.min(sapply(centroids, function(c) sqrt(sum((c - c(0.25, 0.25))^2))))]]
    
    p = fallback_centroid[1]
    q = fallback_centroid[2]
    if (!(p >= 0.25 && p <= 0.5 &&
          q >= 0 && q <= 0.5 &&
          p >= q &&
          (p + q) >= 0.5)) {
      p = 0.5 -p
      q = 0.5 -q
    }
    
    result_df <- tibble::tibble(
      cluster = c("cluster1", 'cluster2'),
      vaf1 = c(p, 0.5-p),
      vaf2 = c(q, 0.5-q), 
      singleton= TRUE
    )
  } else {
    # Use best pair of clusters
    c1 = best_pair[[1]]
    c2 = best_pair[[2]]
    p = round(c1[1], 1)
    q = round(c1[2], 1)
    
    # Swap if cluster1 doesn't meet constraint
    if (!(p >= 0.25 && p <= 0.5 &&
          q >= 0 && q <= 0.5 &&
          p >= q &&
          (p + q) >= 0.5)) {
      tmp = c1
      c1 = c2
      c2 = tmp
      print("swap")
    }
    
    midpoint = (c1 + c2) / 2
    
    result_df = tibble::tibble(
      cluster = c("cluster1", "cluster2", "midpoint"),
      vaf1 = c(c1[1], c2[1], midpoint[1]),
      vaf2 = c(c1[2], c2[2], midpoint[2])
    )
  }
  return(result_df)
}

anchor_df = data.frame(row.names = c('cluster', 'singleton', 'vaf1', 'vaf2', 'family'))

for(family_id in family_ids){
  df <- read_delim(paste0('~/Documents/twin_2025/FINAL_2D_TABLE/', family_id, '.final.vaf.table.tsv'), delim = '\t')
  
  if(family_id != 'ST501'){
    aresult = anchor_pq(df)  
    aresult$family = family_id
    anchor_df = bind_rows(aresult, anchor_df)
    
  }else{
    # for ST501 need to create 3 different dfs for triplet 1vs2 2vs3 and 1vs3
    # also need to remove the EEM1 variants (chr6:166554073 G>A, chr2:27526275 T>C, and chr6:122,471,575 C>A) when calculating the anchor clusters for 1vs2 only but not for 2vs3 or 1vs3 to avoid having an anchor at 0.5, 0.5 which is not what we want in a comparison
    # for 2v3 and 1v3, since the anchor_pq function depend on having a column name as vaf1 and vaf2, need to change the column names in order for the function to process the data
    df1v2 = df %>% mutate(variant=paste0(CHROM, ':', POS, ' ', REF, '>', ALT)) %>% filter(!(variant %in% c('chr6:166554073 G>A', 'chr2:27526275 T>C', 'chr6:122471575 C>A')))
    aresult = anchor_pq(df1v2)  
    aresult$family = 'ST501_1v2'
    anchor_df = bind_rows(aresult, anchor_df)
    
    # change vaf2->vaf1, vaf3 -> vaf2
    df2v3 = df %>% mutate(variant=paste0(CHROM, ':', POS, ' ', REF, '>', ALT))
    df2v3%>% select(CHROM, POS, REF, ALT, vaf1, vaf2) %>% head
    
    names(df2v3) <- sub("^vaf", "vaf_temp", names(df2v3))
    names(df2v3)[names(df2v3) == "vaf_temp2"] <- "vaf1"  # vaf2 -> vaf1
    names(df2v3)[names(df2v3) == "vaf_temp3"] <- "vaf2"  # vaf3 -> vaf2
    df2v3%>% select(CHROM, POS, REF, ALT, vaf1, vaf2) %>% head
    
    
    aresult = anchor_pq(df2v3)  
    aresult$family = 'ST501_2v3'
    anchor_df = bind_rows(aresult, anchor_df)
    
    # change vaf3 -> vaf2
    df1v3 = df %>% mutate(variant=paste0(CHROM, ':', POS, ' ', REF, '>', ALT)) 
    names(df1v3) <- sub("^vaf", "vaf_temp", names(df1v3))
    names(df1v3)[names(df1v3) == "vaf_temp1"] <- "vaf1"  # vaf1 -> vaf1 keep it same
    names(df1v3)[names(df1v3) == "vaf_temp3"] <- "vaf2"  # vaf3 -> vaf2
    
    aresult = anchor_pq(df1v3)  
    aresult$family = 'ST501_1v3'
    anchor_df = bind_rows(aresult, anchor_df)
    
  }
  
  
  
  
  
  
}

write_delim(anchor_df, '~/Documents/twin_2025/TABLE/initial_anchor.tsv', delim='\t')

# anchor_df = read_delim('~/Documents/twin_2025/TABLE/initial_anchor.tsv', delim='\t')


summary_df_1 = anchor_df %>% filter(cluster=='cluster1') %>% mutate(family_type = twin_type_vec(family))
summary_df_2 = anchor_df %>% filter(cluster=='cluster2') %>% mutate(family_type = twin_type_vec(family))
summary_df_mid = anchor_df %>% filter(cluster=='midpoint') %>% mutate(family_type = twin_type_vec(family))
 
# read in a manually phased excel sheet to adjust p value to 0.5 for those that are phased to clonal mutation
phasing_clonal_p = read_xlsx('~/Documents/twin_2025/TABLE/p_clonal_phased.xlsx')

summary_df_1 = summary_df_1 %>% left_join(phasing_clonal_p, by='family') 
summary_df_2 = summary_df_2 %>% left_join(phasing_clonal_p, by='family') 

adjust_p_clonal <- function(vaf1, p_clonal){
  if(p_clonal=='T'){
    return(0.5)
  }else{
    return(vaf1)
  }
}
adjust_p_clonal_vec = Vectorize(adjust_p_clonal)

adjust_q_clonal <- function(vaf1, p_clonal){
  if(p_clonal=='T'){
    return(0)
  }else{
    return(vaf1)
  }
}
adjust_q_clonal_vec = Vectorize(adjust_q_clonal)



summary_df_1 = summary_df_1 %>% mutate(vaf1 = adjust_p_clonal_vec(vaf1, p_clonal))
summary_df_2 = summary_df_2 %>% mutate(vaf1 = adjust_q_clonal_vec(vaf1, p_clonal))

adjusted_df = rbind(summary_df_1, summary_df_2) %>% arrange(family)


write_delim(adjusted_df, '~/Documents/twin_2025/TABLE/p_clonal_adjusted_anchor.tsv', delim='\t')


# col_fun = c('DCDA'='#008080', 'MCDA'='#E69F00', 'MCMA'='#756FB3')
col_fun = c('DCDA'='black', 'MCDA'='black', 'MCMA'='black')
shape_fun = c('DCDA'=16, 'MCDA'=17, 'MCMA' = 15)

anchor_area = data.frame(x=c(0.5, 1, 1), y=c(0.5, 0, 1), group=c(1,1,1), col='1') 
anchor_area = anchor_area %>% mutate(x=0.5*x, y=0.5*y)
geom_polygon(data=anchor_area, aes(x=x, y=y), alpha=0.25, col='black', fill='black')

####
# for this plot, while we calculated anchor clusters for mz twins and also triplet, in the manuscript we only talked about twins upto Fig2 so will exclude any data on ST501 triplet
summary_df_1 = summary_df_1 %>% filter(!str_detect(family, 'ST501'))

p = ggplot() + 
  # Draw the background area where anchor points can be located
  geom_polygon(data=anchor_area, aes(x=x, y=y), alpha=0.1, col='black', fill='black') + 
  # Draw the points for each mz twin
  geom_point(data=summary_df_1, aes(x=vaf1, y=vaf2, shape=family_type, label=family), size=10, alpha=0.95) +
  xlim(0.22,0.55) + ylim(0, 0.55) + 
  # xlab(expression(scriptstyle(frac(1,2)) ~ "p estimate")) +  # Inline small 1/2 without parentheses
  # ylab(expression(scriptstyle(frac(1,2)) ~ "q estimate")) +  # Inline small 1/2 without parentheses
  
  xlab(expression(italic(frac(a,2)) ~ " estimate")) +  # Italicized a/2
  ylab(expression(italic(frac(b,2)) ~ " estimate")) +  # Italicized b/2
  # scale_color_manual(values=col_fun) +
  scale_shape_manual(values=shape_fun) + 
  scale_x_continuous(limit=c(0.25, 0.55),  breaks=c(0.25, 0.375, 0.5), labels=c('0.25', '', '0.5')) + 
  scale_y_continuous(limit=c(0, 0.6), breaks=c(0, 0.25, 0.5), labels=c('0', '', '0.5')) + 
  theme_minimal() + theme(legend.position = "none",axis.text = element_text(size = 80),
                          axis.title = element_text(size = 70, margin = margin(t = 20, r = 10, b = 30, l = 10))  # Push y-axis label downward
  ) + 
  # Add labels for each point
  # Repelling labels for each point
  geom_text_repel(
    data = summary_df_1,
    aes(x = vaf1, y = vaf2, label = family),
    size = 20,
    max.overlaps = Inf,   # allow more labels, but avoid overlaps
    force = 4,             # push labels away more strongly
    box.padding = 2,
    point.padding = 2,
    min.segment.length = 1,
    segment.size = 1
    
  ) 
p
ggsave('~/Documents/twin_2025/FIGURES/clonal_dynamics_summary_20250518_col_v2_label.pdf', width=20, height=20)

p = p +   xlab('0.5 p') +  ylab('0.5 q') 


# zoomed in graph as there are too many points around (0.5, 0)

ggplot() + 
  # Draw the background area where anchor points can be located
  # Draw the points for each mz twin
  geom_point(data=summary_df_1, aes(x=vaf1, y=vaf2, col=family_type, shape=family_type, label=family), size=6, alpha=0.95) +
  # xlab(expression(scriptstyle(frac(1,2)) ~ "p estimate")) +  # Inline small 1/2 without parentheses
  # ylab(expression(scriptstyle(frac(1,2)) ~ "q estimate")) +  # Inline small 1/2 without parentheses
  
  scale_color_manual(values=col_fun) +
  scale_shape_manual(values=shape_fun) + scale_x_continuous(limits = c(0.47, 0.53), breaks=c(0.47, 0.5, 0.53), labels=c('', '0.5', '')) +
  scale_y_continuous(limits = c(0, 0.05),
                     breaks = c(0, 0.05), labels = c('0', '0.05'),     sec.axis = dup_axis(name = NULL)) +  # duplicate axis on right)+
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),       # remove x and y labels
    axis.text = element_text(size = 20),# keep ticks if you want
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.y.left  = element_blank(),
    axis.text.y.right = element_text(size = 20),
  ) + 
  geom_text_repel(data = summary_df_1,
                  aes(x = vaf1, y = vaf2, label = family),
                  size = 6,
                  max.overlaps = Inf,
                  box.padding = 1,
                  point.padding = 0,
                  force = 2)

ggsave('~/Documents/twin_2025/FIGURES/clonal_dynamics_summary_20250518_col_v2_label_zoomed.pdf', width=3, height=3)



# zoomed in graph #2 as there are too many points around (0.5, 0)

ggplot() + 
  # Draw the background area where anchor points can be located
  # Draw the points for each mz twin
  geom_point(data=summary_df_1, aes(x=vaf1, y=vaf2, col=family_type, shape=family_type, label=family), size=6, alpha=0.95) +
  # xlab(expression(scriptstyle(frac(1,2)) ~ "p estimate")) +  # Inline small 1/2 without parentheses
  # ylab(expression(scriptstyle(frac(1,2)) ~ "q estimate")) +  # Inline small 1/2 without parentheses
  
  scale_color_manual(values=col_fun) +
  scale_shape_manual(values=shape_fun) + scale_x_continuous(limits = c(0.47, 0.53), breaks=c(0.47, 0.5, 0.53), labels=c('', '0.5', '')) +
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = c(0, 0.25), labels = c('0', '0.25'),     sec.axis = dup_axis(name = NULL)) +  # duplicate axis on right)+
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),       # remove x and y labels
    axis.text = element_text(size = 20),# keep ticks if you want
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.y.left  = element_blank(),
    axis.text.y.right = element_text(size = 20),
  ) + 
  geom_text_repel(data = summary_df_1,
                  aes(x = vaf1, y = vaf2, label = family),
                  size = 3,
                  max.overlaps = Inf,
                  box.padding = 1,
                  point.padding = 0,
                  force = 2)

ggsave('~/Documents/twin_2025/FIGURES/clonal_dynamics_summary_20250518_col_v2_label_zoomed2.pdf', width=3, height=3)


  
  
figure_path = "~/Documents/twin_2025/FIGURES/clonal_dynamics_summary_20250518_col2_ggplotly_20250518.html"
  ggp = ggplotly(p)
  ggplotly_figure_dir = '~/Documents/twin_2025/FIGURES/'
  htmlwidgets::saveWidget(ggp, figure_path, selfcontained = TRUE)
