# 2022.05.21 cjyoon, monochorionic mixing, dichorionic unmixing figure
# 2025.02.16 cjyoon mle to estimate twin mixing in the blood and plotting.Draws both VAF changes in linear graph and stacked bar plot
# using only the whole genome sequencing ref alt, vaf counts
# 2025.08.10 betabinom with rho WGS=0.0002568

library(tidyverse)
library(readxl)
library(grid)
library(plotly)
library(VGAM)

rho <- 0.0002568 # rho wgs


setwd('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/')

# MLE functions for estimating mixing ratio
likelihood_mixing <- function(alt, depth, buccal_vaf1, buccal_vaf2, x1) {
  alt_allele_fraction <- buccal_vaf1 * x1 + buccal_vaf2 * (1 - x1)
  
  # Convert p, rho to alpha/beta
  alpha <- alt_allele_fraction * (1 - rho) / rho
  beta  <- (1 - alt_allele_fraction) * (1 - rho) / rho
  
  prob <- dbetabinom.ab(alt, size = depth, shape1 = alpha, shape2 = beta) # betabinom
  
  return(log10(prob))
}


mixing_nll_blood1<-function(df, par){
  if(par>0){
    return(-sum(likelihood_mixing(df$twin1_blood_wgs_alt, df$twin1_blood_wgs_depth, df$twin1_buccal_wgs_vaf, df$twin2_buccal_wgs_vaf, par)))
  }else{
    return(10000000000000000)
  }
}

mixing_nll_blood2<-function(df, par){
  if(par>0){
    return(-sum(likelihood_mixing(df$twin2_blood_wgs_alt, df$twin2_blood_wgs_depth, df$twin1_buccal_wgs_vaf, df$twin2_buccal_wgs_vaf, par)))
  }else{
    return(10000000000000000)
  }
}




family_ids = str_extract(list.files('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/', pattern = '^ST[0-9]+.reannotated_manual.xlsx$'), '^ST[0-9]+')
family_ids = c("ST104", "ST105", "ST200", 'ST201', 'ST202', 'ST203', 'ST204', 'ST207', 'ST301', 'ST302', 'ST304')

mixing_ratio = data.frame()

for(family_id in family_ids){
  df = read_xlsx(paste0(family_id, '.reannotated_manual.xlsx'))
  df2 = read_delim(paste0('./FINAL_2D_TABLE/', family_id, '.final.vaf.table.tsv'), delim='\t')
  if('twin1_blood_wgs_vaf' %in% colnames(df)){ # do this analysis only if blood whole genome sequencing was conducted
    print(family_id)
    # filter out manually removed variants
    df2 = df2  %>% mutate(variant=paste0(CHROM, ':', POS,  ' ', REF , '>', ALT)) # other true variants are empty on manual column so NA but should be retained
    
    # buccal1 exclusive variants
    twin1_exclusive = df2 %>% filter(twin1_buccal_tgs_vaf >0 & twin2_buccal_tgs_vaf == 0 & twin2_buccal_wgs_vaf == 0) 
    twin1_exclusive$only = 'twin1'
    
    # buccal2 exclusive variants
    twin2_exclusive = df2 %>% filter(twin2_buccal_tgs_vaf >0 & twin1_buccal_tgs_vaf == 0 & twin1_buccal_wgs_vaf == 0)  
    twin2_exclusive$only = 'twin2'
    exclusive_df = bind_rows(twin1_exclusive, twin2_exclusive)
    
    
      
      x1_blood1 = mle2(mixing_nll_blood1, start=list(par=0), lower=list(par=0), upper=list(par=1), data=list(df=exclusive_df), method='Brent')
      x1_blood2 = mle2(mixing_nll_blood2, start=list(par=0), lower=list(par=0), upper=list(par=1), data=list(df=exclusive_df), method='Brent')
      
      x1_blood1 = coef(x1_blood1)[[1]]
      x1_blood2 = coef(x1_blood2)[[1]]
      
      mixing_ratio_a = data.frame(family_id = c(family_id, family_id), blood = c("Twin1", "Twin2"), x1 = c(x1_blood1, x1_blood2))
      mixing_ratio = rbind(mixing_ratio, mixing_ratio_a)

      if(dim(twin1_exclusive)[1]>0){
        twin1_exclusive %>% select(variant, twin1_buccal_wgs_vaf, twin2_buccal_wgs_vaf, twin1_blood_wgs_vaf, twin2_blood_wgs_vaf) %>% 
          pivot_longer(cols=c('twin1_buccal_wgs_vaf', 'twin2_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin2_blood_wgs_vaf'), names_to = 'source', values_to='vaf') %>% 
          mutate(source = factor(source, levels=c('twin1_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin2_blood_wgs_vaf', 'twin2_buccal_wgs_vaf'))) %>% 
          # ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(aes(col=exclusive), size=4) + 
          ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(col='black', fill='#000080', size=4, alpha=0.8, pch=21) + 
          xlab('') + ylab('VAF') + theme_classic(base_size=24) + ylim(0,1) + 
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
          scale_x_discrete(labels = c('twin1_buccal_wgs_vaf' = 'Twin1 Buccal',
                                      'twin1_blood_wgs_vaf' = 'Twin1 Blood',
                                      'twin2_blood_wgs_vaf' = 'Twin2 Blood',
                                      'twin2_buccal_wgs_vaf' = 'Twin2 Buccal')) +
          ggtitle('')
        # ggsave(paste0('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/monochorionic_mixing_2d/', family_id, '_mixing2.pdf'), width=5, height=5)
      }
      
      if(dim(twin2_exclusive)[1]>0){
        twin2_exclusive %>% select(variant, twin1_buccal_wgs_vaf, twin2_buccal_wgs_vaf, twin1_blood_wgs_vaf, twin2_blood_wgs_vaf) %>% 
          pivot_longer(cols=c('twin1_buccal_wgs_vaf', 'twin2_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin2_blood_wgs_vaf'), names_to = 'source', values_to='vaf') %>% 
          mutate(source = factor(source, levels=c('twin1_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin2_blood_wgs_vaf', 'twin2_buccal_wgs_vaf'))) %>% 
          # ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(aes(col=exclusive), size=4) + 
          ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(col='black', fill='#000080', size=4, alpha=0.8, pch=21) + 
          xlab('') + ylab('VAF') + theme_classic(base_size=24) + ylim(0,1) + 
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
          scale_x_discrete(labels = c('twin1_buccal_wgs_vaf' = 'Twin1 Buccal',
                                      'twin1_blood_wgs_vaf' = 'Twin1 Blood',
                                      'twin2_blood_wgs_vaf' = 'Twin2 Blood',
                                      'twin2_buccal_wgs_vaf' = 'Twin2 Buccal')) +
          ggtitle(family_id)
        # ggsave(paste0('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/monochorionic_mixing_2d/', family_id, '_mixing1.pdf'), width=5, height=5)
      }
      

    }
}
mixing_ratio %>% write_delim('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/mixing_ratio_betabinom.tsv', delim = '\t')
mixing_ratio = mixing_ratio %>% mutate(y1=1-x1)

buccal_ratio = data.frame(family_id = c('Buccal', 'Buccal'), blood=c('Twin1', 'Twin2'), x1=c(1, 0), y1=c(0, 1))
mixing_ratio = rbind(buccal_ratio, mixing_ratio)
# Convert data from wide to long format

df_long <- mixing_ratio %>%
  pivot_longer(cols = c(x1, y1), names_to = "component", values_to = "value")
df_long <- df_long %>%
  mutate(component = factor(component, levels = c("x1", "y1")))  # reverse so x1 ends up at bottom
# df_long <- df_long %>%
#   mutate(component = factor(component, levels = c("y1", "x1")))  # reverse so x1 ends up at bottom

ggplot(df_long, aes(x = blood, y = value, fill = component)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~family_id, nrow = 1) +
  scale_fill_manual(values = c("x1" = "steelblue", "y1" = "orange")) +
  labs(x = "", y = "Blood Contribution", fill = "Component") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 12, face = "bold")
  )


ggsave('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/monochorionic_mixing_2d/mixing_ratios_betabinom.pdf', width=17, height=2.5)



# Fitting x-y for showing mixing ratio
y1_wide <- mixing_ratio %>%
  select(family_id, blood, y1) %>%
  pivot_wider(names_from = blood, values_from = y1, names_prefix = "y1_") %>%
  filter(family_id != "Buccal")  # remove Buccal row if not needed

# View result
print(y1_wide)

y1_wide <- y1_wide %>%
  mutate(
    family_type = case_when(
      grepl("^ST1", family_id) ~ "MCMA",
      grepl("^ST2", family_id) ~ "MCDA",
      grepl("^ST3", family_id) ~ "DCDA",
      TRUE ~ NA_character_
    )
  )


# Example scatter plot
mixing_ratio_plot = ggplot(y1_wide, aes(x = y1_Twin1, y = y1_Twin2)) +
  # All points
  geom_point(size = 12, alpha = 0.9, col = '#888888', aes(shape=family_type)) +
  
  # Fit only on families not starting with ST3, with confidence interval
  geom_smooth(
    data = subset(y1_wide, !grepl("^ST3", family_id)),
    method = "lm", formula = y ~ x, col = "blue", se = TRUE, fill = "lightblue", alpha = 0.3
  ) +
  
  xlab('Twin2 Conbtribution\n in Twin1 Blood') +
  ylab('Twin2 Contribution\n in Twin2 Blood') +
  
  scale_x_continuous(
    limits = c(0, 1.1),
    expand = c(0, 0),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = c('0.0', '', '0.5', '', '1.0')
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    expand = c(0, 0),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = c('', '', '0.5', '', '1.0')
  ) +
  scale_shape_manual(values = c("MCMA" = 15,  # square
                                "MCDA" = 17,  # triangle
                                "DCDA" = 16))+ # circle

  theme_classic() +
  ggtitle('') +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30, hjust = 0.5),
    axis.title.y = element_text(size = 30),
    axis.title.x = element_text(size = 30, vjust = -1),
    axis.text.x = element_text(vjust = 0),
    axis.text.y = element_text(hjust = -0.2),
    axis.line = element_line(size = 1.5), 
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks = element_line(size = 1)
  )
mixing_ratio_plot

fit_nonST3 <- lm(y1_Twin2 ~ y1_Twin1, data = subset(y1_wide, !grepl("^ST3", family_id)))
summary(fit_nonST3)
# lm(formula = y1_Twin2 ~ y1_Twin1, data = subset(y1_wide, !grepl("^ST3", 
#                                                                 family_id)))
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.075366 -0.047313 -0.005776  0.013254  0.150450 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   0.1540     0.0815   1.890  0.10764   
# y1_Twin1      0.6926     0.1398   4.956  0.00256 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07782 on 6 degrees of freedom
# Multiple R-squared:  0.8037,	Adjusted R-squared:  0.771 


pgt = ggplot_gtable(ggplot_build(mixing_ratio_plot))
pgt$layout$clip[pgt$layout$name=="panel"] <- "off"
pdf(paste0('~/Documents/twin_2025/FIGURES/equal_mixing.pdf'), width=7, height=7)
pushViewport(viewport(x=0.5, y=0.5, width=0.93, height=0.93)) # to fit 1.0 axis label
grid.draw(pgt)
dev.off()






  # for triplet
# for triplet
family_id == "ST501"
df2 = read_delim(paste0('./FINAL_2D_TABLE/', family_id, '.final.vaf.table.tsv'), delim='\t')


  twin1_exclusive_vs_twin2 = df2 %>% filter(twin1_buccal_tgs_vaf >0 & ((twin2_buccal_tgs_vaf == 0 & twin2_buccal_wgs_vaf == 0))) 
  twin1_exclusive_vs_twin2$only = 'twin1vs2'
  
  twin1_exclusive_vs_twin3 = df2 %>% filter(twin1_buccal_tgs_vaf >0 & ((twin3_buccal_tgs_vaf == 0 & twin3_buccal_wgs_vaf == 0))) 
  twin1_exclusive_vs_twin3$only = 'twin1vs3'
  
  
  # buccal2 exclusive variants
  twin2_exclusive_vs_twin1 = df2 %>% filter(twin2_buccal_tgs_vaf >0 & ((twin1_buccal_tgs_vaf == 0 & twin1_buccal_wgs_vaf == 0)))
  twin2_exclusive_vs_twin1$only = 'twin2vs1'
  
  twin2_exclusive_vs_twin3 = df2 %>% filter(twin2_buccal_tgs_vaf >0 & ((twin3_buccal_tgs_vaf == 0 & twin3_buccal_wgs_vaf == 0)))
  twin2_exclusive_vs_twin3$only = 'twin2vs3'
  
  # buccal3 exclusive variants
  twin3_exclusive_vs_twin1 = df2 %>% filter(twin3_buccal_tgs_vaf >0 & ((twin1_buccal_tgs_vaf == 0 & twin1_buccal_wgs_vaf == 0)))
  twin3_exclusive_vs_twin1$only = 'twin3vs1'
  
  twin3_exclusive_vs_twin2 = df2 %>% filter(twin3_buccal_tgs_vaf >0 & ((twin2_buccal_tgs_vaf == 0 & twin2_buccal_wgs_vaf == 0)))
  twin3_exclusive_vs_twin2$only = 'twin3vs2'
  
  
  
  # exclusive_df = bind_rows(twin1_exclusive_vs_twin2, twin1_exclusive_vs_twin3, twin2_exclusive_vs_twin1, twin2_exclusive_vs_twin3, twin3_exclusive_vs_twin1, twin3_exclusive_vs_twin2)
  
  
  if(dim(twin1_exclusive_vs_twin2)[1]>0){
    twin1_exclusive_vs_twin2 %>% select(variant, twin1_buccal_wgs_vaf, twin2_buccal_wgs_vaf, twin1_blood_wgs_vaf, twin2_blood_wgs_vaf) %>% 
      pivot_longer(cols=c('twin1_buccal_wgs_vaf', 'twin2_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin2_blood_wgs_vaf'), names_to = 'source', values_to='vaf') %>% 
      mutate(source = factor(source, levels=c('twin1_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin2_blood_wgs_vaf', 'twin2_buccal_wgs_vaf'))) %>% 
      # ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(aes(col=exclusive), size=4) + 
      ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(col='black', fill='#000080', size=4, alpha=0.8, pch=21) + 
      xlab('') + ylab('VAF') + theme_classic(base_size=24) + ylim(0,1) + 
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
      scale_x_discrete(labels = c('twin1_buccal_wgs_vaf' = 'Twin1 Buccal',
                                  'twin1_blood_wgs_vaf' = 'Twin1 Blood',
                                  'twin2_blood_wgs_vaf' = 'Twin2 Blood',
                                  'twin2_buccal_wgs_vaf' = 'Twin2 Buccal')) +
      ggtitle('')
    ggsave(paste0('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/monochorionic_mixing_2d/', family_id, '_mixing_buccal1vs2.pdf'), width=5, height=5)
  }
  

if(dim(twin1_exclusive_vs_twin3)[1]>0){
  twin1_exclusive_vs_twin3 %>% select(variant, twin1_buccal_wgs_vaf, twin3_buccal_wgs_vaf, twin1_blood_wgs_vaf, twin3_blood_wgs_vaf) %>% 
    pivot_longer(cols=c('twin1_buccal_wgs_vaf', 'twin3_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin3_blood_wgs_vaf'), names_to = 'source', values_to='vaf') %>% 
    mutate(source = factor(source, levels=c('twin1_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin3_blood_wgs_vaf', 'twin3_buccal_wgs_vaf'))) %>% 
    # ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(aes(col=exclusive), size=4) + 
    ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(col='black', fill='#000080', size=4, alpha=0.8, pch=21) + 
    xlab('') + ylab('VAF') + theme_classic(base_size=24) + ylim(0,1) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
    scale_x_discrete(labels = c('twin1_buccal_wgs_vaf' = 'Twin1 Buccal',
                                'twin1_blood_wgs_vaf' = 'Twin1 Blood',
                                'twin3_blood_wgs_vaf' = 'Twin3 Blood',
                                'twin3_buccal_wgs_vaf' = 'Twin3 Buccal')) +
    ggtitle('')
  ggsave(paste0('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/monochorionic_mixing_2d/', family_id, '_mixing_buccal1vs3.pdf'), width=5, height=5)
}

if(dim(twin2_exclusive_vs_twin1)[1]>0){
  twin2_exclusive_vs_twin1 %>% select(variant, twin1_buccal_wgs_vaf, twin2_buccal_wgs_vaf, twin1_blood_wgs_vaf, twin2_blood_wgs_vaf) %>% 
    pivot_longer(cols=c('twin1_buccal_wgs_vaf', 'twin2_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin2_blood_wgs_vaf'), names_to = 'source', values_to='vaf') %>% 
    mutate(source = factor(source, levels=c('twin1_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin2_blood_wgs_vaf', 'twin2_buccal_wgs_vaf'))) %>% 
    # ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(aes(col=exclusive), size=4) + 
    ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(col='black', fill='#000080', size=4, alpha=0.8, pch=21) + 
    xlab('') + ylab('VAF') + theme_classic(base_size=24) + ylim(0,1) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
    scale_x_discrete(labels = c('twin1_buccal_wgs_vaf' = 'Twin1 Buccal',
                                'twin1_blood_wgs_vaf' = 'Twin1 Blood',
                                'twin2_blood_wgs_vaf' = 'Twin2 Blood',
                                'twin2_buccal_wgs_vaf' = 'Twin2 Buccal')) +
    ggtitle('')
  ggsave(paste0('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/monochorionic_mixing_2d/', family_id, '_mixing_buccal2vs1.pdf'), width=5, height=5)
}

if(dim(twin2_exclusive_vs_twin3)[1]>0){
  twin2_exclusive_vs_twin3 %>% select(variant, twin2_buccal_wgs_vaf, twin3_buccal_wgs_vaf, twin2_blood_wgs_vaf, twin3_blood_wgs_vaf) %>% 
    pivot_longer(cols=c('twin2_buccal_wgs_vaf', 'twin3_buccal_wgs_vaf', 'twin2_blood_wgs_vaf', 'twin3_blood_wgs_vaf'), names_to = 'source', values_to='vaf') %>% 
    mutate(source = factor(source, levels=c('twin2_buccal_wgs_vaf', 'twin2_blood_wgs_vaf', 'twin3_blood_wgs_vaf', 'twin3_buccal_wgs_vaf'))) %>% 
    # ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(aes(col=exclusive), size=4) + 
    ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(col='black', fill='#000080', size=4, alpha=0.8, pch=21) + 
    xlab('') + ylab('VAF') + theme_classic(base_size=24) + ylim(0,1) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
    scale_x_discrete(labels = c('twin2_buccal_wgs_vaf' = 'Twin2 Buccal',
                                'twin2_blood_wgs_vaf' = 'Twin2 Blood',
                                'twin3_blood_wgs_vaf' = 'Twin3 Blood',
                                'twin3_buccal_wgs_vaf' = 'Twin3 Buccal')) +
    ggtitle('')
  ggsave(paste0('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/monochorionic_mixing_2d/', family_id, '_mixing_buccal2vs3.pdf'), width=5, height=5)
}


if(dim(twin3_exclusive_vs_twin1)[1]>0){
  twin3_exclusive_vs_twin1 %>% select(variant, twin1_buccal_wgs_vaf, twin3_buccal_wgs_vaf, twin1_blood_wgs_vaf, twin3_blood_wgs_vaf) %>% 
    pivot_longer(cols=c('twin1_buccal_wgs_vaf', 'twin3_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin3_blood_wgs_vaf'), names_to = 'source', values_to='vaf') %>% 
    mutate(source = factor(source, levels=c('twin1_buccal_wgs_vaf', 'twin1_blood_wgs_vaf', 'twin3_blood_wgs_vaf', 'twin3_buccal_wgs_vaf'))) %>% 
    # ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(aes(col=exclusive), size=4) + 
    ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(col='black', fill='#000080', size=4, alpha=0.8, pch=21) + 
    xlab('') + ylab('VAF') + theme_classic(base_size=24) + ylim(0,1) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
    scale_x_discrete(labels = c('twin1_buccal_wgs_vaf' = 'Twin1 Buccal',
                                'twin1_blood_wgs_vaf' = 'Twin1 Blood',
                                'twin3_blood_wgs_vaf' = 'Twin3 Blood',
                                'twin3_buccal_wgs_vaf' = 'Twin3 Buccal')) +
    ggtitle('')
  ggsave(paste0('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/monochorionic_mixing_2d/', family_id, '_mixing_buccal3vs1.pdf'), width=5, height=5)
}

if(dim(twin3_exclusive_vs_twin2)[1]>0){
  twin3_exclusive_vs_twin2 %>% select(variant, twin3_buccal_wgs_vaf, twin2_buccal_wgs_vaf, twin3_blood_wgs_vaf, twin2_blood_wgs_vaf) %>% 
    pivot_longer(cols=c('twin3_buccal_wgs_vaf', 'twin2_buccal_wgs_vaf', 'twin3_blood_wgs_vaf', 'twin2_blood_wgs_vaf'), names_to = 'source', values_to='vaf') %>% 
    mutate(source = factor(source, levels=c('twin2_buccal_wgs_vaf', 'twin2_blood_wgs_vaf', 'twin3_blood_wgs_vaf', 'twin3_buccal_wgs_vaf'))) %>% 
    # ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(aes(col=exclusive), size=4) + 
    ggplot(aes(x=source, y=vaf, group=variant)) + geom_line() + geom_point(col='black', fill='#000080', size=4, alpha=0.8, pch=21) + 
    xlab('') + ylab('VAF') + theme_classic(base_size=24) + ylim(0,1) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
    scale_x_discrete(labels = c('twin2_buccal_wgs_vaf' = 'Twin2 Buccal',
                                'twin2_blood_wgs_vaf' = 'Twin2 Blood',
                                'twin3_blood_wgs_vaf' = 'Twin3 Blood',
                                'twin3_buccal_wgs_vaf' = 'Twin3 Buccal')) +
    ggtitle('')
  ggsave(paste0('~/Google Drive/My Drive/mztwin/VARIANT_TABLE/monochorionic_mixing_2d/', family_id, '_mixing_buccal3vs2.pdf'), width=5, height=5)
}





