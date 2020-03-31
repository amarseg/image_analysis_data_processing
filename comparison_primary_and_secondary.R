rm(list = ls())
library(tidyverse)
library(ggpubr)

load_gene_lists <- function()
{
  require(tidyverse)
  
  up_rna <- read_tsv('../../gene_lists//up_rna.tsv') %>%
    add_column( type = 'Up RNA')
  
  down_rna <- read_tsv('../../gene_lists//down_rna.tsv') %>%
    add_column( type = 'Down RNA')
  
  up_prot <- read_tsv('../../gene_lists/up_prot.tsv') %>%
    add_column( type = 'Up prot')
  
  down_prot <- read_tsv('../../gene_lists/down_prot.tsv') %>%
    add_column( type = 'Down prot')
  
  omics_lists <- bind_rows(up_rna, down_rna, up_prot, down_prot) %>%
    rename(Systematic_ID = `Systematic ID`,
           Product = `Product description`) %>%
    mutate(Name = case_when(is.na(Name) ~ Systematic_ID,
                            TRUE ~ Name)) %>%
    return()
  
}

gene_lists <- load_gene_lists()

gene_names <- read_tsv('../Annotation/sysID2product.tsv', skip = 1,
                       col_names = c('Systematic_ID','Name','Other_names','Description')) %>%
  mutate(Name = case_when(is.na(Name) ~ Systematic_ID,
                          TRUE ~ Name)) %>%
  add_row(Systematic_ID = 'wt', Name = 'wt', Other_names = 'wt',Description = 'wt')

no_trace_data <-  read_csv('../processed_data/no_filter_live_imaging_data.csv') %>%
  left_join(gene_names, by = 'Systematic_ID') %>%
  filter(Systematic_ID != 'unknown') %>%
  filter(AreaShape_Area < 5000) %>%
  mutate(Systematic_ID = case_when(Systematic_ID == 'wt' ~ paste0('wt', Metadata_Well, Metadata_Row),
                                   TRUE ~ Systematic_ID),
         Name = case_when(Systematic_ID == 'wt' ~ paste0('wt', Metadata_Well, Metadata_Row),
                          TRUE ~ Name))

primary_screen <- read_csv('../Annotation/hits_i_guess.csv') %>%
  rename(Systematic_ID = `Systematic ID`)

all_screens <- no_trace_data %>%
  inner_join(primary_screen, by = c('Systematic_ID'))  %>%
  group_by(Systematic_ID, Metadata_Time) %>%
  summarise(median_area_live_imaging = median(AreaShape_Area),
            rep_1_area = unique(mean_area_rep1),
            rep_2_area = unique(mean_area_rep2),
            sd_rep_1 = unique(sd_area_rep1),
            sd_rep_2 = unique(sd_area_rep2),
            cv_rep_1 = unique(cv_rep1),
            cv_rep_2 = unique(cv_rep2),
            mean_area_live_imaging = mean(AreaShape_Area),
            sd_live_imaging = sd(AreaShape_Area),
            cv_live_imaging = sd_live_imaging/mean_area_live_imaging) %>%
  write_csv('summary_data_hits.csv')

  

ggplot(all_screens, aes(x = Metadata_Time,  y = median_area_live_imaging,
                        group = Systematic_ID, colour = rep_1_area)) +
  geom_point() +
  geom_line() +
  theme_test()

ggplot(all_screens, aes(x = Metadata_Time,  y = median_area_live_imaging,
                        group = Systematic_ID, colour = rep_2_area)) +
  geom_point() +
  geom_line() +
  theme_test()

all_screens_max <- all_screens %>%
  group_by(Systematic_ID, rep_1_area, rep_2_area) %>%
  summarise(max_area = max(median_area_live_imaging)) %>%
  gather(key = 'variable', value = median_area, -Systematic_ID, -max_area)

ggplot(all_screens_max, aes(x = max_area,
                            y = median_area)) +
  geom_point() +
  facet_wrap(~variable) +
  xlab('Maximal area in live imaging screening') +
  ylab('Area in primary screening') +
  theme_test()

############Join omics and hits################
all_data <- all_screens %>%
  left_join(gene_lists,by = 'Systematic_ID') %>%
  write_csv('../processed_data/annotated_hits_with_omics.csv')

ggplot(filter(all_data, !is.na(type)), 
       aes(x = Metadata_Time, y = median_area_live_imaging,
                     colour = type, group = Systematic_ID)) +
  geom_ribbon(aes(ymax = median_area_live_imaging + sd_live_imaging, 
                  ymin = median_area_live_imaging -sd_live_imaging,
                  alpha = 0.4)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Name)

ggplot(all_data, 
       aes(x = Metadata_Time, y = median_area_live_imaging,
           colour = type, group = Systematic_ID)) +
  geom_line() +
  geom_point() +
  facet_wrap(~type)
