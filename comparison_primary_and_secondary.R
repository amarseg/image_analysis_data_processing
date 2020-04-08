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

hand_data <- read_csv('../Annotation/hand_annotation.csv')

model_data <- read_csv('../processed_data/slopes_maximal_size.csv')

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
######################Comparing models omics########
annot_models <- model_data %>%
  left_join(gene_lists, by = c('Systematic_ID')) %>%
  left_join(hand_data, by = c('Systematic_ID' = 'Systematic ID'))

ggplot(annot_models, aes(x = reorder(Systematic_ID, slope, 'median'), y = slope,
                         fill = type)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(annot_models, aes(x = reorder(Systematic_ID, max_size, 'median'), y = max_size,
                         fill = type)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(annot_models, aes(x = reorder(Systematic_ID, slope, 'median'), y = slope,
                         fill = Grow_not_grow)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(annot_models, aes(x = reorder(Name, max_size, 'median'), y = max_size,
                         fill = Grow_not_grow)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##############all data + hand made#############
hand_all_data <- all_data %>%
  inner_join(hand_data, by = c('Systematic_ID' = 'Systematic ID'))

ggplot(hand_all_data, aes(x = Metadata_Time, y = median_area_live_imaging,
                          group = Systematic_ID, colour = Grow_not_grow)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Grow_not_grow)
###########p value versus slope
ggplot(filter(annot_models, !is.na(Grow_not_grow)), aes(x = p_val, y = slope,
                          colour = Grow_not_grow)) +
  geom_point() +
  facet_wrap(~Grow_not_grow)

ggplot(filter(annot_models, !is.na(Grow_not_grow)), aes(x = Systematic_ID, y = slope,
                                                        colour = Grow_not_grow)) +
  geom_boxplot() 

ggplot(filter(annot_models, !is.na(Grow_not_grow)), aes(x = Systematic_ID, y = p_val,
                                                        colour = Grow_not_grow)) +
  geom_boxplot() 

#########
pval_cutoff <- annot_models %>%
  mutate(growth = case_when(slope < 0.01 ~ 'No',
                            TRUE ~ 'Yes')) %>%
  left_join(gene_names, by = c('Systematic_ID'))  %>%
  mutate(Name = case_when(is.na(Name) ~ Systematic_ID,
                          TRUE ~ Name))%>%
  write_csv('../processed_data/results_secondary_screen.csv')

ggplot(pval_cutoff, aes(x = reorder(Name, slope, 'median'), y = slope,
                         fill = Grow_not_grow)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~growth)

ggplot(pval_cutoff, aes(x = reorder(Name, max_size, 'median'), y = max_size,
                         fill = Grow_not_grow)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~growth)
##########propotion growth and number of cells#########
pval_cutoff %>%
  group_by(Systematic_ID, growth, Grow_not_grow, Name) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Name, y = n,
             colour= Grow_not_grow)) +
  geom_point() +
  facet_wrap(~growth) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


pval_cutoff %>%
  group_by(Systematic_ID,Grow_not_grow, Name) %>%
  summarise(n_growing = sum(growth == 'Yes'),
            n_not_growing = sum(growth == 'No')) %>%
  mutate(proportion = n_growing/(n_growing + n_not_growing)) %>%
  ggplot(aes(x = reorder(Name, proportion), y = proportion,
             colour= Grow_not_grow)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


#########colour according to omics#####
ggplot(pval_cutoff, aes(x = reorder(Name, slope, 'median'), y = slope,
                        fill = type)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~growth)

ggplot(pval_cutoff, aes(x = reorder(Name, max_size, 'median'), y = max_size,
                        fill = type)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~growth)

########Primary screen and omics############
annot_primary <- read_csv('../processed_data/output_rep1/summary_stats_rep1.csv') %>%
  rename('Systematic_ID' = `Systematic ID`) %>%
  left_join(gene_lists, by = c('Systematic_ID')) 

ggplot(annot_primary, aes(x = type, 
                          y = mean_area, 
                          colour = type)) +
  geom_boxplot() +
  geom_jitter()


annot_primary <- read_csv('../processed_data/output_rep2/statistics_rep2.csv') %>%
  rename('Systematic_ID' = `Systematic ID`) %>%
  left_join(gene_lists, by = c('Systematic_ID')) 

ggplot(annot_primary, aes(x = type, 
                          y = mean_area, 
                          colour = type)) +
  geom_boxplot() +
  geom_jitter()
