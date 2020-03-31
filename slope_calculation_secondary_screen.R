rm(list = ls())
library(tidyverse)
library(ggrepel)
library(cowplot)
library(nlme)
library(broom)
library(ggpubr)


gene_names <- read_tsv('../Annotation/sysID2product.tsv', skip = 1,
                       col_names = c('Systematic_ID','Name','Other_names','Description')) %>%
  mutate(Name = case_when(is.na(Name) ~ Systematic_ID,
                          TRUE ~ Name)) %>%
  add_row(Systematic_ID = 'wt', Name = 'wt', Other_names = 'wt',Description = 'wt')

lineage_data <-  read_csv('../processed_data/unfiltered_lineage_cell_data.csv') %>%
  filter(cell_n != 0) %>%
  mutate(cell_n = paste0(Systematic_ID,'_', cell_n)) %>%
  group_by(Systematic_ID, Metadata_Well, cell_n) %>%
  mutate(number_of_timepoints = n()) %>%
  ungroup() %>%
  mutate(Systematic_ID = case_when(str_detect(Systematic_ID, 'wt') ~ paste0('wt',Metadata_Well, Metadata_Row),
                                   TRUE ~ Systematic_ID)) %>%
  unite(new_id,  sep = ';', Systematic_ID, cell_n, Metadata_Well, Metadata_Row, remove= FALSE) %>%
  filter(number_of_timepoints > 5)


############Find time_point where size is maximal and cut there, to then apply lm#############
max_tps <- lineage_data %>%
  group_by(new_id) %>%
  mutate(max_size = max(AreaShape_Area)) %>%
  filter(max_size == AreaShape_Area) %>%
  select(new_id, max_size, Metadata_Time, Systematic_ID) %>%
  rename(max_time = Metadata_Time) 


ggplot(max_tps,aes(x = reorder(Name, max_time, FUN = 'median'), y = max_time)) +
  geom_boxplot() +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 




all_data_unfil <- lineage_data %>%
  inner_join(max_tps, by = c('new_id','Systematic_ID')) 

#############Do lm models in first 10 timepoints##################
# 
early_tps <- lineage_data %>%
  filter(Metadata_Time <= 10)


lm_models_early <- lmList(log(AreaShape_Area) ~ Metadata_Time|new_id,
                          data = early_tps)


model_coeff_early <- lapply(lm_models_early, tidy) %>%
  bind_rows(.id = 'id') %>%
  separate(id, sep = ';', into = c('Systematic_ID','cell_n','well','row', remove = FALSE))

annot_early<- model_coeff_early %>%
  filter(term == 'Metadata_Time') %>%
  mutate(growing_not_growing = case_when(estimate >= 0.03 ~ 'Growing',
                                         estimate < 0.03 ~ 'Not Growing')) %>%
  select(Systematic_ID, cell_n, growing_not_growing) %>%
  write_csv('analysis_live_imaging/output_data/growing_or_not.csv')

##############Growing not growing according to fold change###################

max_min_fold_change <- lineage_data %>%
  group_by(Systematic_ID, new_id) %>%
  summarise(max_size = max(AreaShape_Area), min_size = min(AreaShape_Area)) %>%
  mutate(fold_change = max_size/min_size) %>%
  mutate(growing_not_growing = case_when(max_size > 900 ~ 'Growing',
                                         max_size <  900 ~ 'Not Growing'))


####################merge with all data#######################

all_data <- lineage_data 
#left_join(annot_early, by = c('cell_n','Systematic_ID')) %>%
#group_by(new_id) %>%
# filter((Metadata_Time <= max_time & growing_not_growing == 'Growing') | 
#          growing_not_growing == 'Not Growing')


############lm list for all############################


lm_models <- lmList(log(AreaShape_Area) ~ Metadata_Time|new_id, 
                    data = filter(lineage_data, Metadata_Time < 15))

model_coeff_lm <- lapply(lm_models, tidy) %>%
  bind_rows(.id = 'id') %>%
  separate(id, sep = ';', into = c('Systematic_ID','cell_n','well','row'), remove = FALSE) %>%
  left_join(gene_names, by = 'Systematic_ID') 

write_csv(model_coeff_lm, '../processed_data/linear_models_unfiltered.csv')


slopes <- model_coeff_lm %>%
  filter(term == 'Metadata_Time') 


slopes_for_output <- model_coeff_lm %>%
  filter(term == 'Metadata_Time') %>%
  mutate(p_adj = p.adjust(p.value),
         doubling_rate = exp(estimate) - 1) %>%
  filter(p_adj < 0.05) %>%
  group_by(Systematic_ID) %>%
  summarise(median_slope = median(estimate), sd_slope = sd(estimate, na.rm = T)) %>%
  arrange(desc(median_slope)) 

slopes_p_val_filtered <- model_coeff_lm %>%
  filter(term == 'Metadata_Time') %>%
  mutate(growing_not_growing = case_when(estimate >= 0.03 ~ 'Growing',
                                         estimate < 0.03 ~ 'Not Growing',
                                         str_detect(Systematic_ID, 'wt') ~ 'wt')) %>%
  mutate(p_adj = p.adjust(p.value),
         doubling_rate = exp(estimate) - 1) 


slopes_pval_annot <- slopes_p_val_filtered %>%
  filter(p_adj < 0.05) %>%
  mutate(Name = case_when(str_detect(Systematic_ID, 'wt') ~ Systematic_ID,
                          TRUE ~ Name)) %>%
  mutate(wild_type = case_when(str_detect(Systematic_ID, 'wt') ~ 'wt',
                               TRUE ~ 'not_wt')) 

max_sizes <- all_data %>%
  inner_join(slopes_pval_annot, by = c('Systematic_ID','cell_n')) %>%
  group_by(Systematic_ID, Name, cell_n, growing_not_growing, wild_type, row) %>%
  summarise(max_size = max(AreaShape_Area)) 



max_sizes %>%
  right_join(slopes_p_val_filtered, by = c('cell_n')) %>%
  write_csv('../processed_data/slopes_maximal_size.csv')

