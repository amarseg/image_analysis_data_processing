rm(list = ls())
library(tidyverse)
library(ggpubr)

gene_names <- read_tsv('analysis_live_imaging/annotation/sysID2product.tsv', skip = 1,
                       col_names = c('Systematic_ID','Name','Other_names','Description')) %>%
  mutate(Name = case_when(is.na(Name) ~ Systematic_ID,
                          TRUE ~ Name)) %>%
  add_row(Systematic_ID = 'wt', Name = 'wt', Other_names = 'wt',Description = 'wt')

no_trace_data <-  read_csv('analysis_live_imaging/output_data/no_filter_live_imaging_data.csv') %>%
  left_join(gene_names, by = 'Systematic_ID') %>%
  filter(Systematic_ID != 'unknown') %>%
  filter(AreaShape_Area < 5000) %>%
  mutate(Systematic_ID = case_when(Systematic_ID == 'wt' ~ paste0('wt', Metadata_Well, Metadata_Row),
                                   TRUE ~ Systematic_ID),
         Name = case_when(Systematic_ID == 'wt' ~ paste0('wt', Metadata_Well, Metadata_Row),
                          TRUE ~ Name))

primary_screen <- read_csv('analysis_live_imaging/annotation/hits_i_guess.csv') %>%
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

  

ggplot(all_screens, aes(x = Metadata_Time,  y = median_area,
                        group = Systematic_ID, colour = rep_1_area)) +
  geom_point() +
  geom_line() +
  theme_test()

ggplot(all_screens, aes(x = Metadata_Time,  y = median_area,
                        group = Systematic_ID, colour = rep_2_area)) +
  geom_point() +
  geom_line() +
  theme_test()

all_screens_max <- all_screens %>%
  group_by(Systematic_ID, rep_1_area, rep_2_area, size) %>%
  summarise(max_area -= max(median_area)) %>%
  gather(key = 'variable', value = median_area, -Systematic_ID, -max_area, - size)

ggplot(all_screens_max, aes(x = max_area,
                            y = median_area, colour = size)) +
  geom_point() +
  facet_wrap(~variable) +
  xlab('Maximal area in live imaging screening') +
  ylab('Area in primary screening') +
  theme_test()
  