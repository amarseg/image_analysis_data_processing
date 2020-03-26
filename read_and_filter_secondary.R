library(tidyverse)


#######Load_data############
input_path <- 'Z:/Pers_Amalia/Screening/live_imaging_second_round_8_bit/outlines/'
csv_paths <- list.files(path = input_path,
                        pattern = 'yeastObjects',
                        full.names = TRUE,
                        recursive = TRUE)

data <- csv_paths %>%
  map(read_csv) %>%
  reduce(bind_rows)

data %>%
  select_at(vars(contains('AreaShape'))) %>%
  gather() %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~key, scales = 'free') +
  geom_vline(xintercept = 0.9, colour = 'red')

sys2name <- read_delim('../annotation/sysID2product.tsv',
                       delim = '\t', skip = 1, 
                       col_names = c('Systematic_ID','Name','Synonym')) %>%
  select(Systematic_ID, Name) %>%
  mutate(Name = case_when(is.na(Name) ~ Systematic_ID,
                          TRUE ~ Name))

strain_data <- read_csv('../annotation/strain_code_ready_for_use.csv') %>%
  mutate(Plate_name = str_extract(Plate_name, 'row_[0-9]{1}')) %>%
  rename(Systematic_ID = `Systematic ID`) 

filtered_data <- data %>% 
  inner_join(strain_data, by = c('Metadata_Well' = 'Real_well', 'Metadata_Row' = 'Plate_name')) %>%
  mutate(Aspect_ratio = AreaShape_MinorAxisLength/AreaShape_MajorAxisLength,
         Metadata_Time = as.numeric(Metadata_Time)) %>%
  filter(AreaShape_Area < 4500 & AreaShape_Solidity > 0.55 & AreaShape_MinorAxisLength < 40 
         & AreaShape_MajorAxisLength < 500 & Aspect_ratio < 0.5 & AreaShape_Compactness < 10 &
           AreaShape_Perimeter < 700 & AreaShape_MaximumRadius < 15 & Metadata_Time < 30) %>%
  write_csv('../processed_data/no_filter_live_imaging_data.csv')
