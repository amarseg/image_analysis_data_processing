library('tidyverse')
source('functions_CPA.R')

rep1 <- c('Z:/Pers_Amalia/Screening/Analysed_data/aug_sept_2018/')
rep2 <- c('Z:/Pers_Amalia/Screening/Analysed_data/replicate_2/', 
          'Z:/Pers_Amalia/Screening/Analysed_data/replicate_2_batch_2/',
          'Z:/Pers_Amalia/Screening/Analysed_data/batch_3/',
          'Z:/Pers_Amalia/Screening/Analysed_data/batch_4/')

plates <- c(1:36)
plates <- paste0('Plate',plates)
plates <- paste(plates, collapse = '|')

rep1_folders <- list.files(rep1, full.names = T, recursive = T, pattern = 'Nuclei_AR_Solidity_Filtered.csv')
rep1_extract <- grep(rep1_folders, pattern = plates, value = T)

rep2_folders <- list.files(rep2, full.names = T, recursive = T, pattern = 'Nuclei_AR_Solidity_Filtered.csv')[-2]
rep2_folders <- grep(rep2_folders, pattern = 'Plate', value = T)

cell_number <- read_csv('../Analysed_data/aug_sept_2018/cell_countsImage.csv') %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  mutate(Proportion_AR_filter = Count_Nuclei_AR_filtered/(Count_Nuclei+Count_Nuclei_AR_filtered), 
         Proportion_solidity_filter = Count_Nuclei_AR_Solidity_Filtered/(Count_Nuclei_AR_filtered+ Count_Nuclei_AR_Solidity_Filtered)) 

cell_count_files <- list.files(rep2, full.names = T, recursive = T, pattern = 'cell_countsImage.csv') %>%
  map(read_csv) %>%
  reduce(bind_rows)

cell_number_2 <- cell_count_files %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  mutate(Proportion_AR_filter = Count_Nuclei_AR_filtered/(Count_Nuclei+Count_Nuclei_AR_filtered), 
         Proportion_solidity_filter = Count_Nuclei_AR_Solidity_Filtered/(Count_Nuclei_AR_filtered+ Count_Nuclei_AR_Solidity_Filtered))


cell_number_df <- cell_number %>%
  select(Well,ImageQuality_Correlation_DNA_20:ImageQuality_TotalIntensity_DNA) %>%
  gather(key = 'Parameter', value = value, -Well)

ggplot(cell_number_df, aes(x = value)) +
  geom_histogram()+
  facet_wrap(~Parameter, scales = 'free')

cell_number2_df <- cell_number_2 %>%
  select(Well,ImageQuality_Correlation_DNA_20:ImageQuality_TotalIntensity_DNA) %>%
  gather(key = 'Parameter', value = value, -Well)

ggplot(cell_number2_df, aes(x = value)) +
  geom_histogram()+
  facet_wrap(~Parameter, scales = 'free')

###read dataaaa###
areas1 <- list()
i = 1
for(file in rep1_extract)
{
  plate_n <- str_extract(file, pattern = 'Plate[0-9]{1,2}')
  data <- load_filtered_data(file, plate_n = plate_n) %>%
    left_join(cell_number, by = c('Well' = 'Well','Metadata_Plate_Name','ImageNumber') ) %>%
    filter(ImageQuality_Correlation_DNA_20 < 0.5 & ImageQuality_PowerLogLogSlope_DNA < -2 & AreaShape_Solidity > 0.9) 
  
  areas1[[i]] <- data %>%
    add_count(`Systematic ID`) %>%
    filter(n >=50) %>%
    select(AreaShape_Area,`Systematic ID`,Metadata_Plate_Name, Well)
  i=i+1
}

areas2 <- list()
z_scores <- list()
i = 1
for(file in rep2_folders)
{
  plate_n <- str_extract(file, pattern = 'Plate[0-9]{1,2}')
  data <- load_filtered_data(file, plate_n = plate_n) %>%
    left_join(cell_number_2, by = c('Well' = 'Well','Metadata_Plate_Name','ImageNumber') ) %>%
    filter(ImageQuality_Correlation_DNA_20 < 0.5 & ImageQuality_PowerLogLogSlope_DNA < -2 & AreaShape_Solidity > 0.9) 
  
  areas2[[i]] <- data %>%
    add_count(`Systematic ID`) %>%
    filter(n >=50) %>%
    select(AreaShape_Area,`Systematic ID`,Metadata_Plate_Name, Well)
  
  
  i=i+1
}

areas1_tib <- bind_rows(areas1) %>%
  group_by(`Systematic ID`) %>%
  summarise(mean_area = mean(AreaShape_Area, trim = 0.1))

bind_rows(areas2) %>%
  write_csv('../processed_data/areas_rep1.csv')


areas2_tib <- bind_rows(areas2) %>%
  group_by(`Systematic ID`) %>%
  summarise(mean_area = mean(AreaShape_Area, trim = 0.1)) 

bind_rows(areas2) %>%
  write_csv('../processed_data/areas_rep2.csv')

all_areas <- inner_join(areas1_tib, areas2_tib, by = 'Systematic ID')


load('lists.011018.rda')
hits_summary <- stack(out)


z_score_tib <- bind_rows(areas1) %>%
  group_by(`Systematic ID`) %>%
  summarise(mean_area = mean(AreaShape_Area, trim = 0.1), sd_area = sd(AreaShape_Area)) %>%
  mutate(cv = sd_area/mean_area) %>%
  left_join(hits_summary, by = c('Systematic ID' = 'values')) %>%
  write_csv('../processed_data/statistics_rep1_good.csv')


z_score_tib2 <- bind_rows(areas2) %>%
  group_by(`Systematic ID`) %>%
  summarise(mean_area = mean(AreaShape_Area, trim = 0.1), sd_area = sd(AreaShape_Area)) %>%
  mutate(cv = sd_area/mean_area) %>%
  left_join(hits_summary, by = c('Systematic ID' = 'values')) %>%
  write_csv('../processed_data//statistics_rep2.csv')





