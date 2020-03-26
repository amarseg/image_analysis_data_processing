load_CPA_data <- function(per_object_data, per_image_data, plate_n)
{
  require('tidyverse')
  object_data <- read_csv(per_object_data)
  well_data <- 
    read_csv(per_image_data) %>%
    dplyr::select(ImageNumber, Image_Metadata_QCFlag,Image_FileName_DNA) %>%
    separate(Image_FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
    transform(Well_n = str_sub(Well_n, start = 5, end = 7))
  
 strain_data <- read_csv('library_strains.csv') %>%
   dplyr::select('Ver5.0 position','Systematic ID') %>%
   separate('Ver5.0 position', into = c('Version','Plate','Well')) %>%
   filter(Plate == plate_n)
  
  all_data <- left_join(object_data, well_data, by = c('ImageNumber' = 'ImageNumber'))
  all_data_strain <- left_join(all_data, strain_data, by = c('Well_n' = 'Well'))
  return(all_data_strain)
}

load_filtered_data <- function(csv_path, plate_n)
{
  require(tidyverse)
  
  strain_data <- read_csv('library_strains.csv') %>%
    dplyr::select(`Ver5.0 position`,`Systematic ID`) %>%
    separate(`Ver5.0 position`, into = c('Version','Plate','Well')) %>%
    mutate(Plate = str_extract(Plate, pattern = '[:digit:]{2}')) %>%
    mutate(Plate = paste0('Plate',as.numeric(Plate)), Well = as.numeric(Well)) %>%
    filter(Plate == plate_n)
  
  data <- 
    read_csv(csv_path) %>% 
    mutate(File_name = basename(Metadata_FileLocation)) %>%
    separate(File_name, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
    filter(Metadata_QCFlag != 1) %>%
    left_join(strain_data, by = c( 'Metadata_WellNumber' = 'Well'))
  
  return(data)
}

plot_cell_proportion <- function(cell_count_data, bad_well_thr = 0.2)
{
  proportions <- cell_count_data %>%
    dplyr::select(Count_Nuclei, Count_Nuclei_AR_filtered, Count_Nuclei_AR_Solidity_Filtered, FileName_DNA, Metadata_Plate_Name) %>%
    separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
    mutate(Proportion_AR_filter = Count_Nuclei_AR_filtered/(Count_Nuclei+Count_Nuclei_AR_filtered), 
           Proportion_solidity_filter = Count_Nuclei_AR_Solidity_Filtered/(Count_Nuclei_AR_filtered+ Count_Nuclei_AR_Solidity_Filtered))
  
  test <- gather(cell_number, key = feature, value = cell_count, -plate_name, -Well) %>%
    filter(feature == 'Proportion_AR_filter' | feature == 'Proportion_solidity_filter')
  
  ggplot(test, aes(x = Well, group = feature, fill = feature, y = cell_count)) +
    geom_bar( stat = 'identity') +
    facet_wrap(~plate_name)
  
  ggplot(test, aes(group = feature, fill = feature, x = cell_count)) +
    geom_histogram() +
    facet_wrap(~plate_name)
  
  bad_wells <- filter(proportions, Proportion_solidity_filter < 0.2) %>%
    dplyr::select(plate_name, Well, Picture) %>%
    return()
  
}

draw_control_mean_dist <- function(plates_data, n, size, column)
{
  means <- rerun(n, sample_frac(plates_data[,column], size = size)) %>%
    sapply(function(x){mean(x[[1]])}) %>%
    as_data_frame() %>%
    return()
  
}

load_gene_lists <- function()
{
  require(tidyverse)
  
  up_rna <- read_tsv('up_rna.tsv') %>%
    add_column( type = 'Up RNA')
  
  down_rna <- read_tsv('down_rna.tsv') %>%
    add_column( type = 'Down RNA')
  
  up_prot <- read_tsv('up_prot.tsv') %>%
    add_column( type = 'Up prot')
  
  down_prot <- read_tsv('down_prot.tsv') %>%
    add_column( type = 'Down prot')
  
  omics_lists <- bind_rows(up_rna, down_rna, up_prot, down_prot) %>%
    return()
  
}

load_omics_data <- function()
{
  require(tidyverse)
  require(DESeq2)
  
  #RNa procesing
  rna <- read_csv('rna_seq.csv') 
  size_factors <- rna %>%
    dplyr::select(-ID) %>%
    estimateSizeFactorsForMatrix()
  
  rna[,-1] <- sweep(rna[,-1],MARGIN=2,FUN="/",STATS=size_factors) %>%
    as.data.frame()
  
  final_rna <- rna %>%
    gather(key = sample_name, value  = read_number, -ID) %>%
    separate(sample_name, into = c('not needed','time_point','replicate'), sep = '_') %>%
    dplyr::select(-`not needed`) %>%
    mutate(time_point = str_extract(time_point, pattern = '[:digit:]{1,2}')) %>%
    mutate(replicate = str_extract(replicate, pattern = '[:digit:]{1}')) %>%
    add_column(molecule = 'RNA')
  
  
  #Protein procesing
  prot <- read_delim('SQ_Results_PROTEIN.tsv', delim = '\t') %>%
    dplyr::select(proteinName, b018p004AM_T0_01:b018p004AM_T11_G3_02) %>%
    filter(str_detect(proteinName, pattern = 'SP')) %>%
    gather(key = sample_name, value = read_number, -proteinName) %>%
    separate(sample_name, into = c('no','time_point','replicate','technical_replicate'), sep = '_') %>%
    separate(proteinName, into = c('ID','gene_name','chr','desc'), sep = '\\|') %>%
    mutate(time_point = str_extract(time_point, pattern = '[:digit:]{1,2}')) %>%
    mutate(replicate = str_extract(replicate, pattern = '[:digit:]{1}')) %>%
    dplyr::select(-(gene_name:desc), -no) %>%
    add_column(molecule = 'Protein')
  
  t <- prot[which(is.na(prot$technical_replicate)),]
  t$technical_replicate <- t$replicate  
  t$replicate <-1
  
  prot[which(is.na(prot$technical_replicate)),] <- t
  
  mean_prot <- prot %>%
    group_by(ID, time_point, replicate, molecule) %>%
    summarise_at('read_number', mean)
  
  omics_dataset <- bind_rows(final_rna, mean_prot) %>%
    return()
  
  write_csv(omics_dataset, 'tidy_omics.csv')
}

fypo_database_loading <- function()
{
  require(tidyverse)
  require(ontologyIndex)
  
  fypo_def <- get_OBO('fypo-simple.obo.txt', propagate_relationships = "is_a", extract_tags = "minimal")
  
  fypo_df <- fypo_def$name %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename( 'FYPO ID' = 'rowname', 'Definition' = '.')
  
  fypo_db <- read_tsv('phenotype_annotations.pombase.phaf') %>%
    filter(str_detect(`Allele name`, pattern = 'delta')) %>%
    inner_join(fypo_df, by = 'FYPO ID') %>%
    dplyr::select(`Gene systematic ID`,`FYPO ID`,Definition)
  
  return(fypo_db)
}

load_go <- function(){
  require(AnnotationDbi)
  require(GO.db)
  
  go_data <- read_delim('ftp://ftp.pombase.org/pombe/annotations/Gene_ontology/gene_association.pombase.gz',
                        skip = 44, delim = '\t', col_names = F)
  term2gene <- go_data[,c(5,2)]
  
  term2name <- AnnotationDbi::select(GO.db, keys = keys(GO.db), columns = c('GOID','TERM'))
  
  go_database <- list(term2gene = term2gene, term2name = term2name)
  return(go_database)
}

load_angeli_for_fgsea <- function(){
  require(tidyverse)
  
  angeli_data <- read_delim('AnGeLiDatabase.txt', delim = '\t')
}