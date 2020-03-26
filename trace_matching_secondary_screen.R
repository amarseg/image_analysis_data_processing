rm(list = ls())
library(tidyverse)


filtered_data <- read_csv('../processed_data/no_filter_live_imaging_data.csv') %>%
  add_column(cell_n = 0) 
#filter(Metadata_Time <= 15)
####Create lineages######
center_thr <- 10

final_data = list()

tp_all <- filtered_data 

for(j in 0:29)
{
  
  tp0 <- tp_all %>%
    #group_by(Metadata_Row, Metadata_Position, Systematic_ID) %>% 
    mutate(Metadata_Time = as.numeric(Metadata_Time)) %>%
    filter(Metadata_Time == j & cell_n == 0) %>%
    mutate(cell_n = ObjectNumber)
  
  
  
  
  tp_all <- tp_all %>%
    filter(Metadata_Time != j & cell_n == 0)
  
  # print(dim(tp_all))
  
  for(i in 1:nrow(tp0))
  {
    tp_cell <- tp0[i,]
    
    
    if (nrow(tp_all[which(tp_all$Metadata_Row == tp_cell$Metadata_Row &
                          tp_all$Metadata_Position == tp_cell$Metadata_Position &
                          tp_all$Systematic_ID == tp_cell$Systematic_ID& 
                          tp_all$Location_Center_X > tp_cell$Location_Center_X - center_thr & 
                          tp_all$Location_Center_X < tp_cell$Location_Center_X + center_thr &
                          tp_all$Location_Center_Y > tp_cell$Location_Center_Y - center_thr &
                          tp_all$Location_Center_Y < tp_cell$Location_Center_Y + center_thr &
                          tp_all$Metadata_Well == tp_cell$Metadata_Well),]) >= 1)
    {
      
      
      n = nrow(tp_all[which(tp_all$Metadata_Row == tp_cell$Metadata_Row &
                              tp_all$Metadata_Position == tp_cell$Metadata_Position &
                              tp_all$Systematic_ID == tp_cell$Systematic_ID& 
                              tp_all$Location_Center_X > tp_cell$Location_Center_X - center_thr & 
                              tp_all$Location_Center_X < tp_cell$Location_Center_X + center_thr &
                              tp_all$Location_Center_Y > tp_cell$Location_Center_Y - center_thr &
                              tp_all$Location_Center_Y < tp_cell$Location_Center_Y + center_thr &
                              tp_all$Metadata_Well == tp_cell$Metadata_Well),])
      #print(n)
      tp_all[which(tp_all$Metadata_Row == tp_cell$Metadata_Row &
                     tp_all$Metadata_Position == tp_cell$Metadata_Position &
                     tp_all$Systematic_ID == tp_cell$Systematic_ID& 
                     tp_all$Location_Center_X > tp_cell$Location_Center_X - center_thr & 
                     tp_all$Location_Center_X < tp_cell$Location_Center_X + center_thr &
                     tp_all$Location_Center_Y > tp_cell$Location_Center_Y - center_thr &
                     tp_all$Location_Center_Y < tp_cell$Location_Center_Y + center_thr &
                     tp_all$Metadata_Well == tp_cell$Metadata_Well),]$cell_n <- rep(paste0(tp_cell$ObjectNumber,
                                                                                           '-', j,
                                                                                           '-', tp_cell$ImageNumber,
                                                                                           '-',tp_cell$Metadata_Position), n)
      
      tp0[i,]$cell_n <- tp_cell$ObjectNumber
    }
    
    
  }
  
  final_data[[j+1]] <- tp_all %>%
    filter(cell_n != 0)
}

bop <- do.call(rbind.data.frame, final_data)


############Saving data#######################


track_tp_all <- bop %>%
  write_csv('../processed_data/unfiltered_lineage_cell_data.csv')


