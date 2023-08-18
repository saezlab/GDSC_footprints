library(readr)
library(dplyr)
library(decoupleR)
library(reshape2)
library(readxl)

Cell_line_RMA_proc_basalExp <- as.data.frame(
  read_delim("data/Cell_line_RMA_proc_basalExp.txt", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE))[,-2]


GDSC_RMA <- Cell_line_RMA_proc_basalExp %>% group_by(GENE_SYMBOLS) %>% summarise_each(mean)
GDSC_RMA <- as.data.frame(GDSC_RMA)
GDSC_RMA <- GDSC_RMA[!is.na(GDSC_RMA$GENE_SYMBOLS),]

row.names(GDSC_RMA) <- GDSC_RMA$GENE_SYMBOLS

GDSC_RMA <- GDSC_RMA[,-1]

### progeny

progeny_model <- decoupleR::get_progeny(top = 100)

progeny_activities <- run_ulm(GDSC_RMA, progeny_model)

progeny_activities_wf <- reshape2::dcast(progeny_activities[,c(2,3,4)], source~condition)

write_csv(progeny_activities_wf, file = "results/GDSC_progeny_activities.csv")
write_csv(progeny_activities_wf, file = "results/GDSC_progeny_activities.csv")

### collectrI

collectrI <- decoupleR::get_collectri()

TF_activities <- run_ulm(GDSC_RMA, collectrI)

TF_activities_df <- reshape2::dcast(TF_activities[,c(2,3,4)], source~condition)

write_csv(collectrI, file = "support/collectrI.csv")
write_csv(TF_activities_df, file = "results/GDSC_TF_activities.csv")

#####

GDSC1_fitted_dose_response <- as.data.frame(
  read_excel("data/GDSC1_fitted_dose_response_24Jul22.xlsx"))

cell_name_mapping <- unique(GDSC1_fitted_dose_response[,c(4,5)])
cell_name_mapping$COSMIC_ID <- paste("DATA.",cell_name_mapping$COSMIC_ID, sep = "")

cell_name_mapping <- cell_name_mapping[cell_name_mapping$COSMIC_ID %in% names(progeny_activities_wf),]

write_csv(cell_name_mapping, file = "results/cell_name_mapping.csv")
