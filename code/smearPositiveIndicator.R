##### ABOUT THIS SCRIPT #######################################################
### This script reads in the clean dataset from data folder and then 
### replicates the methods in Horton et. al to determine change.  
### Function takes one stratum argument of: "sex", "rurality", "hiv" or "age"
smrPosIndicator <- function(stratum = "sex", 
                            saveCheckFiles = FALSE){
    ##############################################################################|
    ##### LOAD IN NECESSARY PACKAGES ##############################################
    library(dplyr)
    library(here)
    library(magrittr)
    ##############################################################################|
    ##### LOAD CLEAN DATA SET AND FILTER TO SURVEYS REPORTING STRATUM #############
    cleanStratumDF0 <- readRDS(here("data/fullDataClean.rds")) %>%
        filter(!!sym(paste0("report.",stratum))=="Yes")   
    
    ##############################################################################|
    ##### ADD WHO REGION AND STUDY ID TO DATA  ####################################
    regionWHO <- read.csv("data/who-regions.csv")[,c(1,4)] %>% 
        rename(study.country = Entity,
               WHO.region = World.regions.according.to.WHO)
    
    figID <- read.csv(here("data/titlesForFigures.csv"))
    colnames(figID)[1] <- "covidence.id"
    
    cleanStratumDF0  %<>% left_join(figID[,c("covidence.id", "figure.id")], 
                                    by = "covidence.id") %>% 
        mutate(figure.id.yr = paste(figure.id, study.years))
    
    ### Rename some countries to match the WHO CSV.
    ### For labeling these will remain as extracted. 
    
    cleanStratumDF0[which(cleanStratumDF0$study.country == "Viet Nam"), "study.country"] <- "Vietnam"
    cleanStratumDF0[which(cleanStratumDF0$study.country == "The Gambia"), "study.country"] <- "Gambia"
    cleanStratumDF0[which(cleanStratumDF0$study.country == "United Republic of Tanzania"), "study.country"] <- "Tanzania"
    cleanStratumDF0[which(cleanStratumDF0$study.country == "Lao PDR"), "study.country"] <- "Laos"
    cleanStratumDF0[which(cleanStratumDF0$study.country == "Democratic People's Republic of Korea"), "study.country"] <- "North Korea"
    
    ### Join the WHO data 
    cleanStratumDF0 %<>% left_join(regionWHO)
    
    ##############################################################################|
    ##### CREATE A NEW VARIABLE TO INDICATE WHICH ESTIMATE TO USE             #####
    
    ##### Each stratum will have a different result suffix so this will need to 
    ##### be defined based on the function argument. 
    
    suffix <- switch (stratum,
                      sex = {"male"}, 
                      rurality = {"rural"}, 
                      hiv = {"hiv.positive"},
                      age.grp = {"age.grp.2"},
    )
    
    ##### Filter according to the decision tree algorithm 
    ##### This can likely be simplified but wanted to be extra sure of algorithm at 
    ##### first. 
    newSmrStratumDF <- cleanStratumDF0 %>% mutate(!!paste0(stratum,".smr.analysis.indicator") := case_when(
        ### Smear positive adjusted prevalence with CI reported by survey
            !is.na(get(paste0("adj.prev100k.ci.smear.positive.tb.", suffix))) & 
            ! grepl(paste(c("multiplicati", "inflat"),collapse="|"), 
                    other.study.method.comments, ignore.case = TRUE) ~ "adj.prev100k.ci.smear.positive.tb", 
        ### Smear positive crude prevalence reported with CI by survey
            is.na(get(paste0("adj.prev100k.ci.smear.positive.tb.", suffix))) & 
            ! is.na(get(paste0("prev100k.ci.smear.positive.tb.", suffix))) & 
            !grepl(paste(c("multiplicati", "inflat"),collapse="|"), 
                   other.study.method.comments, ignore.case = TRUE) ~ "prev100k.ci.smear.positive.tb",
        ### Smear positive TB counts and participants reported by survey
            is.na(get(paste0("adj.prev100k.ci.smear.positive.tb.", suffix))) & 
            is.na(get(paste0("prev100k.ci.smear.positive.tb.", suffix))) & 
            !is.na(get(paste0("n.smear.positive.tb.", suffix))) &     
            !is.na(get(paste0("n.participants.", suffix)))  ~ "n.smear.positive.tb", 
        ### Smear positive TB counts and participants reported by survey
            is.na(get(paste0("adj.prev100k.ci.smear.positive.tb.", suffix))) & 
            is.na(get(paste0("prev100k.ci.smear.positive.tb.", suffix))) & 
            is.na(get(paste0("n.smear.positive.tb.", suffix))) &   
            !is.na(get(paste0("adj.prev100k.smear.positive.tb.", suffix))) & 
            !is.na(get(paste0("n.participants.", suffix))) ~ "adj.prev100k.smear.positive.tb",
        ### Smear positive TB counts and participants reported by survey
            is.na(get(paste0("adj.prev100k.ci.smear.positive.tb.", suffix))) & 
            is.na(get(paste0("prev100k.ci.smear.positive.tb.", suffix))) & 
            is.na(get(paste0("n.smear.positive.tb.", suffix))) &   
            is.na(get(paste0("adj.prev100k.smear.positive.tb.", suffix))) & 
            !is.na(get(paste0("prev100k.smear.positive.tb.", suffix))) & 
            !is.na(get(paste0("n.participants.", suffix))) ~ "prev100k.smear.positive.tb",
        .default = "none")) 
    
    print("distribution of samples for calculation:")
    print(table(newSmrStratumDF[,paste0(stratum,".smr.analysis.indicator")]))
    
    ##############################################################################|
    ##### CHECKS TO CONFIRM THE FINAL DATA ########################################
    ### Output the studies which we have no data for bacteriological positive
    ### analysis. 
    excluded <- newSmrStratumDF %>% filter(!!sym(paste0(stratum, ".smr.analysis.indicator")) == "none") %>%
        dplyr::select(title.covidence,
                      contains(suffix),
                      paste0(stratum, ".smr.analysis.indicator"),
                      smear.positive.tb.definition)              
    
    ##### SAVE CHECK FILES IF OPTION IS TRUE ######################################
    if (saveCheckFiles == TRUE){
        write.csv(x = excluded, file = here(paste0("data/checks/excludedSmearSurveys_", suffix, ".csv")))                
    }
    
    return(newSmrStratumDF)
    
}
