###############################################################################
##### This function calculates the logOddsRatio standard error from adjusted
##### and crude confidence intervals for both bacteriological confirmed and 
##### smear positive TB prevalence estimates. Then, taking a variance component
##### inspired approach, it calculates the approximate contribution of the 
##### cluster variance as reported by adjusted confidence intervals for each
##### survey. The script then calculates the mean and median for each of these 
##### and reports the smear positive and bacteriological TB estimates in a 
##### vector that can be incorporated into the meta analysis. 
###############################################################################
############ THIS SCRIPT WILL NOT WORK FOR THE AGE GROUP ANALYSIS ############|
### Technically, it does work for two age groups but will need to be amended
### for use across the non-binary analytical strata. 
###############################################################################
##### Function takes one stratum argument of: "sex", "rurality", "hiv", "age"

varianceAdjust <- function(stratum = "sex"){
##### LOAD IN NECESSARY PACKAGES ##############################################
library(dplyr)
library(tidyverse)
library(here)
library(ggplot2)
library(reshape2)


##### Each stratum will have a different result suffixes so this will need to 
##### be defined based on the function argument. 
##### Suffix needs to match the numerator of the analysis
##### This will need to be updated for age analysis
suffix <- switch (stratum,
                  sex = {"male"}, 
                  rurality = {"urban"}, 
                  hiv = {"hiv.positive"},
                  age.grp = {"age.grp.2"},
)

##### Suffix needs to match the denominator of the analysis
##### This will need to be updated for age analysis
suffix2 <- switch (stratum,
                  sex = {"female"}, 
                  rurality = {"rural"}, 
                  hiv = {"hiv.negative"},
                  age.grp = {"age.grp.5"},
)

###############################################################################
######  BACTERIOLOGICALLY CONFIRMED TB ########################################
###############################################################################

##### LOAD CLEAN DATA SET AND FILTER TO SURVEYS REPORTING SEX #################
source(here("code/bacterialPositiveIndicator.R"))
bactVarDF <- bactPostIndicator(stratum) %>% 
##### Filter to those surveys which have adjusted and crude estimates 
            filter((!is.na(!!sym(paste0("adj.prev100k.ci.bacteriological.tb.", suffix))) &
                   !is.na(!!sym(paste0("prev100k.ci.bacteriological.tb.", suffix)))) == TRUE) %>%
    ### keep only relevant parameters
    select(covidence.id, 
           figure.id.yr, 
           study.geography, 
           title.extracted,
           study.country,
           WHO.region,
           study.start.year, 
           study.end.year, 
           ends_with(suffix), ends_with(suffix2)) %>%
    ### Scale all prevalence estimates back to decimals 
    mutate_at(vars(matches('prev100k.bacteriological.tb.|prev100k.ci.([a-z]{3})er.bacteriological.tb.')),  ~ (. / 1e5)) %>%
    ### Create an ID variable 
    mutate("id" = row_number(),   
           ### Calculate the log odds standard error for the numerator and denominators 
           LogOddsStandardErrorAdjSuffix = (log(get(paste0("adj.prev100k.ci.upper.bacteriological.tb.", suffix))/
                                                    (1-get(paste0("adj.prev100k.ci.upper.bacteriological.tb.", suffix)))) -
                                                log(get(paste0("adj.prev100k.ci.lower.bacteriological.tb.", suffix))/
                                                        (1-get(paste0("adj.prev100k.ci.lower.bacteriological.tb.", suffix))))) /3.92,
           LogOddsStandardErrorAdjSuffix2 =  (log(get(paste0("adj.prev100k.ci.upper.bacteriological.tb.", suffix2))/
                                                      (1-get(paste0("adj.prev100k.ci.upper.bacteriological.tb.", suffix2)))) -
                                                  log(get(paste0("adj.prev100k.ci.lower.bacteriological.tb.", suffix2))/
                                                          (1-get(paste0("adj.prev100k.ci.lower.bacteriological.tb.", suffix2))))) /3.92,
           LogOddsStandardErrorCrudeSuffix = (log(get(paste0("prev100k.ci.upper.bacteriological.tb.", suffix))/
                                                      (1-get(paste0("prev100k.ci.upper.bacteriological.tb.", suffix)))) -
                                                  log(get(paste0("prev100k.ci.lower.bacteriological.tb.", suffix))/
                                                          (1-get(paste0("prev100k.ci.lower.bacteriological.tb.", suffix))))) /3.92,
           LogOddsStandardErrorCrudeSuffix2 = (log(get(paste0("prev100k.ci.upper.bacteriological.tb.", suffix2))/
                                                       (1-get(paste0("prev100k.ci.upper.bacteriological.tb.", suffix2)))) -
                                                   log(get(paste0("prev100k.ci.lower.bacteriological.tb.", suffix2))/
                                                           (1-get(paste0("prev100k.ci.lower.bacteriological.tb.", suffix2))))) /3.92,
           #### Calculate the standard error of the log odds ratio for adjusted and crude
           logOR_Adjse = sqrt(LogOddsStandardErrorAdjSuffix^2 + LogOddsStandardErrorAdjSuffix2^2),
           logOR_Crdse = sqrt(LogOddsStandardErrorCrudeSuffix^2 + LogOddsStandardErrorCrudeSuffix2^2),
           #### Calculate the component of the variance contributed by 
           #### clustering as reported by adjusted estimates
           logOR_ClstVar = (logOR_Adjse^2) - (logOR_Crdse^2)) 

##### Calculate the mean and median of the cluster variance ###################
bactMean <- mean(bactVarDF$logOR_ClstVar, na.rm=TRUE)
bactMedian <- median(bactVarDF$logOR_ClstVar, na.rm=TRUE)

###############################################################################
###### SMEAR POSITIVE TB ######################################################
###############################################################################
source(here("code/smearPositiveIndicator.R"))
cleanStratSmrDF <- smrPosIndicator(stratum) %>% 
               filter(!!sym(paste0(stratum, ".smr.analysis.indicator")) != "none")
if (stratum == "sex"){
##### Survey from Gambia has a zero lower bound; replace with 1 to allow inclusion in the analysis. 
    cleanStratSmrDF[which(cleanStratSmrDF$covidence.id==26920), 
            "adj.prev100k.ci.lower.smear.positive.tb.female"] <- 1/1e5
}

##### Filter to those surveys which have adjusted and crude estimates 
smrVarDF <- cleanStratSmrDF %>% 
    filter((!is.na(!!sym(paste0("adj.prev100k.ci.smear.positive.tb.", suffix))) &
            !is.na(!!sym(paste0("prev100k.ci.smear.positive.tb.", suffix)))) == TRUE) %>%
    ### keep only relevant parameters
    select(covidence.id, 
           figure.id.yr, 
           study.geography, 
           title.extracted,
           study.country,
           WHO.region,
           study.start.year, 
           study.end.year, 
           ends_with(suffix), ends_with(suffix2)) %>%
    ### Scale all prevalence estimates back to decimals 
    mutate_at(vars(matches('prev100k.smear.positive.tb.|prev100k.ci.([a-z]{3})er.smear.positive.tb.')),  ~ (. / 1e5)) %>%
    ### Create an ID variable 
    mutate("id" = row_number(),   
           ### Calculate the log odds standard error for the numerator and denominators 
           LogOddsStandardErrorAdjSuffix = (log(get(paste0("adj.prev100k.ci.upper.smear.positive.tb.", suffix))/
                                                  (1-get(paste0("adj.prev100k.ci.upper.smear.positive.tb.", suffix)))) -
                                              log(get(paste0("adj.prev100k.ci.lower.smear.positive.tb.", suffix))/
                                                      (1-get(paste0("adj.prev100k.ci.lower.smear.positive.tb.", suffix))))) /3.92,
           LogOddsStandardErrorAdjSuffix2 =  (log(get(paste0("adj.prev100k.ci.upper.smear.positive.tb.", suffix2))/
                                                     (1-get(paste0("adj.prev100k.ci.upper.smear.positive.tb.", suffix2)))) -
                                                 log(get(paste0("adj.prev100k.ci.lower.smear.positive.tb.", suffix2))/
                                                         (1-get(paste0("adj.prev100k.ci.lower.smear.positive.tb.", suffix2))))) /3.92,
           LogOddsStandardErrorCrudeSuffix = (log(get(paste0("prev100k.ci.upper.smear.positive.tb.", suffix))/
                                                    (1-get(paste0("prev100k.ci.upper.smear.positive.tb.", suffix)))) -
                                                log(get(paste0("prev100k.ci.lower.smear.positive.tb.", suffix))/
                                                        (1-get(paste0("prev100k.ci.lower.smear.positive.tb.", suffix))))) /3.92,
           LogOddsStandardErrorCrudeSuffix2 = (log(get(paste0("prev100k.ci.upper.smear.positive.tb.", suffix2))/
                                                      (1-get(paste0("prev100k.ci.upper.smear.positive.tb.", suffix2)))) -
                                                  log(get(paste0("prev100k.ci.lower.smear.positive.tb.", suffix2))/
                                                          (1-get(paste0("prev100k.ci.lower.smear.positive.tb.", suffix2))))) /3.92,
           #### Calculate the standard error of the log odds ratio for adjusted and crude
           logOR_Adjse = sqrt(LogOddsStandardErrorAdjSuffix^2 + LogOddsStandardErrorAdjSuffix2^2),
           logOR_Crdse = sqrt(LogOddsStandardErrorCrudeSuffix^2 + LogOddsStandardErrorCrudeSuffix2^2),
           #### Calculate the component of the variance contributed by 
           #### clustering as reported by adjusted estimates. 
           logOR_ClstVar = (logOR_Adjse^2) - (logOR_Crdse^2)) 

##### Calculate the mean and median of the cluster variance ###################
smrMean <- mean(smrVarDF$logOR_ClstVar, na.rm=TRUE)
smrMedian <- median(smrVarDF$logOR_ClstVar, na.rm=TRUE)

adjustment <- c(bactMean, bactMedian, smrMean, smrMedian)
names(adjustment) <- c("Bact Variance Mean", "Bact Variance Median",
                       "Smear Variance Mean", "Smear Variance Median")
return(adjustment)
}