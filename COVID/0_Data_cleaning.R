###############################################################################
# Code originally written by Daniel J. Phillips and Elaine Shuo Feng for a
# previous analysis of the same data, and adapted for this analysis by DJP.
# Functions and packages

require(plyr)
require(stringr)
require(openxlsx)

# directories
directory <- getwd() 
directory_project <- paste0(directory,"/COVID")
output_directory <- paste0(directory_project,"/Output")
plot_directory <- paste0(directory_project,"/Output/Plots")

source(paste0(directory_project,"/functions.R"))

################################################################################
# Load data
################################################################################
data_directory <- paste0(directory_project,"/Data")

eff_excl_COV002 <- read.csv(file = paste0(data_directory,"/COV002InvestigatingA_DATA_2022-04-25_1735.csv"))
colnames(eff_excl_COV002)[1] <- "record_id"

COV002_groups <- read.csv(file = paste0(data_directory,"/cov002_studygroup_030221.csv"))

# Add back in the record_id
columns2 <- c("record_id","sc_repeat_pid","redcap_data_access_group",
              "study_group", "gp4_subgp", "vacc_date", "vacc_deviation",
              "booster_declined", "vacc_date_2", "vacc_deviation_2",
              "sars_cov2_serology", "excl_criterion19",
              "vaccination_complete", "randomisation_complete",
              "randomisation_date_time", "analysis_inclusion_r2",
              "primary_endpoint_r2", "non_primary_symptomatic_r2", 
              "pcr_asymptomatic_r2", "eos_date", "vacc_lotno", "vacc_lotno_2",
              "participant_randomised_to","vaccination_2_complete","start_date_event_r2",
              "primary_symptom_elig_r2","previous_participant_id","pospcr_asymptomatic_r2",
              "discon_auth_date","vacc_date_2_evb","cov_vacc_date","d28_pv_datetime_pb","d28_pv_datetime","d28_pv_datetime_pb_pb28",
              "d56_datetime","pospcr_cough_r2","pospcr_sob_r2",paste0("pospcr_symptoms_r2___",1:5),
              "prevac_hcw", "prevac_hcw_covid_contact", "prevac_hcw_patients","sc_age","sc_eth","sc_bmi_p2","sc_cv","sc_resp","sc_diabetes")

seropositive <- read.csv(paste0(data_directory,"/PPDGCL_Anti-SARS-Cov2_20201228.csv"))
head(seropositive)
seropositive <- seropositive[,c("SPECNUM","CNVRESC")]
head(seropositive)
g <- regexpr("-D",seropositive$SPECNUM)
f <- regexpr("-0",seropositive$SPECNUM)
seropositive$SPECNUM <- substring(seropositive$SPECNUM,first = f + 1, last = g - 1)
seropositive$SPECNUM <- as.character(as.numeric(seropositive$SPECNUM))
seropositive <- seropositive[which(!is.na(seropositive$SPECNUM)),]
head(seropositive)

##########################################
# Read the antibody data
# PPD
PPD <- read.xlsx(paste0(data_directory,"/Antibody_data/AZD1222_UNBLINDED_Final_20210730.xlsx"))
# Remove COV005
PPD <- PPD[which(substring(PPD$MBREFID,1,6)!="COV005"),]

Jenner <- read.xlsx(paste0(data_directory,"/Antibody_data/COV002 Std Spike data QC'd up to 18082021.xlsx"))
Monogram <- read.csv(paste0(data_directory,"/Antibody_data/Monogram_ID50_raw 12AUG21.csv"))
PHE <- read.csv(paste0(data_directory,"/Antibody_data/PHE Merged ND50 12MAY21.csv"))
PHE2 <- read.csv(paste0(data_directory,"/Antibody_data/PHE Merged ND80 12MAY21.csv"))

################################################################################
# COV002
################################################################################

# cut off date for the analysis
cut_off_date <- as.Date("2021-06-30")

all_ids <- eff_excl_COV002[,c("previous_participant_id","previous_participant_id_2","sc_repeat_pid","record_id")]
# eff_excl_COV002 <- eff_excl_COV002[,columns2]

previous_records <- eff_excl_COV002[which(!is.na(eff_excl_COV002$previous_participant_id)),"record_id"]
for (i in 1:length(previous_records)){
  eff_excl_COV002[which(eff_excl_COV002$record_id==previous_records[i]),"previous_participant_id"] <- unique(eff_excl_COV002[which((eff_excl_COV002$record_id==previous_records[i])&
                                                                                                                                     (!is.na(eff_excl_COV002$previous_participant_id)))[1],"previous_participant_id"])
}

eff_excl_COV002 <- eff_excl_COV002[which(eff_excl_COV002$study_group!=""),]
nrow(eff_excl_COV002)

eff_excl_COV002 <- eff_excl_COV002[which(eff_excl_COV002$vaccination_complete=="2"),]
eff_excl_COV002 <- eff_excl_COV002[which(eff_excl_COV002$vacc_date!=""),]

# Rename group 4b
eff_excl_COV002$study_group[which(eff_excl_COV002$gp4_subgp=="4b")] <- "4b"

eff_excl_COV002$vacc_date <- myfix_dates(eff_excl_COV002$vacc_date)
eff_excl_COV002$vacc_date_2 <- myfix_dates(eff_excl_COV002$vacc_date_2)
eff_excl_COV002$eos_date <- myfix_dates(eff_excl_COV002$eos_date)
eff_excl_COV002$start_date_event_r2 <- myfix_dates(eff_excl_COV002$start_date_event_r2)
eff_excl_COV002$discon_auth_date <- myfix_dates(eff_excl_COV002$discon_auth_date)

# Remove the participant who was not vaccinated.
eff_excl_COV002 <- eff_excl_COV002[which((eff_excl_COV002$vacc_deviation!="Patients declined to have the vaccine after it was drawn up.Was very anxious and would not proceed to have vaccine. Disposed of the vaccine in our GMO box 19/06/2020")|
                                           is.na(eff_excl_COV002$vacc_deviation)),]
eff_excl_COV002 <- eff_excl_COV002[which(difftime(cut_off_date,
                                                  as.Date(eff_excl_COV002$vacc_date),
                                                  units = "days") >= 0),]

nrow(eff_excl_COV002)
exclusion_reasons <- c("Non_randomised_groups",
                       "Baseline_serostatus",
                       "Single_dose",
                       "Vaccine_error",
                       "Case_leq_21_days_post_boost",
                       "Follow_up_leq_21_days_post_boost")
exclusion_reasons_cum <- paste0(exclusion_reasons,"_cum")
eff_excl_COV002[,c(exclusion_reasons,exclusion_reasons_cum)] <- NA 
# A value of TRUE means they meet the inclusion criteria, FALSE means they are excluded
####
# Define study group in more detail
eff_excl_COV002$study_group_detailed <- eff_excl_COV002$study_group

# groups_match <- match(eff_excl_COV002$sc_repeat_pid,
#                       COV002_groups$record_id,
#                       incomparables = NA)
# eff_excl_COV002$study_group_detailed[!is.na(groups_match)] <- COV002_groups$study_group_3[groups_match[!is.na(groups_match)]]
# Check the group of those not in the document
# eff_excl_COV002[which(eff_excl_COV002$study_group_detailed==""),]
# eff_excl_COV002$study_group_detailed[which(eff_excl_COV002$study_group_detailed=="")] <- eff_excl_COV002$study_group[which(eff_excl_COV002$study_group_detailed=="")]

# Define D28 and PB28 variables
eff_excl_COV002$d28_pv_datetime_pb <- myfix_dates(eff_excl_COV002$d28_pv_datetime_pb)
eff_excl_COV002$d28_pv_datetime <- myfix_dates(eff_excl_COV002$d28_pv_datetime)
eff_excl_COV002$d28_pv_datetime_pb_pb28 <- myfix_dates(eff_excl_COV002$d28_pv_datetime_pb_pb28)
eff_excl_COV002$d56_datetime <- myfix_dates(eff_excl_COV002$d56_datetime)
eff_excl_COV002$d35_datetime <- myfix_dates(eff_excl_COV002$d35_datetime)
eff_excl_COV002$d42_datetime <- myfix_dates(eff_excl_COV002$d42_datetime)
eff_excl_COV002$d90_datetime <- myfix_dates(eff_excl_COV002$d90_datetime)

eff_excl_COV002$D56 <- eff_excl_COV002$d56_datetime
eff_excl_COV002$vacc_date <- myfix_dates(eff_excl_COV002$vacc_date)
eff_excl_COV002$vacc_date_2 <- myfix_dates(eff_excl_COV002$vacc_date_2)

eff_excl_COV002$D28 <- eff_excl_COV002$d28_pv_datetime
eff_excl_COV002$PB28 <- eff_excl_COV002$d28_pv_datetime_pb_pb28

D28diffs <- round(difftime(eff_excl_COV002$PB28,
                           eff_excl_COV002$vacc_date_2,
                           units = "days"),0)
table(D28diffs) # Not all at PB28 - in fact some before PB21 which could be an issue if their antibody levels have not yet peaked.
eff_excl_COV002[which(D28diffs<21),c("vacc_date_2","PB28","D56","study_group","study_group_detailed")]

#####
# Non-randomised groups
eff_excl_COV002$study_group[which(eff_excl_COV002$study_group=="")] <- NA
eff_excl_COV002[which(is.element(eff_excl_COV002$study_group,c("11","12"))),"Non_randomised_groups"] <- F
eff_excl_COV002[which(!is.element(eff_excl_COV002$study_group,c("11","12"))),"Non_randomised_groups"] <- T
#####
# Baseline serostatus not negative
eff_excl_COV002$previous_participant_id <- myremove_characters(eff_excl_COV002$previous_participant_id)
eff_excl_COV002[which((!is.element(eff_excl_COV002$sc_repeat_pid,seropositive$SPECNUM[which(seropositive$CNVRESC=="Negative")])|
                         (is.element(eff_excl_COV002$previous_participant_id,seropositive$SPECNUM[which(seropositive$CNVRESC=="Negative")])))),"Baseline_serostatus"] <- F
eff_excl_COV002[which(is.element(eff_excl_COV002$sc_repeat_pid,seropositive$SPECNUM[which(seropositive$CNVRESC=="Negative")])|
                        (is.element(eff_excl_COV002$previous_participant_id,seropositive$SPECNUM[which(seropositive$CNVRESC=="Negative")]))),"Baseline_serostatus"] <- T
#####
# Removing participants who are yet to be boosted
eff_excl_COV002[,"Single_dose"] <- F
eff_excl_COV002[which((!is.na(eff_excl_COV002$vacc_date_2))&
                        (difftime(as.Date(cut_off_date),
                                  as.Date(eff_excl_COV002$vacc_date_2),
                                  units = "days") > 0) &
                        eff_excl_COV002$vaccination_2_complete==2),"Single_dose"] <- T
#####
# Incorrectly vaccinated
ChAdOx_lots <- c("A0","AO","K0","K.","01","02","O2","a0","ao","aO","Ao","k0","k.","K,","K-","K/","KO","AZ","az","41")
eff_excl_COV002$vacc_lotno[which(eff_excl_COV002$vacc_lotno=="K/0007")] <- "K0007"
eff_excl_COV002$vacc_lotno[which(eff_excl_COV002$vacc_lotno=="A03373/ 228021")] <- "A03373"
eff_excl_COV002$vacc_lotno[which(eff_excl_COV002$vacc_lotno=="22819-A03373")] <- "A03373"
eff_excl_COV002$vacc_lotno_2[which(eff_excl_COV002$vacc_lotno_2=="A03373/22814")] <- "A03373"
eff_excl_COV002$vacc_lotno_2[which(eff_excl_COV002$vacc_lotno_2=="A03373/22815")] <- "A03373"
eff_excl_COV002$vacc_lotno_2[which(eff_excl_COV002$vacc_lotno_2=="A03373/22817")] <- "A03373"
eff_excl_COV002$vacc_lotno_2[which(eff_excl_COV002$vacc_lotno_2=="A03373/22817")] <- "A03373"

f_lot <- regexpr("/", eff_excl_COV002$vacc_lotno)
f_lot_2 <- regexpr("/", eff_excl_COV002$vacc_lotno_2)
which_slash <- which(f_lot!="-1")
which_slash_2 <- which(f_lot_2!="-1")
vacc_lotno_2_slash <- eff_excl_COV002$vacc_lotno_2[which(f_lot_2!="-1")]
vacc_lotno_slash <- eff_excl_COV002$vacc_lotno[which(f_lot!="-1")]
f_lot <- regexpr("/", vacc_lotno_slash)
f_lot_2 <- regexpr("/", vacc_lotno_2_slash)
g_lot <- nchar(vacc_lotno_slash)
g_lot_2 <- nchar(vacc_lotno_2_slash)

eff_excl_COV002$vacc_lotno[which_slash] <- substring(vacc_lotno_slash, first = f_lot + 1, last = g_lot)
eff_excl_COV002$vacc_lotno <- myremove_characters(eff_excl_COV002$vacc_lotno)
eff_excl_COV002$vacc_lotno_2[which_slash_2] <- substring(vacc_lotno_2_slash, first = f_lot_2 + 1, last = g_lot_2)
eff_excl_COV002$vacc_lotno_2 <- myremove_characters(eff_excl_COV002$vacc_lotno_2)

eff_excl_COV002$vacc_lotno_first_two <- substring(eff_excl_COV002$vacc_lotno, first = 1, last = 2)
eff_excl_COV002$prime <- "Control"
eff_excl_COV002$prime[which(is.element(eff_excl_COV002$vacc_lotno_first_two,ChAdOx_lots))]  <- "ChAdOx1"
eff_excl_COV002$prime[which((is.na(eff_excl_COV002$vacc_date))|
                              (difftime(as.Date(cut_off_date),
                                        as.Date(eff_excl_COV002$vacc_date),
                                        units = "days") < 0))] <- NA
# eff_excl_COV002[which((eff_excl_COV002$prime=="ChAdOx1")&
#                         (eff_excl_COV002$participant_randomised_to=="MenACWY")),c("sc_repeat_pid","participant_randomised_to","vacc_lotno","vacc_deviation")]
eff_excl_COV002$vacc_lotno_2_first_two <- substring(eff_excl_COV002$vacc_lotno_2, first = 1, last = 2)
eff_excl_COV002$boost <- "Control"
eff_excl_COV002$boost[which(is.element(eff_excl_COV002$vacc_lotno_2_first_two,ChAdOx_lots))]  <- "ChAdOx1"
eff_excl_COV002$boost[which(is.na(eff_excl_COV002$vacc_date_2))] <- NA
eff_excl_COV002$boost[which((is.na(eff_excl_COV002$vacc_date_2))|
                              (difftime(as.Date(cut_off_date),
                                        as.Date(eff_excl_COV002$vacc_date_2),
                                        units = "days") < 0))] <- NA
eff_excl_COV002$boost[which(is.na(eff_excl_COV002$vacc_lotno_2))] <- NA
eff_excl_COV002$boost[which(eff_excl_COV002$vaccination_2_complete!=2)] <- NA

all(!eff_excl_COV002[which(is.na(eff_excl_COV002$boost)),exclusion_reasons[4]]) # All those with boost as NA have already been excluded as not boosted, as you'd hope 

eff_excl_COV002$participant_randomised_to[which(is.na(eff_excl_COV002$participant_randomised_to))] <- ""
# Create arm variable of ChAdOx or Control
arm <- eff_excl_COV002$participant_randomised_to
g <- regexpr("_", arm) # Gives a vector of locations of the underscore
f <- nchar(arm)
arm <- substring(arm, first = g + 1, last = f) # Takes all characters after the underscore
eff_excl_COV002$arm <- arm
eff_excl_COV002$arm[which(eff_excl_COV002$arm=="")] <- NA
# Create arm2 variable as above based on env_participant_randomised 
arm2 <- eff_excl_COV002$env_participant_randomised
g <- regexpr("_", arm2) # Gives a vector of locations of the underscore
f <- nchar(arm2)
arm2 <- substring(arm2, first = g + 1, last = f) # Takes all characters after the underscore
eff_excl_COV002$arm2 <- arm2
eff_excl_COV002$arm[which(is.na(eff_excl_COV002$arm))] <- eff_excl_COV002$arm2[which(is.na(eff_excl_COV002$arm))]
# Groups 11 and 12 were non-randomised, all received ChAdOx1
eff_excl_COV002$arm[which(is.element(eff_excl_COV002$study_group,c(11,12)))] <- "ChAdOx1"

eff_excl_COV002$arm[which(eff_excl_COV002$arm=="MenACWY")] <- "Control"

SD_SD <- c("5b","5c","5d","5e","6","7","8","9","10","11","12")

# Removing all those who received a different prime to boost in primary
eff_excl_COV002[which((eff_excl_COV002$prime==eff_excl_COV002$boost)|is.na(eff_excl_COV002$boost)),"Vaccine_error"] <- T
eff_excl_COV002[which((eff_excl_COV002$prime!=eff_excl_COV002$boost)),"Vaccine_error"] <- F

#####
# Removing cases <= 21 days since boost
eff_excl_COV002$start_date_event_r2 <- myfix_dates(eff_excl_COV002$start_date_event_r2)
# eff_excl_COV002$PB28[which(is.na(eff_excl_COV002$d28_pv_datetime_pb_pb28))] <- eff_excl_COV002$d56_datetime[which(is.na(eff_excl_COV002$d28_pv_datetime_pb_pb28))]
eff_excl_COV002$PB28 <- myfix_dates(eff_excl_COV002$PB28)

eff_excl_COV002[which((difftime(as.Date(eff_excl_COV002$start_date_event_r2),
                                as.Date(eff_excl_COV002$vacc_date_2),
                                units = "days") > 21)|
                        (difftime(as.Date(eff_excl_COV002$start_date_event_r2),
                                  as.Date(cut_off_date),
                                  units = "days") > 0)|
                        (is.na(eff_excl_COV002$start_date_event_r2))),"Case_leq_21_days_post_boost"] <- T
eff_excl_COV002[which((difftime(as.Date(eff_excl_COV002$start_date_event_r2),
                                as.Date(eff_excl_COV002$vacc_date_2),
                                units = "days") <= 21)&
                        (difftime(as.Date(eff_excl_COV002$start_date_event_r2),
                                  as.Date(cut_off_date),
                                  units = "days") < 0)),"Case_leq_21_days_post_boost"] <- F
#####
# Exclude any participants who are yet to accrue 21 days post boost
eff_excl_COV002$censoring_date <- pmin(as.Date(cut_off_date),as.Date(eff_excl_COV002$eos_date),as.Date(eff_excl_COV002$discon_auth_date),na.rm=T)
eff_excl_COV002[which(as.numeric(difftime(as.Date(eff_excl_COV002$censoring_date),
                                          as.Date(eff_excl_COV002$vacc_date_2),
                                          units = "days"))>21),"Follow_up_leq_21_days_post_boost"] <- T
eff_excl_COV002[which(as.numeric(difftime(as.Date(eff_excl_COV002$censoring_date),
                                          as.Date(eff_excl_COV002$vacc_date_2),
                                          units = "days"))<=21),"Follow_up_leq_21_days_post_boost"] <- F
#####
# cor2dose_exclusion gives the exclusion reason
# exclusion_reasons_cum gives whether a person is included in all inclusion criterion so far.
eff_excl_COV002$cor2dose_exclusion <- "include"
for(i in 1:length(exclusion_reasons)){
  eff_excl_COV002[,exclusion_reasons_cum[i]] <- apply(as.matrix(eff_excl_COV002[,exclusion_reasons[1:i]],ncol=i),1,all)
}
eff_excl_COV002[!eff_excl_COV002[,exclusion_reasons_cum[1]],"cor2dose_exclusion"] <- exclusion_reasons[1]
for(i in 2:length(exclusion_reasons)){
  eff_excl_COV002[(!eff_excl_COV002[,exclusion_reasons_cum[i]]) & eff_excl_COV002[,exclusion_reasons_cum[i-1]],"cor2dose_exclusion"] <- exclusion_reasons[i]
}
eff_excl_COV002$cor2dose_exclusion <- factor(eff_excl_COV002$cor2dose_exclusion, levels = c("include",exclusion_reasons))
table(eff_excl_COV002[,c("arm","cor2dose_exclusion")],useNA="ifany")
exclusion_levels <- c("include","Non-randomised group","Baseline serostatus not negative","Received single dose of vaccine","Vaccine error","Case before 21 days post second dose","Loss of follow up before 21 days post second dose")
levels(eff_excl_COV002$cor2dose_exclusion) <- exclusion_levels

eff_excl_COV002$cor2dose_prerandomised_exclusion <- NA
eff_excl_COV002$cor2dose_prerandomised_exclusion[which(is.element(eff_excl_COV002$cor2dose_exclusion,exclusion_levels[2]))] <- eff_excl_COV002$cor2dose_exclusion[which(is.element(eff_excl_COV002$cor2dose_exclusion,exclusion_levels[2]))]
eff_excl_COV002$cor2dose_prerandomised_exclusion[which(eff_excl_COV002$cor2dose_prerandomised_exclusion==2)] <- exclusion_levels[2]
eff_excl_COV002$cor2dose_prevaccine_exclusion <- eff_excl_COV002$cor2dose_exclusion
eff_excl_COV002$cor2dose_prevaccine_exclusion[which(is.element(eff_excl_COV002$cor2dose_prevaccine_exclusion,exclusion_levels[c(1,2,6,7)]))] <- NA
eff_excl_COV002$cor2dose_prevaccine_exclusion <-  as.character(eff_excl_COV002$cor2dose_prevaccine_exclusion)
eff_excl_COV002$cor2dose_postvaccine_exclusion <- NA
eff_excl_COV002$cor2dose_postvaccine_exclusion[which(is.element(eff_excl_COV002$cor2dose_exclusion,exclusion_levels[c(6,7)]))] <- eff_excl_COV002$cor2dose_exclusion[which(is.element(eff_excl_COV002$cor2dose_exclusion,exclusion_levels[c(6,7)]))]
eff_excl_COV002$cor2dose_postvaccine_exclusion[which(eff_excl_COV002$cor2dose_postvaccine_exclusion==6)] <- exclusion_levels[6]
eff_excl_COV002$cor2dose_postvaccine_exclusion[which(eff_excl_COV002$cor2dose_postvaccine_exclusion==7)] <- exclusion_levels[7]

eff_excl_COV002$arm_switch <- eff_excl_COV002$arm
eff_excl_COV002$arm_switch <- factor(eff_excl_COV002$arm_switch, levels = c("Control","ChAdOx1"))
levels(eff_excl_COV002$arm_switch) <- c("Allocated to MenACWY control vaccine","Allocated to ChAdOx1 nCoV-19 vaccine")

LD_LD <- c("4b")
LD_SD <- c("1","2","4","5","5f")
SD_SD <- c("5b","5c","5d","5e","6","7","8","9","10","11","12")

eff_excl_COV002$boost_2 <- NA
eff_excl_COV002$boost_2[which((eff_excl_COV002$boost=="Control"))] <- "Control"
eff_excl_COV002$boost_2[which((eff_excl_COV002$boost=="ChAdOx1")&
                                (is.element(eff_excl_COV002$study_group,c(SD_SD,"4","5","5f"))))] <- "SD"
eff_excl_COV002$boost_2[which((eff_excl_COV002$boost=="ChAdOx1")&
                                (is.element(eff_excl_COV002$study_group,LD_LD)))] <- "LD"
eff_excl_COV002$boost_2[which((eff_excl_COV002$boost=="ChAdOx1")&
                                (is.element(eff_excl_COV002$study_group,c("1","2")))&
                                (substring(eff_excl_COV002$participant_randomised_to, first = 1, last = 3)=="Two"))] <- "LD"
eff_excl_COV002$boost_2[which((eff_excl_COV002$boost=="ChAdOx1")&
                                (is.element(eff_excl_COV002$study_group,c("1","2")))&
                                (substring(eff_excl_COV002$participant_randomised_to, first = 1, last = 3)=="Sin"))] <- "SD"
eff_excl_COV002$boost_2[which(is.na(eff_excl_COV002$boost))] <- "Not vaccinated"
eff_excl_COV002[which(is.na(eff_excl_COV002$boost_2)),c("boost","study_group","study_group_detailed","vacc_date_2","participant_randomised_to")]

eff_excl_COV002$prime_2 <- NA
eff_excl_COV002$prime_2[which((eff_excl_COV002$prime=="Control"))] <- "Control"
eff_excl_COV002$prime_2[which((eff_excl_COV002$prime=="ChAdOx1")&
                                (is.element(eff_excl_COV002$study_group,SD_SD)))] <- "SD"
eff_excl_COV002$prime_2[which((eff_excl_COV002$prime=="ChAdOx1")&
                                (is.element(eff_excl_COV002$study_group,c(LD_LD,LD_SD))))] <- "LD"
eff_excl_COV002$prime_2[which(is.na(eff_excl_COV002$prime))] <- "Not vaccinated"

eff_excl_COV002$As_vaccinated_group[which(is.na(eff_excl_COV002$boost_2))] <- eff_excl_COV002$prime_2[which(is.na(eff_excl_COV002$boost_2))]
eff_excl_COV002$As_vaccinated_group[which(!is.na(eff_excl_COV002$boost_2))] <- paste0(eff_excl_COV002$prime_2[which(!is.na(eff_excl_COV002$boost_2))],"/",eff_excl_COV002$boost_2[which(!is.na(eff_excl_COV002$boost_2))])

LD_LD <- c("4b")
LD_SD <- c("1","2","4","5","5f")
SD_SD <- c("5b","5c","5d","5e","6","7","8","9","10","11","12")

eff_excl_COV002$LDSD <- 0
eff_excl_COV002$LDSD[which((!is.na(eff_excl_COV002$boost))&
                             (is.element(eff_excl_COV002$study_group,c("4","5","5f"))))] <- 1
eff_excl_COV002$LDSD[which((!is.na(eff_excl_COV002$boost))&
                             (is.element(eff_excl_COV002$study_group,c("1","2")))&
                             (substring(eff_excl_COV002$participant_randomised_to, first = 1, last = 3)=="Sin"))] <- 1
eff_excl_COV002$SDSD <- 0
eff_excl_COV002$SDSD[which((!is.na(eff_excl_COV002$boost))&
                             (is.element(eff_excl_COV002$study_group,SD_SD)))] <- 1
eff_excl_COV002$LDLD <- 0
eff_excl_COV002$LDLD[which((!is.na(eff_excl_COV002$boost))&
                             (is.element(eff_excl_COV002$study_group,LD_LD)))] <- 1
eff_excl_COV002$LDLD[which((!is.na(eff_excl_COV002$boost))&
                             (is.element(eff_excl_COV002$study_group,c("1","2")))&
                             (substring(eff_excl_COV002$participant_randomised_to, first = 1, last = 3)=="Two"))] <- 1
# Everyone who received a boost fits exactly one of these categories
table(rowSums(eff_excl_COV002[which(!is.na(eff_excl_COV002$boost)),c("LDLD","SDSD","LDSD")]), useNA="ifany")
eff_excl_COV002$cor2dose_schedule <- NA
eff_excl_COV002$cor2dose_schedule[which(eff_excl_COV002$LDLD==1)] <- "LDLD"
eff_excl_COV002$cor2dose_schedule[which(eff_excl_COV002$LDSD==1)] <- "LDSD"
eff_excl_COV002$cor2dose_schedule[which(eff_excl_COV002$SDSD==1)] <- "SDSD"

table(eff_excl_COV002$boost_2,useNA="ifany")
table(eff_excl_COV002$boost,useNA="ifany")
table(eff_excl_COV002$prime_2,useNA="ifany")
table(eff_excl_COV002$prime,useNA="ifany")

###############################################################################
# Merge the data
###############################################################################
columns <-
  c("Study","sc_repeat_pid","Serostatus",
    "participant_randomised_to",
    "vacc_lotno","vacc_lotno_2",
    "As_vaccinated_group",
    "LDSD",
    "SDSD",
    "D28",
    "PB28",
    "D56",
    "cor2dose_cohort",
    "cor2dose_outcome",
    "cor2dose_start_date",
    "cor2dose_end_date",
    "cor2dose_offset_days",
    "cor2dose_exclusion",
    "cor2dose_positive",
    "cor2dose_primary",
    "cor2dose_sob",
    "cor2dose_nosob",
    "cor2dose_sob_cough",
    "cor2dose_nosob_cough",
    "cor2dose_no_symp",
    "cor2dose_symp_3ormore",
    "cor2dose_symp_3less",
    "cor2dose_schedule",
    exclusion_reasons,
    exclusion_reasons_cum,
    "study_group", "prime", "boost",
    "primary_endpoint_r2","start_date_event_r2","non_primary_symptomatic_r2",
    "pcr_asymptomatic_r2","eos_date","vacc_date","vacc_date_2","primary_symptom_elig_r2","pospcr_asymptomatic_r2",
    "analysis_inclusion_r2","primary.exclusion", "ITT.exclusion",
    "vaccination_2_complete","discon_auth_date","vacc_date_2_evb","cov_vacc_date","pospcr_cough_r2","pospcr_sob_r2",paste0("pospcr_symptoms_r2___",1:5),
    "prevac_hcw", "prevac_hcw_covid_contact", "prevac_hcw_patients","sc_age","sc_eth","sc_bmi_p2","sc_cv","sc_resp","sc_diabetes")

excl_data <- eff_excl_COV002
excl_data[,columns[which(!is.element(columns,colnames(excl_data)))]] <- NA

# Fix dates and remove characters
excl_data$start_date_event_r2 <- myfix_dates(excl_data$start_date_event_r2)
excl_data$discon_auth_date <- myfix_dates(excl_data$discon_auth_date)
excl_data$vacc_date_2_evb <- myfix_dates(excl_data$vacc_date_2_evb)
excl_data$cov_vacc_date <- myfix_dates(excl_data$cov_vacc_date)
excl_data$primary_endpoint_r2 <- myremove_characters(excl_data$primary_endpoint_r2)
excl_data$primary_symptom_elig_r2 <- myremove_characters(excl_data$primary_endpoint_r2)
excl_data$pcr_asymptomatic_r2 <- myremove_characters(excl_data$pcr_asymptomatic_r2)

# Two-dose correlates

# cohort
excl_data$cor2dose_cohort <- excl_data$cor2dose_exclusion=="include"

# First define if they are a case or not
# (if event occurs strictly after 7 days post PB28 and before or equal to end date)
# COV002 only included in correlates
excl_data$cor2dose_outcome <- "Negative"
excl_data$cor2dose_positive <- 0
excl_data$cor2dose_positive[which(difftime(as.Date(excl_data$start_date_event_r2),
                                           as.Date(excl_data$vacc_date_2) + 21,
                                           units = "days") > 0 &
                                    difftime(pmin(as.Date(excl_data$eos_date),
                                                  as.Date(cut_off_date),
                                                  as.Date(excl_data$discon_auth_date),
                                                  na.rm=TRUE),
                                             as.Date(excl_data$start_date_event_r2),
                                             units = "days") >= 0)] <- 1
excl_data$cor2dose_outcome[which(excl_data$cor2dose_positive==1)] <- "Unknown"
# Then define type
excl_data$cor2dose_outcome[which(excl_data$cor2dose_positive==1 &
                                   excl_data$primary_endpoint_r2==1)] <- "Primary"
excl_data$cor2dose_outcome[which(excl_data$cor2dose_positive==1 & 
                                   excl_data$non_primary_symptomatic_r2==1 &
                                   excl_data$analysis_inclusion_r2==1 &
                                   excl_data$pospcr_asymptomatic_r2!=2)] <- "Non-Primary"
excl_data$cor2dose_outcome[which(excl_data$cor2dose_positive==1 & 
                                   excl_data$pcr_asymptomatic_r2==1)] <- "Asymptomatic"
# Ensure those who aren't positive in the at-risk period are counted as negative.
excl_data$cor2dose_outcome[which(excl_data$cor2dose_positive==0)] <- "Negative"
excl_data$cor2dose_primary <- as.numeric(excl_data$cor2dose_outcome=="Primary")

# Offset days - at risk from 21 days post boost.
excl_data$cor2dose_start_date <- as.Date(excl_data$vacc_date_2) + 21
excl_data$cor2dose_end_date <- pmin(as.Date(excl_data$start_date_event_r2),
                                    as.Date(excl_data$eos_date),
                                    as.Date(cut_off_date),
                                    as.Date(excl_data$discon_auth_date),
                                    na.rm = TRUE)
excl_data$cor2dose_offset_days <- as.numeric(difftime(excl_data$cor2dose_end_date,
                                                      excl_data$cor2dose_start_date,
                                                      units = "days"))

# cor2dose other variables
# sob
excl_data$cor2dose_sob <- excl_data$pospcr_sob_r2
excl_data$cor2dose_nosob[which(excl_data$pospcr_sob_r2==1)] <- 0
excl_data$cor2dose_nosob[which(excl_data$pospcr_asymptomatic_r2==1)] <- 1
excl_data$cor2dose_nosob[which(excl_data$pospcr_asymptomatic_r2==0 &
                                 excl_data$pospcr_sob_r2 == 0)] <- 1
excl_data$cor2dose_sob[which(excl_data$cor2dose_outcome!="Primary")] <- 0
excl_data$cor2dose_nosob[which(excl_data$cor2dose_outcome!="Primary")] <- 0

# sob_cough
excl_data$cor2dose_sob_cough <- as.numeric(excl_data$pospcr_sob_r2 + excl_data$pospcr_cough_r2 == 2)
excl_data$cor2dose_nosob_cough <- 0
excl_data$cor2dose_nosob[which(excl_data$pospcr_sob_r2==1 & excl_data$pospcr_cough_r2==1)] <- 0
excl_data$cor2dose_nosob_cough[which(excl_data$pospcr_asymptomatic_r2==1)] <- 1
excl_data$cor2dose_nosob_cough[which(excl_data$pospcr_asymptomatic_r2==0 &
                                       excl_data$pospcr_sob_r2 == 0)] <- 1
excl_data$cor2dose_nosob_cough[which(excl_data$pospcr_asymptomatic_r2==0 &
                                       excl_data$pospcr_cough_r2 == 0)] <- 1
excl_data$cor2dose_sob_cough[which(excl_data$cor2dose_outcome!="Primary")] <- 0
excl_data$cor2dose_nosob_cough[which(excl_data$cor2dose_outcome!="Primary")] <- 0

# no_symp number of symptoms
excl_data$cor2dose_no_symp <- rowSums(excl_data[,paste0("pospcr_symptoms_r2___",1:5)])
excl_data$cor2dose_no_symp[which(excl_data$cor2dose_outcome!="Primary")] <- NA

# symp_3ormore
excl_data$cor2dose_symp_3ormore <- as.numeric(excl_data$cor2dose_no_symp >= 3)
excl_data$cor2dose_symp_3less <- as.numeric(excl_data$cor2dose_no_symp < 3)
excl_data$cor2dose_symp_3ormore[which(excl_data$cor2dose_outcome!="Primary")] <- 0
excl_data$cor2dose_symp_3less[which(excl_data$cor2dose_outcome!="Primary")] <- 0

excl_data[!excl_data$cor2dose_cohort,c("cor2dose_positive",
                                       "cor2dose_sob",
                                       "cor2dose_nosob",
                                       "cor2dose_sob_cough",
                                       "cor2dose_nosob_cough",
                                       "cor2dose_no_symp",
                                       "cor2dose_symp_3ormore",
                                       "cor2dose_symp_3less",
                                       "cor2dose_outcome",
                                       "cor2dose_start_date",
                                       "cor2dose_end_date",
                                       "cor2dose_offset_days",
                                       "cor2dose_primary")] <- NA
# excl_data <- excl_data[excl_data$Any_dose_for_safety,]

# age_group
excl_data$age_group <- NA
excl_data$age_group[which(excl_data$sc_age>=18)] <- "18-55"
excl_data$age_group[which(excl_data$sc_age>=56)] <- "56-69"
excl_data$age_group[which(excl_data$sc_age>=70)] <- ">=70"

# prime-boost interval
excl_data$cor2dose_interval_cont <- NA
excl_data$cor2dose_interval <- NA
excl_data$cor2dose_interval_cont <- as.numeric(difftime(as.Date(excl_data$vacc_date_2),
                                                        as.Date(excl_data$vacc_date),
                                                        units = "days"))
excl_data$cor2dose_interval[which(excl_data$cor2dose_interval_cont < 6*7)] <- "<6"
excl_data$cor2dose_interval[which(excl_data$cor2dose_interval_cont >= 6*7)] <- "6-8"
excl_data$cor2dose_interval[which(excl_data$cor2dose_interval_cont >= 9*7)] <- "9-11"
excl_data$cor2dose_interval[which(excl_data$cor2dose_interval_cont >= 12*7)] <- ">=12"

# ethnicity
excl_data$cor2dose_non_white <- NA
excl_data$cor2dose_non_white[which(is.element(excl_data$sc_eth,c(1,2,3,4)))] <- 0
excl_data$cor2dose_non_white[which(is.element(excl_data$sc_eth,5:19))] <- 1

# bmi_geq_30
# Note one individual has missing BMI. We impute theirs with the mean value
excl_data$sc_bmi_p2[which(is.na(excl_data$sc_bmi_p2))] <- mean(excl_data$sc_bmi_p2[which(excl_data$cor2dose_cohort)],na.rm=T)
excl_data$cor2dose_bmi_geq_30 <- NA
excl_data$cor2dose_bmi_geq_30[which(excl_data$sc_bmi_p2<30)] <- 0
excl_data$cor2dose_bmi_geq_30[which(excl_data$sc_bmi_p2>=30)] <- 1

# cor2dose_comorbidities
excl_data$cor2dose_comorbidities <- NA
excl_data$cor2dose_comorbidities[which(rowSums(excl_data[,c("sc_cv","sc_resp","sc_diabetes")])==0)] <- 0
excl_data$cor2dose_comorbidities[which(rowSums(excl_data[,c("sc_cv","sc_resp","sc_diabetes")])>0)] <- 1

# cor2dose_hcw_status
excl_data$cor2dose_hcw_status <- NA
excl_data$cor2dose_hcw_status[which(excl_data$prevac_hcw==0)] <- "Not_hcw"
excl_data$cor2dose_hcw_status[which(excl_data$prevac_hcw==1 & is.element(excl_data$prevac_hcw_covid_contact,c(1,2)))] <- "Less_than_1_covid"
excl_data$cor2dose_hcw_status[which(excl_data$prevac_hcw==1 & is.element(excl_data$prevac_hcw_covid_contact,c(3,4,5)))] <- "More_than_1_covid"
# How to treat the missing data - have assumed those with missing healthcare worker status are not healthcare workers.
excl_data$cor2dose_hcw_status[which(is.na(excl_data$prevac_hcw))] <- "Not_hcw"
excl_data$cor2dose_hcw_status[which(excl_data$prevac_hcw==1 & is.na(excl_data$prevac_hcw_covid_contact))] <- "Less_than_1_covid"
excl_data$cor2dose_hcw_status[which(excl_data$prevac_hcw==1 & excl_data$prevac_hcw_covid_contact==6)] <- "Less_than_1_covid"

excl_data[is.na(excl_data$cor2dose_hcw_status[excl_data$cor2dose_cohort]),c("prevac_hcw","prevac_hcw_covid_contact")]

# As_vaccinated_arm_2
excl_data$As_vaccinated_arm_2 <- "Control"
excl_data$As_vaccinated_arm_2[which(excl_data$prime=="ChAdOx1")] <- "ChAdOx1"
excl_data$As_vaccinated_arm_2[which(excl_data$boost=="ChAdOx1")] <- "ChAdOx1"
excl_data <- excl_data[order(excl_data$sc_repeat_pid, decreasing = FALSE),]
dim(excl_data)
# Write name with current commit name on from Github.
# write.csv(excl_data,
#            paste0(directory,"/2_Data_for_Analysis/COV002_joint_data_",system("git rev-parse HEAD", intern=TRUE),".csv"),
#            row.names = FALSE)

primary_exclusions <- rbind(c("Total vaccinated", nrow(excl_data)),
                            cbind(levels(excl_data$cor2dose_exclusion),table(excl_data$cor2dose_exclusion)))
primary_exclusions[which(primary_exclusions=="include")] <- "Total included"

###############################################################################
# Add in antibodies
###############################################################################

##########################################
# Data cleaning

# Check the ids in the PPD data
PPD$pid_unreliable <- PPD$SUBJID
# Remove the first digit if it is 0
PPD$pid_unreliable[which(substring(PPD$SUBJID,first = 1, last = 1)=="0")] <- substring(PPD$pid_unreliable[which(substring(PPD$SUBJID,first = 1, last = 1)=="0")], first = 2, last = nchar(PPD$pid_unreliable[which(substring(PPD$SUBJID,first = 1, last = 1)=="0")]))
# Remove the first digit if it is 0
PPD$pid_unreliable[which(substring(PPD$pid_unreliable,first = 1, last = 1)=="0")] <- substring(PPD$pid_unreliable[which(substring(PPD$pid_unreliable,first = 1, last = 1)=="0")], first = 2, last = nchar(PPD$pid_unreliable[which(substring(PPD$pid_unreliable,first = 1, last = 1)=="0")]))

# Create the pid from the other column
PPD$sc_repeat_pid <- PPD$MBREFID
initial_dash <- regexpr("-",PPD$sc_repeat_pid)
PPD$sc_repeat_pid <- substring(PPD$sc_repeat_pid,initial_dash+1,nchar(PPD$sc_repeat_pid))
second_dash <- regexpr("-",PPD$sc_repeat_pid)
PPD$sc_repeat_pid <- substring(PPD$sc_repeat_pid,first = 1, last = second_dash - 1)
PPD$sc_repeat_pid[substring(PPD$sc_repeat_pid,1,1)=="0"] <- substring(PPD$sc_repeat_pid[substring(PPD$sc_repeat_pid,1,1)=="0"],2,nchar(PPD$sc_repeat_pid[substring(PPD$sc_repeat_pid,1,1)=="0"]))
PPD$sc_repeat_pid[substring(PPD$sc_repeat_pid,1,1)=="0"] <- substring(PPD$sc_repeat_pid[substring(PPD$sc_repeat_pid,1,1)=="0"],2,nchar(PPD$sc_repeat_pid[substring(PPD$sc_repeat_pid,1,1)=="0"]))
PPD$sc_repeat_pid[which(!is.element(substring(PPD$sc_repeat_pid,1,1),1:3))] <- NA
PPD$pid_unreliable[which(substring(PPD$pid_unreliable,1,1)=="0")] <- NA
PPD$pid_unreliable[which(nchar(PPD$pid_unreliable)<7)] <- NA

PPD <- PPD[which(substring(PPD$sc_repeat_pid,1,1) == "2"),]

# The following data need to be queried:
PPD[which((PPD$pid_unreliable)!=(PPD$sc_repeat_pid)),]

#########
# Merge correlates and PPD data
joint_correlates <- excl_data
# Keep only those in 2 dose cohort
joint_correlates <- joint_correlates[which(joint_correlates$cor2dose_cohort),]
# Only keep certain columns
# Define start and end_time variables
time_zero <- min(as.Date(excl_data$vacc_date_2[which(excl_data$cor2dose_cohort)])) + 21
joint_correlates$start_time <- as.numeric(difftime(as.Date(joint_correlates$cor2dose_start_date),
                                                   as.Date(time_zero),
                                                   units = "days"))
joint_correlates$end_time <- as.numeric(difftime(as.Date(joint_correlates$cor2dose_end_date),
                                                 as.Date(joint_correlates$cor2dose_start_date),
                                                 units = "days"))
joint_correlates$site <- joint_correlates$redcap_data_access_group

timepoint_names <- c("PB+28","PB+90","PB+182","D35","D42","D56","D90","D182")
antibody_times <- c("d28_pv_datetime_pb_pb28","d90_datetime_pb90","d182_datetime_pb182","d35_datetime","d42_datetime","d56_datetime","d90_datetime","d182_datetime")
for (time_var in antibody_times){
  joint_correlates[,time_var] <- myfix_dates(joint_correlates[,time_var])
}
# Define joint_correlates variables and set factors
joint_correlates$age_group <- factor(joint_correlates$age_group, levels = c("18-55","56-69",">=70"))
joint_correlates$cor2dose_schedule <- factor(joint_correlates$cor2dose_schedule, levels = c("SDSD","LDLD","LDSD")) # See below
# Time in weeks between first and second dose of vaccine
# Set baseline to be >=12 as this is what has mainly been used in practise.
joint_correlates$cor2dose_interval <- factor(joint_correlates$cor2dose_interval, levels = c(">=12","9-11","6-8","<6"))

# Use only the prime schedule (first dose given LD = low dose, SD = standard dose) instead of overall schedule
# (both doses, LDSD = low dose followed by standard dose etc.)
# We only use the prime schedule as there are only two antibody observations for LDLD individuals, which
# leads to serious problems with instability in estimation.
joint_correlates$cor2dose_primeschedule <- as.character(joint_correlates$cor2dose_schedule)
joint_correlates$cor2dose_primeschedule[which(is.element(joint_correlates$cor2dose_primeschedule,c("LDLD","LDSD")))] <- "LD"
joint_correlates$cor2dose_primeschedule[which(joint_correlates$cor2dose_primeschedule=="SDSD")] <- "SD"
joint_correlates$cor2dose_primeschedule <- factor(joint_correlates$cor2dose_primeschedule,levels=c("SD","LD"))

# Define the combined interval variable, splitting the time between first and second dose into only <8 weeks and >=9 weeks.
joint_correlates$cor2dose_interval_comb <- as.character(joint_correlates$cor2dose_interval)
joint_correlates$cor2dose_interval_comb[which(is.element(joint_correlates$cor2dose_interval,c("<6","6-8")))] <- "<8"
joint_correlates$cor2dose_interval_comb[which(is.element(joint_correlates$cor2dose_interval,c("9-11",">=12")))] <- ">=9"
joint_correlates$cor2dose_interval_comb <- factor(joint_correlates$cor2dose_interval_comb, levels = c(">=9","<8"))

# Convert variables to factor
# gender
joint_correlates$sc_gender <- factor(joint_correlates$sc_gender, levels = c(1,2))
levels(joint_correlates$sc_gender) <- c("Male","Female")
# non-white
joint_correlates$cor2dose_non_white <- factor(joint_correlates$cor2dose_non_white, levels = c(0,1))
levels(joint_correlates$cor2dose_non_white) <- c("White","Other")
# comorbidities
joint_correlates$cor2dose_comorbidities <- factor(joint_correlates$cor2dose_comorbidities, levels = c(0,1))
levels(joint_correlates$cor2dose_comorbidities) <- c("None","Comorbidity")
# BMI
joint_correlates$cor2dose_bmi_geq_30 <- factor(joint_correlates$cor2dose_bmi_geq_30, levels = c(0,1))
levels(joint_correlates$cor2dose_bmi_geq_30) <- c("BMI_less_than_30","BMI_geq_30")
# Outcome
joint_correlates$cor2dose_outcome <- factor(joint_correlates$cor2dose_outcome, levels = c("Negative","Primary","Non-Primary","Asymptomatic","Unknown"))
# Positive
# Indicator variable
joint_correlates$cor2dose_positive_ind <- factor(joint_correlates$cor2dose_positive, levels = c(0,1))
# and text variable
joint_correlates$cor2dose_positive <- factor(joint_correlates$cor2dose_positive, levels = c(0,1))
levels(joint_correlates$cor2dose_positive) <- c("Negative","Positive")
# Primary
# Indicator variable
joint_correlates$cor2dose_primary_ind <- factor(joint_correlates$cor2dose_primary, levels = c(0,1))
# And text variable
joint_correlates$cor2dose_primary <- factor(joint_correlates$cor2dose_primary, levels = c(0,1))
levels(joint_correlates$cor2dose_primary) <- c("Non-case","Case")
# Healthcare worker
joint_correlates$cor2dose_hcw_status <- factor(joint_correlates$cor2dose_hcw_status, levels = c("Not_hcw","Less_than_1_covid","More_than_1_covid"))

###
# Define geographic regions
unique(joint_correlates$redcap_data_access_group)
joint_correlates$geographic_region <- NA
# South East
joint_correlates$geographic_region[which(is.element(joint_correlates$redcap_data_access_group,c("oxford__ccvtm","southampton__wtcrf")))] <- "SouthEast"
# London
joint_correlates$geographic_region[which(is.element(joint_correlates$redcap_data_access_group,c("london__hammersmit","london__guys","northwick_park","london__st_georges","london__uclh")))] <- "London"
# West Midlands
joint_correlates$geographic_region[which(is.element(joint_correlates$redcap_data_access_group,c("birmingham__wtcrf","")))] <- "WestMidlands"
# Yorkshire
joint_correlates$geographic_region[which(is.element(joint_correlates$redcap_data_access_group,c("sheffield","")))] <- "Yorkshire"
# Wales
joint_correlates$geographic_region[which(is.element(joint_correlates$redcap_data_access_group,c("cardiff","")))] <- "Wales"
# South West
joint_correlates$geographic_region[which(is.element(joint_correlates$redcap_data_access_group,c("university_hospita","north_bristol")))] <- "SouthWest"
# East Midlands
joint_correlates$geographic_region[which(is.element(joint_correlates$redcap_data_access_group,c("nottingham","")))] <- "EastMidlands"
# Scotland
joint_correlates$geographic_region[which(is.element(joint_correlates$redcap_data_access_group,c("glasgow","edinburgh")))] <- "Scotland"
# North West
joint_correlates$geographic_region[which(is.element(joint_correlates$redcap_data_access_group,c("liverpool")))] <- "NorthWest"
# North East
joint_correlates$geographic_region[which(is.element(joint_correlates$redcap_data_access_group,c("newcastle","hull__castle_hill")))] <- "NorthEast"
# East
joint_correlates$geographic_region[which(is.element(joint_correlates$redcap_data_access_group,c("cambridge")))] <- "East"
table(joint_correlates[,c("geographic_region","cor2dose_positive")])
table(joint_correlates[,c("redcap_data_access_group","cor2dose_positive")])
joint_correlates$geographic_region <- factor(joint_correlates$geographic_region,levels=c("London","SouthEast","SouthWest","East","EastMidlands","WestMidlands","Wales","Yorkshire","NorthWest","NorthEast","Scotland"))

# Some people would find this so offensive! :P
joint_correlates$geographic_big_region <- NA
joint_correlates$geographic_big_region[which(is.element(joint_correlates$geographic_region,c("NorthEast","NorthWest","Scotland")))] <- "NorthEnglandandScotland"
joint_correlates$geographic_big_region[which(is.element(joint_correlates$geographic_region,c("WestMidlands","EastMidlands","Yorkshire")))] <- "MidlandsandYorkshireEngland"
joint_correlates$geographic_big_region[which(is.element(joint_correlates$geographic_region,c("SouthEast","East","Wales","SouthWest")))] <- "SouthEnglandandWales"
joint_correlates$geographic_big_region[which(is.element(joint_correlates$geographic_region,c("London")))] <- "London"
joint_correlates$geographic_big_region <- factor(joint_correlates$geographic_big_region,levels=c("London","SouthEnglandandWales","MidlandsandYorkshireEngland","NorthEnglandandScotland"))
table(joint_correlates$geographic_big_region,useNA="ifany")
table(joint_correlates[,c("geographic_big_region","cor2dose_positive")])
# All have at least 1500 individuals and at least 100 cases.

# Create blank columns for antibody values - to make a wide dataset which we will later convert to long
joint_correlates[,paste0(rep(c("","S.PPD.","RBD.PPD.","S.Jenner.","neuts.","lneuts."),
                             each=length(timepoint_names)),timepoint_names)] <- NA
joint_correlates[,timepoint_names] <- rep(timepoint_names,each=nrow(joint_correlates))
long_correlates <- reshape(joint_correlates,varying = list(S.PPD = paste0("S.PPD.",timepoint_names),
                                                           RBD.PPD = paste0("RBD.PPD.",timepoint_names),
                                                           S.Jenner = paste0("S.Jenner.",timepoint_names),
                                                           neuts = paste0("neuts.",timepoint_names),
                                                           lneuts = paste0("lneuts.",timepoint_names),
                                                           antibody_date = antibody_times,
                                                           timepoint = timepoint_names),
                           direction = "long", v.names = c("S.PPD","RBD.PPD","S.Jenner","neuts","lneuts","antibody_date","timepoint"), idvar = "sc_repeat_pid", times = timepoint_names)
long_correlates <- long_correlates[,setdiff(colnames(long_correlates),"time")]

###############################################################################
# Matching the antibody values

# Only keep certain rows of PPD
PPD <- PPD[which(is.element(PPD$VISIT,timepoint_names)),]
PPD <- PPD[which(is.element(PPD$MBTEST,c("VAC67_ECL_CoV-2_RBD", "VAC67_ECL_CoV-2_S"))),]
PPD <- PPD[which(!is.na(PPD$MBORRES)),]

PPD <- PPD[,c("MBTEST","MBORRES","VISIT","MBDTC","sc_repeat_pid")]
colnames(PPD) <- c("MBTEST","Result","VISIT","MBDTC","sc_repeat_pid")
#joint_correlates[which(joint_correlates$sc_repeat_pid=="20100005"),]

# First create a record_id variable and match onto this to improve the proportion who will be matched
PPD$record_id <- NA
PPD_big_match <- match(PPD$sc_repeat_pid,
                       all_ids$previous_participant_id,
                       incomparables = NA)
PPD$sc_repeat_pid[which(!is.na(PPD_big_match))] <- all_ids$sc_repeat_pid[PPD_big_match[which(!is.na(PPD_big_match))]]
PPD$record_id[which(!is.na(PPD_big_match))] <- all_ids$record_id[PPD_big_match[which(!is.na(PPD_big_match))]]
PPD_big_match <- match(PPD$sc_repeat_pid,
                       all_ids$sc_repeat_pid,
                       incomparables = NA)
PPD$record_id[which(!is.na(PPD_big_match))] <- all_ids$record_id[PPD_big_match[which(!is.na(PPD_big_match))]]

##########

# Then match the S results
PPD_S <- PPD[which(PPD$MBTEST=="VAC67_ECL_CoV-2_S"),]

correlates_rid_date <- paste0(long_correlates$record_id,long_correlates$timepoint)
PPD_S_rid_date <- paste0(PPD_S$record_id,PPD_S$VISIT)

COV002_PPD_match <- match(correlates_rid_date,
                          PPD_S_rid_date,
                          incomparables = NA)

long_correlates$S.PPD[which(!is.na(COV002_PPD_match))] <- PPD_S$Result[COV002_PPD_match[which(!is.na(COV002_PPD_match))]]

# Then match the RBD results
PPD_RBD <- PPD[which(PPD$MBTEST=="VAC67_ECL_CoV-2_RBD"),]
sum(is.na(long_correlates$timepoint))
sum(is.na(long_correlates$record_id))

correlates_rid_date <- paste0(long_correlates$record_id,long_correlates$timepoint)
PPD_RBD_rid_date <- paste0(PPD_RBD$record_id,PPD_RBD$VISIT)

COV002_PPD_match <- match(correlates_rid_date,
                          PPD_RBD_rid_date,
                          incomparables = NA)
long_correlates$RBD.PPD[which(!is.na(COV002_PPD_match))] <- PPD_RBD$Result[COV002_PPD_match[which(!is.na(COV002_PPD_match))]]

#########
# Jenner
table(Jenner$Visit.Day)
unique(substring(Jenner$Vol.no,1,6))
Jenner$sc_repeat_pid <- substring(Jenner$Vol.no,6,nchar(Jenner$Vol.no))
head(Jenner)
# First create a record_id variable and match onto this to improve the proportion who will be matched
Jenner$record_id <- NA
Jenner_big_match <- match(Jenner$sc_repeat_pid,
                          all_ids$previous_participant_id,
                          incomparables = NA)
Jenner$sc_repeat_pid[which(!is.na(Jenner_big_match))] <- all_ids$sc_repeat_pid[Jenner_big_match[which(!is.na(Jenner_big_match))]]
Jenner$record_id[which(!is.na(Jenner_big_match))] <- all_ids$record_id[Jenner_big_match[which(!is.na(Jenner_big_match))]]
Jenner_big_match <- match(Jenner$sc_repeat_pid,
                          all_ids$sc_repeat_pid,
                          incomparables = NA)
Jenner$record_id[which(!is.na(Jenner_big_match))] <- all_ids$record_id[Jenner_big_match[which(!is.na(Jenner_big_match))]]

Jenner$Visit.Day[which(Jenner$Visit.Day=="B+28")] <- "PB+28"
Jenner$Visit.Day[which(Jenner$Visit.Day=="B+90")] <- "PB+90"
Jenner$Visit.Day[which(Jenner$Visit.Day=="B+182")] <- "PB+182"
Jenner$Visit.Day[which(Jenner$Visit.Day=="B+364")] <- "PB+364"

# Match to correlates
correlates_rid_date <- paste0(long_correlates$record_id,long_correlates$timepoint)
Jenner_rid_date <- paste0(Jenner$record_id,Jenner$Visit.Day)

COV002_PPD_match <- match(correlates_rid_date,
                          Jenner_rid_date,
                          incomparables = NA)

long_correlates$S.Jenner[which(!is.na(COV002_PPD_match))] <- Jenner$Std.Spike.ELISA.EU[COV002_PPD_match[which(!is.na(COV002_PPD_match))]]


##########
# Monogram
head(Monogram)
table(Monogram$timepoint)
# Keep COV002 only
Monogram$sc_repeat_pid <- Monogram$record_id
unique(substring(Monogram$sc_repeat_pid,1,1))
Monogram[which(is.na(substring(Monogram$sc_repeat_pid,1,1))),]
Monogram <- Monogram[which(substring(Monogram$SampleID,1,6)=="COV-02" | (substring(Monogram$sc_repeat_pid,1,1)=="2") & nchar(Monogram$sc_repeat_pid)==8),]
Monogram$sc_repeat_pid2 <- substring(Monogram$SampleID,6,13)
table(Monogram$sc_repeat_pid2==Monogram$sc_repeat_pid)
# First create a record_id variable and match onto this to improve the proportion who will be matched
Monogram$record_id <- NA
Monogram_big_match <- match(Monogram$sc_repeat_pid,
                            all_ids$previous_participant_id,
                            incomparables = NA)
Monogram$sc_repeat_pid[which(!is.na(Monogram_big_match))] <- all_ids$sc_repeat_pid[Monogram_big_match[which(!is.na(Monogram_big_match))]]
Monogram$record_id[which(!is.na(Monogram_big_match))] <- all_ids$record_id[Monogram_big_match[which(!is.na(Monogram_big_match))]]
Monogram_big_match <- match(Monogram$sc_repeat_pid,
                            all_ids$sc_repeat_pid,
                            incomparables = NA)
Monogram$record_id[which(!is.na(Monogram_big_match))] <- all_ids$record_id[Monogram_big_match[which(!is.na(Monogram_big_match))]]
table(is.na(Monogram$record_id))
Monogram[which(is.na(Monogram$record_id)),]

# Fix timepoint names
Monogram$timepoint[which(Monogram$timepoint=="PB28")] <- "PB+28"
Monogram$timepoint[which(Monogram$timepoint=="PB90")] <- "PB+90"
Monogram$timepoint[which(Monogram$timepoint=="PB182")] <- "PB+182"

# Match to correlates
correlates_rid_date <- paste0(long_correlates$record_id,long_correlates$timepoint)
Monogram_rid_date <- paste0(Monogram$record_id,Monogram$timepoint)

COV002_PPD_match <- match(correlates_rid_date,
                          Monogram_rid_date,
                          incomparables = NA)

long_correlates$neuts[which(!is.na(COV002_PPD_match))] <- Monogram$SARS_CoV_2[COV002_PPD_match[which(!is.na(COV002_PPD_match))]]
# May be better to use ic50 here instead

##########
# PHE
head(PHE)
table(PHE$timepoint)
PHE$sc_repeat_pid <- PHE$record_id
table(substring(PHE$sc_repeat_pid,1,1))
PHE <- PHE[which(substring(PHE$sc_repeat_pid,1,1)=="2"),]
head(PHE)
# First create a record_id variable and match onto this to improve the proportion who will be matched
PHE$record_id <- NA
PHE_big_match <- match(PHE$sc_repeat_pid,
                       all_ids$previous_participant_id,
                       incomparables = NA)
PHE$sc_repeat_pid[which(!is.na(PHE_big_match))] <- all_ids$sc_repeat_pid[PHE_big_match[which(!is.na(PHE_big_match))]]
PHE$record_id[which(!is.na(PHE_big_match))] <- all_ids$record_id[PHE_big_match[which(!is.na(PHE_big_match))]]
PHE_big_match <- match(PHE$sc_repeat_pid,
                       all_ids$sc_repeat_pid,
                       incomparables = NA)
PHE$record_id[which(!is.na(PHE_big_match))] <- all_ids$record_id[PHE_big_match[which(!is.na(PHE_big_match))]]

# Match to correlates
correlates_rid_date <- paste0(long_correlates$record_id,long_correlates$timepoint)
PHE_rid_date <- paste0(PHE$record_id,PHE$timepoint)

COV002_PPD_match <- match(correlates_rid_date,
                          PHE_rid_date,
                          incomparables = NA)
rm(correlates_rid_date)
long_correlates$lneuts[which(!is.na(COV002_PPD_match))] <- PHE$sample_nd50[COV002_PPD_match[which(!is.na(COV002_PPD_match))]]
# Check sample_nd50 is the right thing to use

# Convert all assays to numeric
# Note I should come back and check the issue here - n/a missing in the lab data may be different to not having the lab data and informative. Ignoring could introduce bias
long_correlates$S.PPD[which(long_correlates$S.PPD=="<33")] <- 33/2
long_correlates$S.PPD <- as.numeric(long_correlates$S.PPD)

long_correlates$RBD.PPD[which(long_correlates$RBD.PPD=="<204")] <- 204/2
long_correlates$RBD.PPD <- as.numeric(long_correlates$RBD.PPD)

# long_correlates$S.Jenner[which(long_correlates$S.Jenner=="<33")] <- 33/2
long_correlates$S.Jenner <- as.numeric(long_correlates$S.Jenner)

long_correlates$neuts[which(long_correlates$neuts=="<40")] <- 40/2
long_correlates$neuts <- as.numeric(long_correlates$neuts)

long_correlates$lneuts[which(long_correlates$lneuts=="<LLOQ")] <- 4.3 # Check this one
long_correlates$lneuts <- as.numeric(long_correlates$lneuts) # Check how to deal with >ULOQ

# Create a visit2 variable
long_correlates$visit2 <- long_correlates$timepoint
long_correlates$visit2[which(is.element(long_correlates$visit2,c("D35","D42","D56","D90")))] <- "PB+28"
long_correlates$visit2 <- factor(long_correlates$visit2, levels = c("PB+28","PB+90","PB+182"))
sum(difftime(eff_excl_COV002$d28_pv_datetime_pb_pb28,eff_excl_COV002$vacc_date_2,units="days") < 20, na.rm=T)
head(difftime(eff_excl_COV002$vacc_date_2,eff_excl_COV002$d28_pv_datetime_pb_pb28,units="days"))
#########
# Summaries
table(PHE$timepoint)
table(Monogram$timepoint)
table(Jenner$Visit.Day)

# Create the variables in long_correlates too
long_covariates <- c("age_group","sc_gender","cor2dose_non_white","cor2dose_comorbidities","cor2dose_bmi_geq_30","cor2dose_interval","cor2dose_primeschedule","cor2dose_positive")
# Check there are sufficient longitudinal observations for each covariate
for(i in 1:length(long_covariates)){
  print(table(joint_correlates[which(joint_correlates$cor2dose_positive==1),c(long_covariates[i],"As_vaccinated_arm_2")]))
}
match_long_joint <- match(long_correlates$sc_repeat_pid,joint_correlates$sc_repeat_pid)
long_correlates[,long_covariates] <- joint_correlates[match_long_joint,long_covariates]
# Include only some columns
explanatory <- c(long_covariates,"cor2dose_hcw_status","cor2dose_primary")
correlates_cols <- c("sc_repeat_pid", "record_id","previous_participant_id","redcap_data_access_group",
                     "sc_age", "sc_gender", "vacc_date", "vacc_date_2", "cor2dose_start_date","cor2dose_end_date",
                     "As_vaccinated_arm_2","cor2dose_outcome","cor2dose_primary","cor2dose_positive","cor2dose_cohort",
                     "study_group","cor2dose_bmi_geq_30","cor2dose_comorbidities","cor2dose_hcw_status",
                     "cor2dose_interval","cor2dose_non_white","cor2dose_schedule","age_group",
                     "start_time","end_time","sc_eth","cor2dose_interval_cont","redcap_data_access_group",
                     "cor2dose_primeschedule","cor2dose_interval_comb","cor2dose_positive_ind","cor2dose_primary_ind",
                     "geographic_region","geographic_big_region","sc_bmi_p2")
long_correlates_cols <- unique(c(correlates_cols,"antibody_date","timepoint","visit2","S.PPD","RBD.PPD","S.Jenner","neuts","lneuts"))
long_correlates <- long_correlates[,c(long_correlates_cols)]
# Further adjustments
# Need times to be time since vacc.
long_correlates$antibody_time <- as.numeric(difftime(as.Date(long_correlates$antibody_date),
                                                     as.Date(long_correlates$cor2dose_start_date),
                                                     units = "days"))
sum(is.na(long_correlates$antibody_time))

write.csv(long_correlates[long_correlates$timepoint %in% c("PB+28","PB+90","PB+182"),
                          c("sc_repeat_pid","antibody_time","antibody_date","timepoint")],
          paste0(output_directory,"/visit_times.csv"),
          row.names = FALSE)
# Include only those with dates
long_correlates <- long_correlates[which(!is.na(long_correlates$antibody_date)),]
# Only those whose dates are after their second dose + 14 days
long_correlates <- long_correlates[which(as.Date(long_correlates$antibody_date)+7>as.Date(long_correlates$cor2dose_start_date)),]
# Only those whose dates are on or before their censoring/event date
long_correlates <- long_correlates[which(as.Date(long_correlates$antibody_date)+7<=as.Date(long_correlates$cor2dose_end_date)),]
# joint_correlates <- joint_correlates[,correlates_cols]
joint_correlates <- joint_correlates[which(joint_correlates$cor2dose_cohort),]
# joint_correlates <- joint_correlates[which(joint_correlates$As_vaccinated_arm_2=="ChAdOx1"),]
# Only keep certain columns
long_correlates$ll <- match(long_correlates$sc_repeat_pid,
                            joint_correlates$sc_repeat_pid,
                            incomparables = NA)

#########
# Create long_correlates_S
# Include only S
long_correlates_S <- long_correlates[which(!is.na(long_correlates$S.PPD)),]
# Define log_S
long_correlates_S$log_S <- log(long_correlates_S$S.PPD)

#########
# Create long_correlates_neuts
# Include only S
long_correlates_neuts <- long_correlates[which(!is.na(long_correlates$neuts)),]
# Define log_neuts
long_correlates_neuts$log_neuts <- log(long_correlates_neuts$neuts)

# Create survival dataset

long_correlates_neuts$ll <- match(long_correlates_neuts$sc_repeat_pid,
                                  joint_correlates$sc_repeat_pid,
                                  incomparables = NA)


correlates_cols <- unique(c(correlates_cols,explanatory))

joint_correlates <- joint_correlates[,c(correlates_cols,antibody_times)]
# system("git rev-parse HEAD", intern=TRUE) gives the name of the latest commit from this GitHub project
# So we name the file by the latest commit.
write.csv(long_correlates_S,
          paste0(data_directory,"/long_correlates_S_0c126098abccf39b82d0d8285b8659cd10c5b355.csv"),
          row.names = FALSE)
write.csv(long_correlates_neuts,
          paste0(data_directory,"/long_correlates_neuts_0c126098abccf39b82d0d8285b8659cd10c5b355.csv"),
          row.names = FALSE)
write.csv(joint_correlates,
          paste0(data_directory,"/joint_correlates_0c126098abccf39b82d0d8285b8659cd10c5b355.csv"),
          row.names = FALSE)