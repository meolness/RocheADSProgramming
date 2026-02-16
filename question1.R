#==============================================================================
# Question 1: SDTM DS Domain Creation
# ==============================================================================
# Purpose: Create ADSL (Subject Level Analysis) dataset from SDTM data
# Author: Aness Chafer
# Date: 2026-02-16
# ==============================================================================



#install.packages("stdm.oak") if needed

#install.packages("pharmaverseraw")

library(sdtm.oak)
library(dplyr)
library(pharmaverseraw)

# Load raw data
ds_raw <- pharmaverseraw::ds_raw

study_ct <- data.frame(
  stringsAsFactors = FALSE,
  codelist_code = rep("C66727", 10),
  term_code = c("C41331","C25250","C28554","C48226","C48227","C48250","C142185","C49628","C49632","C49634"),
  term_value = c("ADVERSE EVENT","COMPLETED","DEATH","LACK OF EFFICACY","LOST TO FOLLOW-UP",
                 "PHYSICIAN DECISION","PROTOCOL VIOLATION","SCREEN FAILURE",
                 "STUDY TERMINATED BY SPONSOR","WITHDRAWAL BY SUBJECT"),
  collected_value = c("Adverse Event","Complete","Dead","Lack of Efficacy","Lost To Follow-Up",
                      "Physician Decision","Protocol Violation","Trial Screen Failure",
                      "Study Terminated By Sponsor","Withdrawal by Subject"),
  term_preferred_term = c("AE","Completed","Died",NA,NA,NA,"Violation",
                          "Failure to Meet Inclusion/Exclusion Criteria",NA,"Dropout"),
  term_synonyms = c("ADVERSE EVENT","COMPLETE","Death",NA,NA,NA,NA,NA,NA,"Discontinued Participation")
)

# Step 1: Start with raw data and create basic SDTM variables
ds <- ds_raw

# STUDYID already exists as STUDY column
# DOMAIN
ds$DOMAIN <- "DS"

# USUBJID: combine STUDY-SITENM-PATNUM  
ds$USUBJID <- paste(ds$STUDY, ds$SITENM, ds$PATNUM, sep = "-")

# DSSEQ: sequence number per subject
ds <- ds %>%
  group_by(USUBJID) %>%
  mutate(DSSEQ = row_number()) %>%
  ungroup()

# Step 2: Map IT.DSTERM to standardized CDISC terms
ds <- left_join(ds, study_ct, by = c("IT.DSTERM" = "collected_value"))

# DSTERM: use mapped value if exists, otherwise keep original
ds$DSTERM <- ifelse(is.na(ds$term_value), ds$IT.DSTERM, ds$term_value)

# DSDECOD: use IT.DSDECOD if available, otherwise use DSTERM
ds$DSDECOD <- ifelse(is.na(ds$IT.DSDECOD) | ds$IT.DSDECOD == "", 
                     ds$DSTERM, 
                     ds$IT.DSDECOD)

# Step 3: Add category
ds$DSCAT <- "DISPOSITION EVENT"

# Step 4: Create visit variables from INSTANCE column
ds$VISIT <- ds$INSTANCE

# Map visit text to visit numbers
ds$VISITNUM <- case_when(
  grepl("Baseline", ds$INSTANCE, ignore.case = TRUE) ~ 1,
  grepl("Week 2", ds$INSTANCE, ignore.case = TRUE) ~ 2,
  grepl("Week 4", ds$INSTANCE, ignore.case = TRUE) ~ 3,
  grepl("Week 8", ds$INSTANCE, ignore.case = TRUE) ~ 4,
  grepl("Week 12", ds$INSTANCE, ignore.case = TRUE) ~ 5,
  grepl("Week 16", ds$INSTANCE, ignore.case = TRUE) ~ 6,
  grepl("Week 20", ds$INSTANCE, ignore.case = TRUE) ~ 7,
  grepl("Week 24", ds$INSTANCE, ignore.case = TRUE) ~ 8,
  grepl("Week 26", ds$INSTANCE, ignore.case = TRUE) ~ 9,
  TRUE ~ NA_real_
)

# Step 5: Add date variables from IT.DSSTDAT
ds$DSDTC <- as.character(ds$IT.DSSTDAT)
ds$DSSTDTC <- as.character(ds$IT.DSSTDAT)
ds$IT.DSSTDAT <- as.Date(trimws(ds$IT.DSSTDAT), format = "%m-%d-%Y")
# Calculate study day relative to first date
ds <- ds %>%
  group_by(USUBJID) %>%
  mutate(
    first_date = min(IT.DSSTDAT, na.rm = TRUE),
    DSSTDY = as.integer(difftime(IT.DSSTDAT, first_date, units = "days")) + 1
  ) %>%
  ungroup() %>%
  select(-first_date, -IT.DSSTDAT)

# Step 6: Select final SDTM variables and sort
ds_final <- ds %>%
  select(STUDYID = STUDY, DOMAIN, USUBJID, DSSEQ, DSTERM, DSDECOD, DSCAT,
         VISITNUM, VISIT, DSDTC, DSSTDTC, DSSTDY) %>%
  arrange(USUBJID, DSSEQ)

# Save outputs
write.csv(ds_final, "question_1_sdtm_ds_domain.csv", row.names = FALSE)
save(ds_final, file = "question_1_sdtm_ds_domain.RData")

print("DS domain created successfully!")
print(paste("Records:", nrow(ds_final)))
print(paste("Subjects:", length(unique(ds_final$USUBJID))))
print(head(ds_final))


# ==============================================================================
#  END OF PROGRAM
# ==============================================================================

