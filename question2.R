#==============================================================================
# Question 2: ADaM ADSL Dataset Creation
# ==============================================================================
# Purpose: Create ADSL (Subject Level Analysis) dataset from SDTM data
# Author: Aness Chafer
# Date: 2026-02-16
# ==============================================================================

# Load required libraries
library(admiral)
library(admiraldev)
library(dplyr)
library(stringr)
library(lubridate)
library(pharmaversesdtm)


# ==============================================================================
# 1. Load Input SDTM Datasets
# ==============================================================================


dm <- pharmaversesdtm::dm
ds <- pharmaversesdtm::ds
ex <- pharmaversesdtm::ex
ae <- pharmaversesdtm::ae
vs <- pharmaversesdtm::vs

dm <- convert_blanks_to_na(dm)
ds <- convert_blanks_to_na(ds)
ex <- convert_blanks_to_na(ex)
ae <- convert_blanks_to_na(ae)
vs <- convert_blanks_to_na(vs)

# ==============================================================================
# 2. Start with DM as Base for ADSL
# ==============================================================================

adsl <- dm %>%
  select(STUDYID, USUBJID, SUBJID, RFSTDTC, RFENDTC, 
         SITEID, AGE, AGEU, SEX, RACE, ETHNIC, 
         ARM, ARMCD, ACTARM, ACTARMCD, COUNTRY, DMDTC, DMDY)


# ==============================================================================
# 3. Derive Age Groups (AGEGR9, AGEGR9N)
# ==============================================================================

adsl <- adsl %>%
  mutate(
    AGEGR9 = case_when(
      AGE < 18 ~ "<18",
      AGE >= 18 & AGE <= 50 ~ "18 - 50",
      AGE > 50 ~ ">50"
    ),
    AGEGR9N = case_when(
      AGE < 18 ~ 1,
      AGE >= 18 & AGE <= 50 ~ 2,
      AGE > 50 ~ 3
    )
  )

# ==============================================================================
# 4. Derive Treatment Start Date/Time (TRTSDTM, TRTSTMF)
# ==============================================================================
# (the first exposure record is in the raw dataset EX).
# derive trt start date exist only if (exdose > 0) or (exdose = 0 and extrt = placebo)

# Filter exposure records for valid doses
ex_filtered <- ex %>%
  filter(
    EXDOSE > 0 |
      (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))
  )

# (derive_vars_dtm() is a function from the admiral package.
# (It converts a character date/time variable (like "2024-01-15T08:30") into a proper R datetime variable.)
ex_valid <- ex_filtered %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST",
    time_imputation = "first",
    flag_imputation = "time",
    ignore_seconds_flag = TRUE
  )

# Derive first exposure date/time for each subject
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_valid,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM),
    mode = "first",
    filter_add = !is.na(EXSTDTM)
  )


# ==============================================================================
# 5. Derive Treatment End Date/Time (TRTEDTM)
# ==============================================================================

# Derive last exposure date/time
ex_end <- ex_filtered %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last",
    flag_imputation = "time",
    ignore_seconds_flag = TRUE
  )


adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_end,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(TRTEDTM = EXENDTM),
    order = exprs(EXENDTM),
    mode = "last",
    filter_add = !is.na(EXENDTM)
  )

#(format date for first and last exposure)
adsl <- adsl %>%
  mutate(
    TRTSDT = if_else(!is.na(TRTSDTM), as.Date(TRTSDTM), as.Date(NA)),
    TRTEDT = if_else(!is.na(TRTEDTM), as.Date(TRTEDTM), as.Date(NA))
  )


# ==============================================================================
# 6. Derive ITT Flag (ITTFL)
# ==============================================================================


# (ITTFL Set to "Y" if [DM.ARM] not equal to missing Else set to "N" )
adsl <- adsl %>%
  mutate(
    ITTFL = if_else(!is.na(ARM), "Y", "N")
  )


# ==============================================================================
# 7. Derive Last Known Alive Date (LSTAVLDT)
# ==============================================================================


# (1) Last vital signs date with valid result
# (1.  last complete date of vital assessment with a valid test result ([VS.VSSTRESN] and [VS.VSSTRESC] not both missing) 
# and datepart of [VS.VSDTC] not missing. )
vs_dates <- vs %>%
  filter(!is.na(VSSTRESN) | !is.na(VSSTRESC)) %>%
  derive_vars_dt(
    dtc = VSDTC,
    new_vars_prefix = "VS"
  ) %>%
  filter(!is.na(VSDT)) %>%
  group_by(STUDYID, USUBJID) %>%
  summarise(
    LSTVS_DT = max(VSDT, na.rm = TRUE),
    .groups = "drop"
  )

# (2) Last AE onset date
# ( (2) last complete onset date of AEs (datepart of Start Date/Time of Adverse Event [AE.AESTDTC]))
ae_dates <- ae %>%
  derive_vars_dt(
    dtc = AESTDTC,
    new_vars_prefix = "AEST"
  ) %>%
  filter(!is.na(AESTDT)) %>%
  group_by(STUDYID, USUBJID) %>%
  summarise(
    LSTAE_DT = max(AESTDT, na.rm = TRUE),
    .groups = "drop"
  )

# (3) Last disposition date
# (3) last complete disposition date (datepart of Start Date/Time of Disposition Event [DS.DSSTDTC]))
ds_dates <- ds %>%
  derive_vars_dt(
    dtc = DSSTDTC,
    new_vars_prefix = "DSST"
  ) %>%
  filter(!is.na(DSSTDT)) %>%
  group_by(STUDYID, USUBJID) %>%
  summarise(
    LSTDS_DT = max(DSSTDT, na.rm = TRUE),
    .groups = "drop"
  )

# (4) -- information already in the dataset ADSL - no filter needed but need to transform as a date
# Extract date from TRTEDTM
adsl <- adsl %>%
  mutate(
    LSTEX_DT = if_else(!is.na(TRTEDTM), as.Date(TRTEDTM), as.Date(NA))
  )

# (5) Merge all dates back to ADSL 
# (to find the overall last dateof exposture by subjid)
adsl <- adsl %>%
  left_join(vs_dates, by = c("STUDYID", "USUBJID")) %>%
  left_join(ae_dates, by = c("STUDYID", "USUBJID")) %>%
  left_join(ds_dates, by = c("STUDYID", "USUBJID"))



# Calculate LSTAVLDT as max of all dates
adsl <- adsl %>%
  rowwise() %>%
  mutate(
    LSTAVLDT = max(c(LSTVS_DT, LSTAE_DT, LSTDS_DT, LSTEX_DT), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    LSTAVLDT = if_else(is.infinite(LSTAVLDT), as.Date(NA), LSTAVLDT)
  )

# Clean up temporary variables
# (on supprime les variables qu'on a utilise mais qui ne doivent pas figurer dans le dataset adsl)
adsl <- adsl %>%
  select(-c(LSTVS_DT, LSTAE_DT, LSTDS_DT, LSTEX_DT))

# ==============================================================================
# 9. Additional Standard ADSL Variables
# ==============================================================================

# Safety population flag
adsl <- adsl %>%
  mutate(
    SAFFL = if_else(!is.na(TRTSDT), "Y", "N")
  )

# Completed flag - based on disposition
ds_comp <- ds %>%
  filter(DSDECOD == "COMPLETED") %>%
  select(STUDYID, USUBJID) %>%
  distinct() %>%
  mutate(COMPFL = "Y")

adsl <- adsl %>%
  select(-COMPFL) %>%  # remove old one if it exists
  left_join(ds_comp, by = c("STUDYID", "USUBJID")) %>%
  mutate(COMPFL = if_else(is.na(COMPFL), "N", "Y"))

# Treatment duration
adsl <- adsl %>%
  mutate(
    TRTDURD = if_else(
      !is.na(TRTSDT) & !is.na(TRTEDT),
      as.numeric(difftime(TRTEDT, TRTSDT, units = "days")) + 1,
      NA_real_
    )
  )

# ==============================================================================
# 10. Select and Order Final Variables
# ==============================================================================

adsl_final <- adsl %>%
  select(
    # Identifiers
    STUDYID, USUBJID, SUBJID, SITEID,
    # Treatment variables
    ARM, ARMCD, ACTARM, ACTARMCD,
    # Demographics
    AGE, AGEU, AGEGR9, AGEGR9N, SEX, RACE, ETHNIC, COUNTRY,
    # Treatment dates and times
    TRTSDT, TRTEDT, TRTSDTM, TRTEDTM, TRTSTMF, TRTDURD,
    # Reference dates
    RFSTDTC, RFENDTC,
    # Derived dates
    LSTAVLDT,
    # Population flags
    ITTFL, SAFFL, COMPFL,
    # Original DM dates
    DMDTC, DMDY
  ) %>%
  arrange(USUBJID)


# ==============================================================================
#  END OF PROGRAM
# ==============================================================================

