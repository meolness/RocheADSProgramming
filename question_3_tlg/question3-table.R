# ==============================================================================
# Question 3 : Adverse Events Summary Table using {gtsummary}
# ==============================================================================
# Purpose: Create summary table of treatment-emergent adverse events 
# Author: Aness Chafer
# Date: 2026-02-16
# ==============================================================================

# Libraries
library(pharmaverseadam)
library(gtsummary)
library(dplyr)
library(tidyr)
library(gt)

# ==============================================================================
# 1. Load Input Datasets
# ==============================================================================

adae <- pharmaverseadam::adae
adsl <- pharmaverseadam::adsl

cat("ADAE:", nrow(adae), "records\n")
cat("ADSL:", nrow(adsl), "subjects\n\n")

# ==============================================================================
# 2. Filter Treatment-Emergent Adverse Events
# ==============================================================================

#Filtering treatment-emergent adverse events (TRTEMFL = 'Y')...\n

adae_teae <- adae %>%
  filter(TRTEMFL == "Y")

cat("Treatment-emergent AEs:", nrow(adae_teae), "records\n\n")

# ==============================================================================
# 3. Prepare Data for Summary Table by Preferred Term
# ==============================================================================


# Get denominators (total subjects per treatment arm)
denom <- adsl %>%
  group_by(ACTARM) %>%
  summarise(N = n(), .groups = "drop")

cat("Treatment arm denominators:\n")
print(denom)

# Count subjects with at least one AE by term and treatment
ae_by_term <- adae_teae %>%
  group_by(ACTARM, AETERM) %>%
  summarise(
    n = n_distinct(USUBJID),
    .groups = "drop"
  ) %>%
  left_join(denom, by = "ACTARM") %>%
  mutate(
    pct = (n / N) * 100,
    display = sprintf("%d (%.1f%%)", n, pct)
  )

# Pivot to wide format for table
ae_table_wide <- ae_by_term %>%
  select(AETERM, ACTARM, display) %>%
  pivot_wider(
    names_from = ACTARM,
    values_from = display,
    values_fill = "0 (0.0%)"
  )

# Calculate total column
ae_total <- adae_teae %>%
  group_by(AETERM) %>%
  summarise(
    n = n_distinct(USUBJID),
    .groups = "drop"
  ) %>%
  mutate(
    N = sum(denom$N),
    pct = (n / N) * 100,
    Total = sprintf("%d (%.1f%%)", n, pct)
  ) %>%
  select(AETERM, Total)

ae_table_wide <- ae_table_wide %>%
  left_join(ae_total, by = "AETERM")

# Sort by descending total frequency
ae_table_wide <- ae_table_wide %>%
  mutate(
    sort_n = as.numeric(sub(" .*", "", Total))
  ) %>%
  arrange(desc(sort_n)) %>%
  select(-sort_n)

cat(nrow(ae_table_wide), "unique preferred terms\n\n")

# ==============================================================================
# 4. Create Enhanced Table using {gt}
# ==============================================================================


treat_arms <- setdiff(names(ae_table_wide), c("AETERM", "Total"))

print("Treatment arms for table:")
print(treat_arms)

# Create gt table
ae_summary_table <- ae_table_wide %>%
  gt() %>%
  tab_header(
    title = md("**Summary of Treatment-Emergent Adverse Events**"),
    subtitle = md("*Number (Percentage) of Subjects with TEAEs by Preferred Term*")
  ) %>%
  cols_label(AETERM = md("**Preferred Term**")) %>%
  tab_spanner(label = md("**Treatment Group**"), columns = all_of(treat_arms)) %>%
  tab_source_note(source_note = "TEAE = Treatment-Emergent Adverse Event (TRTEMFL = 'Y')") %>%
  tab_source_note(source_note = paste0("Analysis based on ", nrow(adsl), " subjects")) %>%
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.border.bottom.color = "black",
    column_labels.border.bottom.color = "black",
    column_labels.font.weight = "bold"
  ) %>%
  cols_align(align = "left", columns = AETERM) %>%
  cols_align(align = "center", columns = -AETERM)


# Create alternative summary by System Organ Class
ae_by_soc <- adae_teae %>%
  group_by(ACTARM, AESOC) %>%
  summarise(n = n_distinct(USUBJID), .groups = "drop") %>%
  left_join(denom, by = "ACTARM") %>%
  mutate(pct = (n / N) * 100,
         display = sprintf("%d (%.1f%%)", n, pct))

ae_soc_wide <- ae_by_soc %>%
  select(AESOC, ACTARM, display) %>%
  pivot_wider(names_from = ACTARM, values_from = display, values_fill = "0 (0.0%)")

# Calculate total for SOC
ae_soc_total <- adae_teae %>%
  group_by(AESOC) %>%
  summarise(n = n_distinct(USUBJID), .groups = "drop") %>%
  mutate(N = sum(denom$N),
         pct = (n / N) * 100,
         Total = sprintf("%d (%.1f%%)", n, pct)) %>%
  select(AESOC, Total)

ae_soc_wide <- ae_soc_wide %>%
  left_join(ae_soc_total, by = "AESOC") %>%
  mutate(sort_n = as.numeric(sub(" .*", "", Total))) %>%
  arrange(desc(sort_n)) %>%
  select(-sort_n)

# Create SOC table
ae_soc_table <- ae_soc_wide %>%
  gt() %>%
  tab_header(
    title = md("**Summary of Treatment-Emergent Adverse Events**"),
    subtitle = md("*Number (Percentage) of Subjects with TEAEs by System Organ Class*")
  ) %>%
  cols_label(AESOC = md("**System Organ Class**")) %>%
  tab_spanner(label = md("**Treatment Group**"), columns = all_of(treat_arms)) %>%
  tab_source_note(source_note = "TEAE = Treatment-Emergent Adverse Event (TRTEMFL = 'Y')") %>%
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    column_labels.border.bottom.color = "black",
    column_labels.font.weight = "bold"
  ) %>%
  cols_align(align = "left", columns = AESOC) %>%
  cols_align(align = "center", columns = -AESOC)

#View(ae_soc_table)

# Save outputs
gtsave(ae_summary_table, "question_3_tlg/ae_summary_table.html")
print("✓ HTML table saved: question_3_tlg/ae_summary_table.html")

gtsave(ae_soc_table, "question_3_tlg/ae_summary_table_by_soc.html")
print("✓ HTML SOC table saved: question_3_tlg/ae_summary_table_by_soc.html")

gtsave(ae_summary_table, "question_3_tlg/ae_summary_table.docx")
print("✓ DOCX table saved: question_3_tlg/ae_summary_table.docx")

write.csv(ae_table_wide, "question_3_tlg/ae_summary_data.csv", row.names = FALSE)
print("✓ Data CSV saved: question_3_tlg/ae_summary_data.csv")