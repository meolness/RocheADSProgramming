# ==============================================================================
# Question 3 Part 2: Adverse Events Visualizations using {ggplot2}
# ==============================================================================
# Purpose: Create visualizations for adverse events analysis
# Author: Aness CHAFER
# Date: 2026-02-16
# ==============================================================================

# Libraries
library(pharmaverseadam)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(forcats)

# ==============================================================================
# 1. Load Input Datasets
# ==============================================================================

adae <- pharmaverseadam::adae
adsl <- pharmaverseadam::adsl

# Filter treatment-emergent AEs
adae_teae <- adae %>%
  filter(TRTEMFL == "Y")

cat("ADAE TEAE records:", nrow(adae_teae), "\n")
cat("ADSL subjects:", nrow(adsl), "\n\n")

# ==============================================================================
# 2. PLOT 1: AE Severity Distribution by Treatment
# ==============================================================================

# Prepare data for severity plot
severity_data <- adae_teae %>%
  filter(!is.na(AESEV) & !is.na(ACTARM)) %>%
  group_by(ACTARM, AESEV) %>%
  summarise(
    n_events = n(),
    n_subjects = n_distinct(USUBJID),
    .groups = "drop"
  )

# Get denominators
denom <- adsl %>%
  group_by(ACTARM) %>%
  summarise(N = n(), .groups = "drop")

severity_data <- severity_data %>%
  left_join(denom, by = "ACTARM") %>%
  mutate(
    pct_subjects = (n_subjects / N) * 100
  )

# Severity levels
severity_data <- severity_data %>%
  mutate(
    AESEV = factor(AESEV, levels = c("MILD", "MODERATE", "SEVERE"))
  )

# Create bar chart
plot1_bar <- ggplot(severity_data, aes(x = ACTARM, y = pct_subjects, fill = AESEV)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_subjects)),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 3
  ) +
  scale_fill_manual(
    values = c("MILD" = "#90EE90", "MODERATE" = "#FFA500", "SEVERE" = "#FF6B6B"),
    name = "Severity"
  ) +
  labs(
    title = "Adverse Event Severity Distribution by Treatment Group",
    subtitle = "Percentage of Subjects with Treatment-Emergent Adverse Events",
    x = "Treatment Arm",
    y = "Percentage of Subjects (%)",
    caption = "Based on treatment-emergent adverse events (TRTEMFL = 'Y')"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# Save bar chart
ggsave(
  "question_3_tlg/plot1_ae_severity_bar.png",
  plot = plot1_bar,
  width = 10,
  height = 7,
  dpi = 300
)

# Create alternative heatmap visualization
severity_heatmap <- severity_data %>%
  select(ACTARM, AESEV, pct_subjects)

plot1_heatmap <- ggplot(severity_heatmap, aes(x = ACTARM, y = AESEV, fill = pct_subjects)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_subjects)),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_gradient(
    low = "#FFFFCC",
    high = "#FF4444",
    name = "% of\nSubjects"
  ) +
  labs(
    title = "Adverse Event Severity Heatmap by Treatment Group",
    subtitle = "Percentage of Subjects with Treatment-Emergent Adverse Events",
    x = "Treatment Arm",
    y = "Severity Grade",
    caption = "Based on treatment-emergent adverse events (TRTEMFL = 'Y')"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

# Save heatmap
ggsave(
  "question_3_tlg/plot1_ae_severity_heatmap.png",
  plot = plot1_heatmap,
  width = 10,
  height = 6,
  dpi = 300
)


# ==============================================================================
# 4. PLOT BONUS: AE Onset Timeline
# ==============================================================================
# Prepare data for AE onset over time
ae_timeline <- adae_teae %>%
  filter(!is.na(ASTDY)) %>%
  mutate(
    week = ceiling(ASTDY / 7)
  ) %>%
  filter(week > 0 & week <= 52) %>%  # Limit to 1 year
  group_by(week, ACTARM) %>%
  summarise(
    n_events = n(),
    .groups = "drop"
  )

if(nrow(ae_timeline) > 0) {
  plot_bonus <- ggplot(ae_timeline, aes(x = week, y = n_events, color = ACTARM, group = ACTARM)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_brewer(palette = "Set1", name = "Treatment Arm") +
    labs(
      title = "Adverse Event Onset Over Study Time",
      subtitle = "Number of Treatment-Emergent AEs by Study Week",
      x = "Study Week",
      y = "Number of Adverse Events",
      caption = "Based on treatment-emergent adverse events (TRTEMFL = 'Y')"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      axis.title = element_text(size = 11, face = "bold")
    ) +
    scale_x_continuous(breaks = seq(0, 52, by = 4))
}

ggsave(
  "question_3_tlg/plot_bonus_ae_timeline.png",
  plot = plot_bonus,
  width = 12,
  height = 7,
  dpi = 300
)

# ==============================================================================
#  END OF PROGRAM
# ==============================================================================

