library(tidyverse)
library(readxl)
library(ggtext)


# *** Physiological phenotyping ***

# ENV - soil water content
tmp.file <- paste(DATA_ROOT, "Experiments", "Phenotyping.soil.tsv", sep = '/')
phenotyping.env.water <- read_tsv(tmp.file, col_names = read_tsv(tmp.file) %>% names(), col_select = Label:'Theta(9)', skip = 2, locale = locale(decimal_mark = ",")) %>%
  mutate(DateTime =  dmy_hms(Label), Label = NULL) %>%
  pivot_longer(cols = `Theta(1)`:`Theta(9)`, names_to = "Sensor", values_to = "WaterContent") %>%
  mutate(
    Species = case_when(
      Sensor %in% c("Theta(1)", "Theta(2)", "Theta(3)") ~ "Clusia rosea",
      Sensor %in% c("Theta(4)", "Theta(5)", "Theta(6)") ~ "Clusia minor s.l.",
      Sensor %in% c("Theta(7)", "Theta(8)", "Theta(9)") ~ "Clusia major"),
    Individual = case_when(
      Sensor %in% c("Theta(1)", "Theta(4)", "Theta(7)") ~ "A",
      Sensor %in% c("Theta(2)", "Theta(5)", "Theta(8)") ~ "B",
      Sensor %in% c("Theta(3)", "Theta(6)", "Theta(9)") ~ "C"),
  .before = Sensor) %>%
  mutate(
    Species = fct_relevel(Species, c("Clusia major", "Clusia minor s.l.", "Clusia rosea")),
    Individual = fct_relevel(Individual, c("A", "B", "C"))) %>%
  relocate(Power, .after = Sensor) %>%
  arrange(DateTime, Species, Individual)

pdf(paste(RESULTS_DIR, "figs", "phenotyping.env.water.pdf", sep = '/'), width = 7, height = 3)
phenotyping.env.water %>%
  ggplot(aes(x = DateTime, y = WaterContent, color = Species)) +
  geom_vline(xintercept = ymd_hms("2025-11-17 21:00:00"), linetype="dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = ymd_hms("2025-11-18 21:00:00"), linetype="dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = ymd_hms("2025-12-03 21:00:00"), linetype="dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = ymd_hms("2025-12-04 21:00:00"), linetype="dashed", color = "black", alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, span = 0.1) +
  scale_color_manual(values = c("#A73130", "grey50", "#133F66")) +
  scale_x_datetime(name = "Date", limits = c(ymd("2025-11-01", tz = "UTC"), ymd("2025-12-05", tz = "UTC"))) +
  scale_y_continuous(name = "**Soil water volume**<span style='font-size:12px'> [%]</span>", limits = c(NA, NA)) +
  coord_cartesian(ylim=c(8, 32)) +
  theme_bw() +
  theme(axis.title.y = element_markdown(), axis.title.x = element_text(face = "bold"), axis.line = element_line(), axis.ticks = element_line())
dev.off()

# ENV - light scheme (PAR)
phenotyping.env.par <- read_csv(paste(DATA_ROOT, "Experiments", "Phenotyping.light.csv", sep = '/')) %>%
  mutate(DateTime = ymd(Date) + hms(Timestamp), .before = Date)

pdf(paste(RESULTS_DIR, "figs", "phenotyping.env.light.pdf", sep = '/'), width = 5, height = 1.5)
phenotyping.env.par %>%
  ggplot(aes(x = DateTime, y = μmoles)) +
  annotate("rect", xmin = ymd_hms("2025-11-15 21:00:00"), xmax = ymd_hms("2025-11-16 06:00:00"), ymin = -Inf, ymax = Inf,  fill = "grey30", alpha = 0.1) +
  geom_hline(yintercept=0, linetype="solid", color = "gray") +
  geom_line(color = "#DBA800") +
  scale_y_continuous(name = "**PAR**") +
  scale_x_datetime(name = "Time", date_labels="%H:%M", limits = c(ymd_hms("2025-11-15 21:00:00"), ymd_hms("2025-11-16 21:10:00"))) +
  theme_bw() +
  theme(axis.title.y = element_markdown(), axis.title.x = element_text(face = "bold"), axis.line = element_line(), axis.ticks = element_line())
dev.off()

# ENV - temperature & rel. humidity
phenotyping.env.fytotron <- read_tsv(paste(DATA_ROOT, "Experiments", "Phenotyping.env.tsv", sep = '/')) %>%
  mutate(DateTime = dmy_hms(Date), .before = Date) %>%
  mutate(T_Actual = T_Actual / 10)

pdf(paste(RESULTS_DIR, "figs", "phenotyping.env.temp.pdf", sep = '/'), width = 5, height = 1.5)
phenotyping.env.fytotron %>%
  ggplot(aes(x = DateTime, y = T_Actual)) +
  annotate("rect", xmin = ymd_hms("2025-11-15 21:00:00"), xmax = ymd_hms("2025-11-16 06:00:00"), ymin = -Inf, ymax = Inf,  fill = "grey30", alpha = 0.1) +
  geom_hline(yintercept=0, linetype="solid", color = "gray") +
  geom_line(color = "#A73130") +
  scale_y_continuous(name = "**Temp.**<span style='font-size:12px'> [°C]</span>", breaks = c(20, 24, 28), limits = c(19, 29)) +
  scale_x_datetime(name = "Time", date_labels="%H:%M", limits = c(ymd_hms("2025-11-15 21:00:00"), ymd_hms("2025-11-16 21:10:00"))) +
  theme_bw() +
  theme(axis.title.y = element_markdown(), axis.title.x = element_text(face = "bold"), axis.line = element_line(), axis.ticks = element_line())
dev.off()

pdf(paste(RESULTS_DIR, "figs", "phenotyping.env.rh.pdf", sep = '/'), width = 5, height = 1.5)
phenotyping.env.fytotron %>%
  ggplot(aes(x = DateTime, y = Rh_Actual)) +
  annotate("rect", xmin = ymd_hms("2025-11-15 21:00:00"), xmax = ymd_hms("2025-11-16 06:00:00"), ymin = -Inf, ymax = Inf,  fill = "grey30", alpha = 0.1) +
  geom_hline(yintercept=0, linetype="solid", color = "gray") +
  geom_line(color = "#133F66") +
  scale_y_continuous(name = "**Humidity**<span style='font-size:12px'> [%]</span>", limits = c(35, 70)) +
  scale_x_datetime(name = "Time", date_labels="%H:%M", limits = c(ymd_hms("2025-11-15 21:00:00"), ymd_hms("2025-11-16 21:10:00"))) +
  theme_bw() +
  theme(axis.title.y = element_markdown(), axis.title.x = element_text(face = "bold"), axis.line = element_line(), axis.ticks = element_line())
dev.off()



# Titratable acidity measurements
phenotyping.acidity <- read_excel(paste(DATA_ROOT, "Experiments", "Phenotyping.acidity.xlsx", sep = '/'), sheet = "Acidity", range = "A1:T91", col_names = TRUE) %>%
  mutate(Group = substr(SampleID, 1, 3), .after = SampleID) %>%
  relocate(Treatment, .after = Group) %>%
  mutate(BioRep = substr(Individual, 2, 2), .after = Individual) %>%
  mutate(DateTime = case_when(
    Timepoint == 1 ~ ymd_hms("2025-11-17 21:00:00"),
    Timepoint == 2 ~ ymd_hms("2025-11-18 06:00:00"),
    Timepoint == 3 ~ ymd_hms("2025-11-18 09:00:00"),
    Timepoint == 4 ~ ymd_hms("2025-11-18 18:00:00"),
    Timepoint == 5 ~ ymd_hms("2025-11-18 21:00:00"),
    .default = NA),
  .after = Timepoint) %>%
  mutate_at(c("Group", "Treatment", "Species", "Individual", "BioRep", "Timepoint"), as.factor) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Well-watered", "Drought"))) %>%
  rename(TA1 = "TA [μmol H+/g FW]") %>%
  rowwise(SampleID) %>%
  mutate(TA.mean = mean(c_across(starts_with("TA"))), TA.sd = sd(c_across(starts_with("TA")))) %>%
  ungroup()


# pH
pdf(paste(RESULTS_DIR, "figs", "phenotyping.ph.pdf", sep = '/'), width = 10, height = 6)
phenotyping.acidity %>%
  ggplot(aes(x = DateTime, y = -pH, color = Treatment, shape = Treatment, group = Group)) +
  annotate("rect", xmin = ymd_hms("2025-11-17 21:00:00"), xmax = ymd_hms("2025-11-18 06:00:00"), ymin = -Inf, ymax = Inf,  fill = "grey30", alpha = 0.1) +
  annotate("rect", xmin = ymd_hms("2025-11-18 21:00:00"), xmax = ymd_hms("2025-11-18 22:00:00"), ymin = -Inf, ymax = Inf,  fill = "grey30", alpha = 0.1) +
  geom_line(na.rm=TRUE) +
  geom_point(na.rm=TRUE) +
  scale_y_continuous(name = "**-pH**", limits = c(NA, NA)) +
  scale_x_datetime(name = "Time", date_labels="%H:%M") +
  scale_color_manual(breaks = c("Well-watered", "Drought"), values = c("#133F66", "#A73130")) +
  facet_grid(Species ~ BioRep, switch = "y") +
  theme_minimal() +
  theme(axis.title.y = element_markdown(), axis.title.x = element_text(face = "bold"), axis.line = element_line(), axis.ticks = element_line()) +
  theme(strip.text = element_text(face = "bold", size = rel(1)), strip.background = element_rect(fill = "gray95", colour = "black", linewidth = 0))

dev.off()


# TA over 24-hours
pdf(paste(RESULTS_DIR, "figs", "phenotyping.ta.24h.pdf", sep = '/'), width = 10, height = 6)
phenotyping.acidity %>%
  ggplot(aes(x = DateTime, y = TA.mean, color = Treatment, shape = Treatment, group = Group)) +
  annotate("rect", xmin = ymd_hms("2025-11-17 21:00:00"), xmax = ymd_hms("2025-11-18 06:00:00"), ymin = -Inf, ymax = Inf,  fill = "grey30", alpha = 0.1) +
  annotate("rect", xmin = ymd_hms("2025-11-18 21:00:00"), xmax = ymd_hms("2025-11-18 22:00:00"), ymin = -Inf, ymax = Inf,  fill = "grey30", alpha = 0.1) +
  geom_line(na.rm=TRUE) +
  geom_errorbar(aes(ymin = TA.mean-TA.sd, ymax = TA.mean+TA.sd), size = 0.2) +
  geom_point(na.rm=TRUE) +
  scale_y_continuous(name = "**Titratable Acidity** <span style='font-size:12px'>[µmol H<sup>+</sup>  g FW<sup>-1</sup>]</span>", limits = c(NA, NA)) +
  scale_x_datetime(name = "Time", date_labels="%H:%M") +
  scale_color_manual(breaks = c("Well-watered", "Drought"), values = c("#133F66", "#A73130")) +
  facet_grid(Species ~ BioRep, switch = "y") +
  theme_minimal() +
  theme(axis.title.y = element_markdown(), axis.title.x = element_text(face = "bold"), axis.line = element_line(), axis.ticks = element_line()) +
  theme(strip.text = element_text(face = "bold", size = rel(1)), strip.background = element_rect(fill = "gray95", colour = "black", linewidth = 0))

dev.off()


# delta TA
phenotyping.acidity.delta <- phenotyping.acidity %>%
  pivot_longer(cols = c(TA1, TA2, TA3), names_to = "TechRep", values_to = "TA", values_drop_na = TRUE) %>%
  pivot_wider(id_cols = c(Treatment, Species, Individual, BioRep, TechRep), names_from = Timepoint, names_prefix = "T", values_from = TA) %>%
  mutate(TA.delta = T2-T1) %>%
  group_by(Treatment, Species, Individual) %>%
  summarise(TA.TechRep.mean = mean(TA.delta)) %>%
  mutate(TA.BioRep.n = n(), TA.BioRep.mean = mean(TA.TechRep.mean), TA.BioRep.sd = sd(TA.TechRep.mean), TA.BioRep.se = TA.BioRep.sd / sqrt(TA.BioRep.n))

# Statistical test via two-way ANOVA and Tukey HSD
phenotyping.acidity.delta.aov <-
  TukeyHSD(
    aov(
      TA.TechRep.mean ~ Species * Treatment,
      data = phenotyping.acidity.delta %>% filter(Species != "Clusia minor")),
    which = "Species:Treatment")$`Species:Treatment` %>%
  as_tibble(rownames = "Species:Treatment")


pdf(paste(RESULTS_DIR, "figs", "phenotyping.ta.delta.main.pdf", sep = '/'), width = 5, height = 7)
phenotyping.acidity.delta %>%
  filter(Species != "Clusia minor") %>%
  ggplot(aes(x = Treatment, y = TA.BioRep.mean, fill = Treatment)) +
    geom_hline(yintercept=0, linetype="longdash", color = "black") +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), color = "black", width = 0.6) +
    geom_errorbar(aes(ymin = TA.BioRep.mean - TA.BioRep.sd, ymax = TA.BioRep.mean + TA.BioRep.sd), width = 0.3, size = 0.3, alpha = 0.7) +
    geom_jitter(aes(y = TA.TechRep.mean), position = position_jitterdodge(jitter.width = 0.1, dodge.width  = 0.7), size = 2, alpha = 0.7, shape = 21, color = "black") +
    scale_fill_manual(values = c("Well-watered" = "steelblue", "Drought" = "firebrick")) +
    labs(x = "Treatment", y = expression(Delta~"TA ["*mu*"mol H"^"+"*" g FW"^{-1}*"] (Morning – Evening)"), fill = "Treatment") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(strip.text = element_text(face = "bold", size = rel(0.7)), strip.background = element_rect(fill = "gray95", colour = "black", linewidth = 0)) +
    facet_wrap(~ Species, scales = "free_x")

dev.off()

pdf(paste(RESULTS_DIR, "figs", "phenotyping.ta.delta.supplement.pdf", sep = '/'), width = 4, height = 7)
phenotyping.acidity.delta %>%
  ggplot(aes(x = Species, y = TA.BioRep.mean, fill = Treatment)) +
    geom_hline(yintercept = 0, linetype = "longdash", color = "black") +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), color = "black", width = 0.6) +
    geom_errorbar(aes(ymin = TA.BioRep.mean - TA.BioRep.sd, ymax = TA.BioRep.mean + TA.BioRep.sd), width = 0.3, size = 0.3, alpha = 0.7) +
    geom_jitter(aes(y = TA.TechRep.mean), position = position_jitterdodge(jitter.width = 0.1, dodge.width  = 0.7), size = 2, alpha = 0.7, shape = 21, color = "black") +
    scale_fill_manual(breaks = c("Well-watered", "Drought"), values = c("steelblue", "#A73130")) +
    labs(x = "Species", y = expression(Delta~"TA ["*mu*"mol H"^"+"*" g FW"^{-1}*"] (Morning – Evening)"), fill = "Treatment") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(strip.text = element_text(face = "bold", size = rel(0.7)), strip.background = element_rect(fill = "gray95", colour = "black", linewidth = 0)) +
    facet_grid(. ~ Treatment)

dev.off()



# Gas exchange measurement
tmp.file <- paste(DATA_ROOT, "Experiments", "Phenotyping.gasexchange.csv", sep = '/')
phenotyping.gasexchange <- read_delim(tmp.file, delim = ';', col_names = read_delim(tmp.file, delim = ';') %>% names(), col_select = Date:wa, skip = 2) %>%
  mutate(DateTime = ymd(Date) + hms(Time), .before = Date) %>%
  mutate_at(c("Code", "Object", "Status"), as.factor) %>%
  filter(Code == "MP_001") %>%
  slice(which(row_number() %% 2 == 1)) %>%
  mutate(A = A-0.7, GH2O = GH2O-2) %>%  # Zero point calibration
  mutate(GH2O = if_else(GH2O < 0, 0, GH2O)) # Issue of negative values, see Methods


# Clusia major (Figure 1b)
start <- ymd_hms("2025-05-15 21:00:00")
end <- ymd_hms("2025-05-16 21:00:00")
pdf(paste(RESULTS_DIR, "figs", "phenotyping.gasexchange.major.pdf", sep = '/'), width = 7, height = 4)

phenotyping.gasexchange %>%
  mutate(x = ymd_hms(cut(DateTime, breaks = "25 mins"), quiet = TRUE)) %>%
  group_by(x) %>%
  summarise(A = mean(A), GH2O = mean(GH2O)) %>%
  ggplot(aes(x = x)) +
  geom_hline(yintercept=0, linetype="longdash", color = "gray") +
  geom_line(aes(y = GH2O/6), na.rm=TRUE, color = "#133F66") +
  geom_point(aes(y = na_if(GH2O, 0)/6), na.rm=TRUE, color = "#133F66") +
  geom_line(aes(y = A), na.rm=TRUE, color = "#A73130") +
  geom_point(aes(y = A), na.rm=TRUE, color = "#A73130") +
  scale_y_continuous(name = "<span style='color:#A73130'>**Net carbon uptake (CO<sub>2</sub>)**</span><br>µM m<sup>-2</sup> s<sup>-1</sup>", limits = c(-1,10), breaks = c(0,2,4,6,8,10), sec.axis = sec_axis(trans=~.*6, name="<span style='color:#133F66'>**Leaf conductance (H<sub>2</sub>O)**</span><br>mM m<sup>-2</sup> s<sup>-1</sup>")) +
  scale_x_datetime(name = "Time", limits = c(start,end + hms("00:10:00")), date_labels="%H:%M") +
  theme_classic() +
  theme(axis.title.y = element_markdown(), legend.position="right")

dev.off()

# Clusia rosea (Figure 1b)
start <- ymd_hms("2025-05-17 21:00:00")
end <- ymd_hms("2025-05-18 21:00:00")
pdf(paste(RESULTS_DIR, "figs", "phenotyping.gasexchange.rosea.pdf", sep = '/'), width = 7, height = 4)

phenotyping.gasexchange %>%
  mutate(x = ymd_hms(cut(DateTime, breaks = "25 mins"), quiet = TRUE)) %>%
  group_by(x) %>%
  summarise(A = mean(A), GH2O = mean(GH2O)) %>%
  ggplot(aes(x = x)) +
  geom_hline(yintercept=0, linetype="longdash", color = "gray") +
  geom_line(aes(y = GH2O/6), na.rm=TRUE, color = "#133F66") +
  geom_point(aes(y = na_if(GH2O, 0)/6), na.rm=TRUE, color = "#133F66") +
  geom_line(aes(y = A), na.rm=TRUE, color = "#A73130") +
  geom_point(aes(y = A), na.rm=TRUE, color = "#A73130") +
  scale_y_continuous(name = "<span style='color:#A73130'>**Net carbon uptake (CO<sub>2</sub>)**</span><br>µM m<sup>-2</sup> s<sup>-1</sup>", limits = c(-1,10), breaks = c(0,2,4,6,8,10), sec.axis = sec_axis(trans=~.*6, name="<span style='color:#133F66'>**Leaf conductance (H<sub>2</sub>O)**</span><br>mM m<sup>-2</sup> s<sup>-1</sup>")) +
  scale_x_datetime(name = "Time", limits = c(start,end + hms("00:10:00")), date_labels="%H:%M") +
  theme_classic() +
  theme(axis.title.y = element_markdown(), legend.position="right")

dev.off()

# PAR
phenotyping.gasexchange %>%
  ggplot(aes(x = DateTime)) +
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_line(aes(y = PARamb), na.rm=TRUE) +
  scale_x_datetime(limits = c(start,end), date_breaks = "3 hours", labels=scales::date_format("%H:%M")) +
  theme_classic()



# *** OpenGreenhouse experiment ***

# Environmental parameters of OpenGreenhouse experiment (used for physiological phenotyping)
read_csv(paste(DATA_ROOT, "Experiments", "OpenGreenhouse.arduino.env.csv", sep = '/')) %>%
  mutate(datetime = ymd_hms(paste(date, timeofday)), date = NULL, timeofday = NULL, condition = as.factor(condition), .before = condition) %>%
  mutate(cycle = if_else(datetime %within% interval("2021-07-14 21:00:00", "2021-07-15 05:59:59"), "Night", if_else(datetime %within% interval("2021-07-15 06:00:00", "2021-07-15 20:59:59"), "Day", NA))) %>%
  group_by(condition, cycle) %>%
  summarise(across(where(is.numeric), mean)) %>%
  drop_na(cycle)

# Water availability (Figure 5a)
opengreenhouse.arduino.soil <- read_csv(paste(DATA_ROOT, "Experiments", "OpenGreenhouse.arduino.soil.csv", sep = '/')) %>%
  mutate(datetime = ymd(date) + hms(timeofday), .before = date) %>%
  mutate(species = case_when(
    str_starts(rep, "F") ~ "Cmajor",
    #str_starts(rep, "M") ~ "Cminor",
    str_starts(rep, "R") ~ "Crosea",
    .default = NA)) %>%
  drop_na(species) %>%
  group_by(datetime, condition, species) %>%
  summarise(SoilMoisture = mean(SoilMoisture, na.rm = TRUE))

pdf(paste(RESULTS_DIR, "figs", "opengreenhouse.water.pdf", sep = '/'), width = 7, height = 4)
opengreenhouse.arduino.soil %>%
  ggplot(aes(x = datetime, y = SoilMoisture, color = species, linetype = condition)) +
  geom_line() +
  scale_color_manual(breaks=c("Cmajor", "Crosea"), values = c("#A73130", "#133D66")) +
  scale_x_datetime(name = "Time") +
  theme_classic()

dev.off()


# Read PhotosynQ data from Excel file and filter for cuttings only
opengreenhouse.photosynq <- read_excel(paste(DATA_ROOT, "Experiments", "OpenGreenhouse.photosynQ.xlsx", sep = '/'), sheet = "All_curated") %>%
  filter(Kind == "cutting") %>%
  mutate(Species = fct_relevel(Species, c("Cmajor", "Cminor", "Crosea"))) %>%
  mutate(Condition = fct_relevel(Condition, c("shaded", "exposed")))

# Light intensity - PAR (Figure 5bc)
opengreenhouse.photosynq.par <- opengreenhouse.photosynq %>%
  group_by(Species, Condition, Timepoint) %>%
  summarise(
    mean_PAR = mean(Light_Intensity_PAR, na.rm = TRUE),   # Mean PAR
    sd_PAR = sd(Light_Intensity_PAR, na.rm = TRUE),       # Standard deviation
    n = sum(!is.na(Light_Intensity_PAR)),                 # Number of observations
    se_PAR = sd_PAR / sqrt(n)                             # Standard error
  )

pdf(paste(RESULTS_DIR, "figs", "opengreenhouse.light.shaded.pdf", sep = '/'), width = 7, height = 4)
opengreenhouse.photosynq.par %>%
  filter(Species != "Cminor") %>%
  filter(Condition == "shaded") %>%
  ggplot(aes(x = Timepoint, y = mean_PAR, color = Species, group = Species)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_PAR - se_PAR, ymax = mean_PAR + se_PAR), width = 0.2, alpha = 0.6) +  # SE error bars
  scale_color_manual(breaks=c("Cmajor", "Crosea"), values = c("#A73130", "#133D66")) +
  theme_classic() +
  labs(title = "PAR shaded", x = "Time of Day", y = "PAR", color = "Species")

dev.off()

pdf(paste(RESULTS_DIR, "figs", "opengreenhouse.light.exposed.pdf", sep = '/'), width = 7, height = 4)
opengreenhouse.photosynq.par %>%
  filter(Species != "Cminor") %>%
  filter(Condition == "exposed") %>%
  ggplot(aes(x = Timepoint, y = mean_PAR, color = Species, group = Species)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_PAR - se_PAR, ymax = mean_PAR + se_PAR), width = 0.2, alpha = 0.6) +  # SE error bars
  scale_color_manual(breaks=c("Cmajor", "Crosea"), values = c("#A73130", "#133D66")) +
  theme_classic() +
  labs(title = "PAR exposed", x = "Time of Day", y = "PAR", color = "Species")

dev.off()


# Plot PhiNPQ over time (Figure 5d)
pdf(paste(RESULTS_DIR, "figs", "opengreenhouse.phinpq.pdf", sep = '/'), width = 7, height = 4)
opengreenhouse.photosynq %>%
  filter(Species != "Cminor") %>%
  ggplot(aes(x = Timepoint, y = PhiNPQ, color = Species)) +
  geom_smooth(method = "loess", span = 0.3) +
  facet_wrap(~ Condition) +
  scale_color_manual(breaks=c("Cmajor", "Crosea"), values = c("#A73130", "#133D66")) +
  theme_minimal() +
  labs(title = "PhiNPQ", x = "Time of Day", y = "PhiNPQ")

dev.off()
