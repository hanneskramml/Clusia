library(tidyverse)
library(readxl)
library(ggtext)


# *** Climate chamber experiment ***

# Gas exchange measurement
tmp.file <- paste(DATA_ROOT, "Physiology", "ClimateChamber.gasexchange.csv", sep = '/')
src.physiology <- read_delim(tmp.file, delim = ';', col_names = read_delim(tmp.file, delim = ';') %>% names(), col_select = Date:wa, skip = 2) %>%
  mutate(DateTime = ymd(Date) + hms(Time), .before = Date) %>%
  mutate_at(c("Code", "Object", "Status"), as.factor) %>%
  filter(Code == "MP_001") %>%
  slice(which(row_number() %% 2 == 1)) %>%
  mutate(A = A-0.7, GH2O = GH2O-2) %>%  # Zero point calibration
  mutate(GH2O = if_else(GH2O < 0, 0, GH2O)) # Issue of negative values, see Methods


# Clusia major (Figure 1b)
start <- ymd_hms("2025-05-15 21:00:00")
end <- ymd_hms("2025-05-16 21:00:00")
pdf(paste(RESULTS_DIR, "figs", "gasexchange.major.pdf", sep = '/'), width = 7, height = 3)

src.physiology %>%
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

# Clusia rosea (Figure 1c)
start <- ymd_hms("2025-05-17 21:00:00")
end <- ymd_hms("2025-05-18 21:00:00")
pdf(paste(RESULTS_DIR, "figs", "gasexchange.rosea.pdf", sep = '/'), width = 7, height = 3)

src.physiology %>%
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
src.physiology %>%
  ggplot(aes(x = DateTime)) +
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_line(aes(y = PARamb), na.rm=TRUE) +
  scale_x_datetime(limits = c(start,end), date_breaks = "3 hours", labels=scales::date_format("%H:%M")) +
  theme_classic()



# *** OpenGreenhouse experiment ***

# Environmental parameters of OpenGreenhouse experiment (used for climate chambers)
read_csv(paste(DATA_ROOT, "Physiology", "OpenGreenhouse.arduino.env.csv", sep = '/')) %>%
  mutate(datetime = ymd_hms(paste(date, timeofday)), date = NULL, timeofday = NULL, condition = as.factor(condition), .before = condition) %>%
  mutate(cycle = if_else(datetime %within% interval("2021-07-14 21:00:00", "2021-07-15 05:59:59"), "Night", if_else(datetime %within% interval("2021-07-15 06:00:00", "2021-07-15 20:59:59"), "Day", NA))) %>%
  group_by(condition, cycle) %>%
  summarise(across(where(is.numeric), mean)) %>%
  drop_na(cycle)

# Water availability (Figure 5a)
pdf(paste(RESULTS_DIR, "figs", "opengreenhouse.soilwater.pdf", sep = '/'), width = 7, height = 4)
read_csv(paste(DATA_ROOT, "Physiology", "OpenGreenhouse.arduino.soil.csv", sep = '/')) %>%
  mutate(datetime = ymd(date) + hms(timeofday), .before = date) %>%
  mutate(species = case_when(
    str_starts(rep, "F") ~ "Cmajor",
    #str_starts(rep, "M") ~ "Cminor",
    str_starts(rep, "R") ~ "Crosea",
    .default = NA)) %>%
  drop_na(species) %>%
  group_by(datetime, condition, species) %>%
  summarise(SoilMoisture = mean(SoilMoisture, na.rm = TRUE)) %>%
  ggplot(aes(x = datetime, y = SoilMoisture, color = species, linetype = condition)) +
  geom_line() +
  scale_color_manual(breaks=c("Cmajor", "Crosea"), values = c("#A73130", "#133D66")) +
  scale_x_datetime(name = "Time") +
  theme_classic()

dev.off()


# Read PhotosynQ data from Excel file and filter for cuttings only
src.photosynq <- read_excel(paste(DATA_ROOT, "Physiology", "OpenGreenhouse.photosynQ.xlsx", sep = '/'), sheet = "All_curated") %>%
  filter(Kind == "cutting") %>%
  mutate(Species = fct_relevel(Species, c("Cmajor", "Cminor", "Crosea"))) %>%
  mutate(Condition = fct_relevel(Condition, c("shaded", "exposed")))

# Light intensity - PAR (Figure 5bc)
feature.par <- src.photosynq %>%
  group_by(Species, Condition, Timepoint) %>%
  summarise(
    mean_PAR = mean(Light_Intensity_PAR, na.rm = TRUE),   # Mean PAR
    sd_PAR = sd(Light_Intensity_PAR, na.rm = TRUE),       # Standard deviation
    n = sum(!is.na(Light_Intensity_PAR)),                 # Number of observations
    se_PAR = sd_PAR / sqrt(n)                             # Standard error
  )

pdf(paste(RESULTS_DIR, "figs", "opengreenhouse.light.shaded.pdf", sep = '/'), width = 7, height = 4)
feature.par %>%
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
feature.par %>%
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
pdf(paste(RESULTS_DIR, "figs", "opengreenhouse.phiNPQ.pdf", sep = '/'), width = 7, height = 4)
src.photosynq %>%
  filter(Species != "Cminor") %>%
  ggplot(aes(x = Timepoint, y = PhiNPQ, color = Species)) +
  geom_smooth(method = "loess", span = 0.3) +
  facet_wrap(~ Condition) +
  scale_color_manual(breaks=c("Cmajor", "Crosea"), values = c("#A73130", "#133D66")) +
  theme_minimal() +
  labs(title = "PhiNPQ", x = "Time of Day", y = "PhiNPQ")

dev.off()
