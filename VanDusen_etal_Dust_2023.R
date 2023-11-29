#code associated with analyses in Van Dusen et al 2023 "Identification of SARS-CoV-2 variants in indoor dust"

#load R packages
library(tidyverse)
library(patchwork)
library(nortest)
library(MixviR)


#Download mutations associated with lineages of interest from outbreak.info database. Filtering for mutations that occur with a frequency of 0.75 in the lineage/sublineage.
outbreakinfo::authenticateUser()

# only considering Omicron lineage BA.1 here as the dust sampling only goes through early 2022.
lineages <- c("B.1.1.7", "B.1.617.2", "BA.1")

#obtain mutations associated with lineages of interest
target_lineage_muts <- outbreakinfo::getMutationsByLineage(pangolin_lineage = lineages, logInfo = FALSE)
target_lineage_muts <- target_lineage_muts %>% 
  filter(prevalence > 0.75)

#rename the lineages and reformat mutation identifiers to match those from MixviR
target_lineage_muts <- target_lineage_muts %>%
  rowwise() %>%
  dplyr::mutate(lineage = case_when(lineage == "B.1.1.7" ~ "Alpha",
                                    lineage == "B.1.617.2" ~ "Delta",
                                    lineage == "BA.1" ~ "Omicron")) %>%
  ungroup() %>%
  dplyr::mutate(ALT_ID = stringr::str_replace_all(mutation, pattern = ":", "_"),
                ALT_ID = toupper(ALT_ID),
                ALT_ID = stringr::str_replace_all(ALT_ID, "DEL", "del"),
                ALT_ID = stringr::str_replace_all(ALT_ID, "INS", "ins")) %>% 
  tidyr::separate(col = ALT_ID, into = c("GENE", "MUTATION"), sep = "_", remove = FALSE) %>%
  select(gene, MUTATION, lineage) %>%
  dplyr::rename("Gene" = "gene", 
                "Mutation" = "MUTATION",
                "Lineage" = "lineage") %>%
  dplyr::distinct()


## Obtain mutations associated with sublineages of interest and format for MixviR

sublineages <- c("AY.1",
                 "AY.3",
                 "AY.4",
                 "AY.5",
                 "AY.6",
                 "AY.7",
                 "AY.9",
                 "AY.13",
                 "AY.14",
                 "AY.16",
                 "AY.17",
                 "AY.20",
                 "AY.21",
                 "AY.23",
                 "AY.24",
                 "AY.25",
                 "AY.26",
                 "AY.27",
                 "AY.28",
                 "AY.29",
                 "AY.30",
                 "AY.32",
                 "AY.33",
                 "AY.34",
                 "AY.35",
                 "AY.36",
                 "AY.37",
                 "AY.38",
                 "AY.39",
                 "AY.41",
                 "AY.42",
                 "AY.43",
                 "AY.44",
                 "AY.45",
                 "AY.46",
                 "AY.47",
                 "AY.48",
                 "AY.51",
                 "AY.53",
                 "AY.54",
                 "AY.55",
                 "AY.56",
                 "AY.57",
                 "AY.58",
                 "AY.59",
                 "AY.60",
                 "AY.61",
                 "AY.62",
                 "AY.64",
                 "AY.65",
                 "AY.66",
                 "AY.67",
                 "AY.68",
                 "AY.69",
                 "AY.70",
                 "AY.71",
                 "AY.72",
                 "AY.73",
                 "AY.74",
                 "AY.75",
                 "AY.76",
                 "AY.77",
                 "AY.78",
                 "AY.79",
                 "AY.80",
                 "AY.81",
                 "AY.83",
                 "AY.84",
                 "AY.85",
                 "AY.86",
                 "AY.87",
                 "AY.88",
                 "AY.90",
                 "AY.91",
                 "AY.92",
                 "AY.93",
                 "AY.94",
                 "AY.95",
                 "AY.98",
                 "AY.99",
                 "AY.100",
                 "AY.101",
                 "AY.102",
                 "AY.103",
                 "AY.104",
                 "AY.105",
                 "AY.106",
                 "AY.107",
                 "AY.108",
                 "AY.109",
                 "AY.110",
                 "AY.111",
                 "AY.112",
                 "AY.113",
                 "AY.114",
                 "AY.116",
                 "AY.117",
                 "AY.118",
                 "AY.119",
                 "AY.120",
                 "BA.1.1.1",
                 "BA.1.1.2",
                 "BA.1.1.3",
                 "BA.1.1.4",
                 "BA.1.1.5",
                 "BA.1.1.6",
                 "BA.1.1.7",
                 "BA.1.1.8",
                 "BA.1.1.9",
                 "BA.1.1.10",
                 "BA.1.1.11",
                 "BA.1.1.12",
                 "BA.1.1.13",
                 "BA.1.1.14",
                 "BA.1.1.15",
                 "BA.1.1.16",
                 "BA.1.1.17",
                 "BA.1.1.18",
                 "BA.1.2",
                 "BA.1.3",
                 "BA.1.4",
                 "BA.1.5",
                 "BA.1.6",
                 "BA.1.7",
                 "BA.1.8",
                 "BA.1.9",
                 "BA.1.10",
                 "BA.1.12",
                 "BA.1.13",
                 "BA.1.14",
                 "BA.1.15",
                 "BA.1.16",
                 "BA.1.17",
                 "BD.1",
                 "BA.1.18",
                 "BA.1.19",
                 "BA.1.20",
                 "BA.1.21",
                 "BA.1.22",
                 "BA.1.23",
                 "BA.1.24")


#obtain mutations associated with sublineages of interest
target_sublineage_muts <- outbreakinfo::getMutationsByLineage(pangolin_lineage = sublineages, logInfo = FALSE)
target_sublineage_muts <- target_sublineage_muts %>% 
  filter(prevalence > 0.75)

#rename the lineages and reformat mutation identifiers to match those from MixviR
target_sublineage_muts_clean <- target_sublineage_muts %>%
  dplyr::mutate(ALT_ID = stringr::str_replace_all(mutation, pattern = ":", "_"),
                ALT_ID = toupper(ALT_ID),
                ALT_ID = stringr::str_replace_all(ALT_ID, "DEL", "del"),
                ALT_ID = stringr::str_replace_all(ALT_ID, "INS", "ins")) %>% 
  tidyr::separate(col = ALT_ID, into = c("GENE", "MUTATION"), sep = "_", remove = FALSE) %>%
  select(gene, MUTATION, query_key) %>%
  dplyr::rename("Gene" = "gene", 
                "Mutation" = "MUTATION",
                "Lineage" = "query_key") %>%
  dplyr::mutate(Sublineage = Lineage) %>%
  rowwise() %>%
  dplyr::mutate(Lineage = case_when(str_detect(Sublineage, "^AY") ~ "Delta",
                                    str_detect(Sublineage, "^(BA\\.1|BD)") ~ "Omicron"))


#combine main lineage mutations with sublineage mutations
target_lineage_muts <- target_lineage_muts %>%
  mutate("Sublineage" = NA)

lineage_muts <- bind_rows(target_lineage_muts, target_sublineage_muts_clean)

#write file with mutations for lineages of interest
write.table(lineage_muts, "lineage_muts_dust.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)



#Analyze with MixviR v3.5.2.1.

#vcf files used as input were generated with Illumina Basespace DRAGEN Covid Lineage App v3.5.1 with the Min Aligner Score set to 50.
dust_muts_clean <- call_mutations(sample.dir = "vcfs/",
                                  reference = "Wuhan", 
                                  name.sep = ".vcf.gz")

#Call lineages with MixviR based on mutation file generated from outbreakinfo above
dust_lineages <- estimate_lineages(muts.df = dust_muts_clean,
                                   lineage.muts = "lineage_muts_dust.csv",
                                   outfile.name = "dust_lineages_out.csv")



dir.create("manuscript_figures")


# Figure 1: Lineage Estimation From A Single Building
dust_lineages_all <- dust_lineages %>%
  mutate(Date = str_replace(Sample, "(.+_)(.+)(_.+)", "\\2")) %>%
  mutate(Location = str_replace(Sample, "(.+)(_.+)(_.+)", "\\1"))

fig1_data <- dust_lineages_all %>%
  filter(Location == "Building01") %>%
  select(Lineage, Estimated_Freq, Date)

fig1_data_sums <- fig1_data %>%
  group_by(Date) %>%
  summarise(pct_sum = sum(Estimated_Freq)) %>%
  select(Date, pct_sum)

fig1_data_proportions <- left_join(x = fig1_data, y = fig1_data_sums, by= "Date")

fig1_data_proportions <- fig1_data_proportions %>%
  mutate(Estimated_Freq = Estimated_Freq/pct_sum) %>%
  select(-pct_sum) %>%
  mutate("Date" = as.Date(Date, "%Y%m%d"))

ggplot() +
  geom_col(data = fig1_data_proportions,
           mapping = aes(x = Date, 
                         y = Estimated_Freq*100,
                         fill = Lineage),
           position = "stack",
           width = 4) +
  labs(x = element_blank(),
       y = "Estimated Percentage") +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, 
                                    linewidth=1)) +
  scale_fill_manual(values = c("Blue", "Orange"))

ggsave(filename = "manuscript_figures/Figure1.pdf", device = "pdf")


# Figure 2: Lineage Estimation From Multiple Floors Of A Single Building
fig2_data <- dust_lineages_all %>%
  #mutate(Location = str_replace(Location, "Library", "Building35-02")) %>%
  filter(str_detect(Location, "Building35")) %>%
  select(Lineage, Estimated_Freq, Date, Location) %>%
  mutate(Floor = str_replace(Location, "(Building35-)(..)", "\\2")) %>%
  mutate(Date = as.Date(Date, "%Y%m%d"))

min_date <- min(fig2_data$Date)
max_date <- max(fig2_data$Date)


floor1 <- fig2_data %>% filter(Floor == "01") %>%
  ggplot() +
  geom_col(aes(x = Date,
               y = Estimated_Freq*100,
               fill = Lineage),
           width = 4,
           show.legend = FALSE) +
  coord_cartesian(xlim = c(min_date, max_date)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, 
                                    linewidth=1)) +
  labs(y = NULL) +
  scale_fill_manual(values = "orange") +
  geom_label(data = tibble("Date" = as.Date("20210615", "%Y%m%d"),
                           "Floor" = "Floor 1"),
             aes(x = Date,
                 y = 85,
                 label = Floor),
             fill = "lightgray")

floor2 <- fig2_data %>% filter(Floor == "02") %>%
  ggplot() +
  geom_col(aes(x = Date,
               y = Estimated_Freq*100,
               fill = Lineage),
           width = 4,
           show.legend = FALSE) +
  coord_cartesian(xlim = c(min_date, max_date)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, 
                                    linewidth=1)) +
  labs(y = NULL) +
  #scale_fill_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("blue")) +
  geom_label(data = tibble("Date" = as.Date("20210615", "%Y%m%d"),
                           "Floor" = "Floor 2"),
             aes(x = Date,
                 y = 85,
                 label = Floor),
             fill = "lightgray")

floor3 <- fig2_data %>% filter(Floor == "03") %>%
  ggplot() +
  geom_col(aes(x = Date,
               y = Estimated_Freq*100,
               fill = Lineage),
           width = 4) +
  coord_cartesian(xlim = c(min_date, max_date)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, 
                                    linewidth=1)) +
  labs(y = "Estimated Percent") +
  scale_fill_manual(values = c("red", "blue", "orange")) +
  geom_col(data = fig2_data,
           aes(x = Date,
               y = 0,
               fill = Lineage),
           width = 4) +
  geom_label(data = tibble("Date" = as.Date("20210615", "%Y%m%d"),
                           "Floor" = "Floor 3"),
             aes(x = Date,
                 y = 85,
                 label = Floor),
             fill = "lightgray")


floor4 <- fig2_data %>% filter(Floor == "04") %>%
  ggplot() +
  geom_col(aes(x = Date,
               y = Estimated_Freq*100,
               fill = Lineage),
           width = 4,
           show.legend = FALSE) +
  coord_cartesian(xlim = c(min_date, max_date)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, 
                                    linewidth=1)) +
  labs(y = NULL) +
  scale_fill_manual(values = c("blue", "orange")) +
  geom_label(data = tibble("Date" = as.Date("20210615", "%Y%m%d"),
                           "Floor" = "Floor 4"),
             aes(x = Date,
                 y = 85,
                 label = Floor),
             fill = "lightgray")


floor11 <- fig2_data %>% filter(Floor == "11") %>%
  ggplot() +
  geom_col(aes(x = Date,
               y = Estimated_Freq*100,
               fill = Lineage),
           width = 4,
           show.legend = FALSE) +
  coord_cartesian(xlim = c(min_date, max_date)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, 
                                    linewidth=1)) +
  labs(y = NULL) +
  scale_fill_manual(values = c("blue", "orange")) +
  geom_label(data = tibble("Date" = as.Date("20210615", "%Y%m%d"),
                           "Floor" = "Floor 11"),
             aes(x = Date,
                 y = 85,
                 label = Floor),
             fill = "lightgray")


floor11 / floor4 / floor3 / floor2 / floor1
ggsave(filename = "manuscript_figures/Figure2.pdf", device = "pdf")




# Clinical (Saliva) vs Dust Data Comparison

#Get estimated weekly lineage frequencies from saliva data (obtained from GISAID)

saliva_data <- read_csv("OSU_saliva.csv") %>%
  select(-1) %>%
  mutate("Lineage" = pangolin_lineage) %>%
  mutate("Lineage" = str_replace(Lineage, "B.1.1.7", "Alpha")) %>%
  mutate("Lineage" = str_replace(Lineage, "AY.+", "Delta")) %>%
  mutate("Lineage" = str_replace(Lineage, "BA\\.1.*", "Omicron")) %>%
  mutate("Lineage" = str_replace(Lineage, "BA\\.2.*", "Omicron")) %>%
  mutate("Lineage" = str_replace(Lineage, "BA\\.3.*", "Omicron")) %>%
  mutate("Lineage" = str_replace(Lineage, "BA\\.4.*", "Omicron")) %>%
  mutate("Lineage" = str_replace(Lineage, "BA\\.5.*", "Omicron"))

target_lineages <- c("Alpha", "Delta", "Omicron")

saliva_data <- saliva_data %>%
  rowwise() %>%
  mutate(Lineage = ifelse(Lineage %in% target_lineages, Lineage, "Other"))

#group by week
saliva_data <- saliva_data %>%
  mutate("Week" = lubridate::ceiling_date(date, "week"))

#Get weekly sample sizes and add in missing lineage/week combinations for plots.
week_cts <- saliva_data %>% 
  group_by(Week) %>%
  summarise("week_n" = n())

#add in missing lineage/week combinations
complete_df <- complete(saliva_data, Lineage, Week)

#merge in the weekly sample sizes - these are for any lineage sampled during that week.
complete_df <- left_join(x = complete_df, y = week_cts, by = "Week")

#calculate frequencies of each lineage each week (only calculate when the week had at least 5 lineages sampled; N >= 5)
saliva_freqs <- complete_df %>%
  group_by(Lineage, Week) %>%
  summarise("pct" = ifelse(unique(week_n) > 5, sum(!is.na(gisaid_epi_isl))/unique(week_n), NA),
            "n" = unique(week_n)) %>%
  filter(!is.na(pct)) %>%
  mutate("Dataset" = "OSU_saliva")




#Get estimated weekly lineage frequencies from dust data for the OSU campus community

#Bin sample dates into weeks (week ending), average over lineage frequencies for each week (to get estimate of the OSU community freq of the lineage for the week), and filter for Alpha, Delta, and Omicron.

dust_lineages_weekly <- dust_lineages_all %>%
  mutate(Date = as.Date(Date,  "%Y%m%d")) %>%
  mutate(Week = lubridate::ceiling_date(Date, "week"))

dust_lineages_weekly <- dust_lineages_weekly %>%
  group_by(Week, Lineage) %>%
  summarise(pct = mean(Estimated_Freq),
            n = n()) %>%
  ungroup() %>%
  select(Week, Lineage, pct, n) %>%
  mutate("Dataset" = "Dust")

#convert dust values to proportions within each week. For example, if we had a 
#building that had 100% Delta in a given week, and a separate building that had 100% omicron
#for that same week, we'd call this 50% Delta and 50% Omicron for the population.

dust_lineages_proportions_sum <- dust_lineages_weekly %>%
  group_by(Week) %>%
  summarise(pct_sum = sum(pct)) %>%
  select(Week, pct_sum)

dust_lineages_proportions <- left_join(x = dust_lineages_weekly, y = dust_lineages_proportions_sum, by= "Week")

dust_lineages_proportions <- dust_lineages_proportions %>%
  mutate(pct = pct/pct_sum) %>%
  select(-pct_sum)

#combine the dust data with the OSU saliva data
saliva_dust_data <- bind_rows(saliva_freqs, dust_lineages_proportions)

#filter for dates covered by both datasets and target lineages
saliva_dust_data <- saliva_dust_data %>% 
  filter(Week < as.Date("20220308", "%Y%m%d")) %>%
  filter(Week > as.Date("20210301", "%Y%m%d")) %>%
  filter(Lineage != "Other")



# Figure 3: Estimates Of Community-Level SARS-CoV-2 Lineage Frequencies From Saliva And Dust

ggplot() +
  geom_col(data = saliva_dust_data %>% filter(Dataset == "Dust"), 
           mapping = aes(x =  Week, 
                         y = pct*100, 
                         fill = Lineage), 
           alpha = 0.8, 
           color = NA, 
           width = 4,
           show.legend = TRUE) +
  theme_classic() +
  scale_fill_manual(values = c("Red", "Blue", "Orange")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  geom_line(data = saliva_dust_data %>% filter(Dataset == "OSU_saliva"), 
            mapping = aes(x = Week, 
                          y = pct*100, 
                          color = Lineage), 
            linetype = "21", 
            size = 1,
            show.legend = FALSE) +
  scale_color_manual(values = c("Red", "Blue", "Orange")) +
  labs(x = element_blank(),
       y = "Estimated Percentage") +
  theme(axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 12))


#Statistical Analyses
#Evaluate whether a non-random relationship exists between the estimates from the saliva data and the estimates from the dust data. 

# Test For Normality

#apply Kolmogorov-Smirnov test with Lilliefors correction to test for normality

#all Delta and Omicron
all_for_ks <- saliva_dust_data %>%
  filter(Lineage == "Delta" | Lineage == "Omicron") %>%
  pull(pct)

lillie.test(all_for_ks)


#plot distribution of all Delta and Omicron data

saliva_dust_data %>%
  filter(Lineage == "Omicron" | Lineage == "Delta") %>%
  ggplot(aes(x = pct)) +
  geom_histogram()


#The Lilliefor's test suggests the data are not normally distributed, which is consistent with the observation that most observations cluster around the extremes of 0 or 1.

# Spearman Analysis: Delta Lineage
#Filter for Delta observations 
delta_spearman <- saliva_dust_data %>%
  filter(Lineage == "Delta")

#select only weeks that have Delta frequency estimates from both the saliva and dust datasets
delta_duplicated <- delta_spearman$Week[duplicated(delta_spearman$Week)]

delta_spearman <- delta_spearman %>%
  filter(Week %in% delta_duplicated)

delta_spearman_saliva <- delta_spearman %>%
  ungroup() %>%
  filter(Dataset == "OSU_saliva") %>%
  select(Week, pct) %>%
  rename("pct_saliva" = "pct")

delta_spearman_dust <- delta_spearman %>%
  ungroup() %>%
  filter(Dataset == "Dust") %>%
  select(Week, pct) %>%
  rename("pct_dust" = "pct")

delta_spearman <- inner_join(x = delta_spearman_saliva, y = delta_spearman_dust, by = "Week")

#add in nominal small value to all observations to break ties

vals_to_add_saliva <- rnorm(20, mean = 0, sd = .00001)
vals_to_add_dust <- rnorm(20, mean = 0, sd = .00001)

delta_spearman <- delta_spearman %>%
  mutate(pct_saliva = pct_saliva+vals_to_add_saliva) %>%
  mutate(pct_dust = pct_dust+vals_to_add_dust)

cor(x = delta_spearman$pct_saliva, 
    y = delta_spearman$pct_dust,
    method = "spearman")

cor.test(x = delta_spearman$pct_saliva, 
         y = delta_spearman$pct_dust,
         method = "spearman")



# Spearman Analysis: Omicron Lineage
#Filter for Omicron observations 
omicron_spearman <- saliva_dust_data %>%
  filter(Lineage == "Omicron")

#select only weeks that have Omicron frequency estimates from both the saliva and dust datasets
omicron_duplicated <- omicron_spearman$Week[duplicated(omicron_spearman$Week)]

omicron_spearman <- omicron_spearman %>%
  filter(Week %in% omicron_duplicated)

omicron_spearman_saliva <- omicron_spearman %>%
  ungroup() %>%
  filter(Dataset == "OSU_saliva") %>%
  select(Week, pct) %>%
  rename("pct_saliva" = "pct")

omicron_spearman_dust <- omicron_spearman %>%
  ungroup() %>%
  filter(Dataset == "Dust") %>%
  select(Week, pct) %>%
  rename("pct_dust" = "pct")

omicron_spearman <- inner_join(x = omicron_spearman_saliva, y = omicron_spearman_dust, by = "Week")

#add in nominal small value to all observations to break ties

vals_to_add_saliva <- rnorm(8, mean = 0, sd = .00001)
vals_to_add_dust <- rnorm(8, mean = 0, sd = .00001)

omicron_spearman <- omicron_spearman %>%
  mutate(pct_saliva = pct_saliva+vals_to_add_saliva) %>%
  mutate(pct_dust = pct_dust+vals_to_add_dust)


cor(x = omicron_spearman$pct_saliva, 
    y = omicron_spearman$pct_dust,
    method = "spearman")

cor.test(x = omicron_spearman$pct_saliva, 
         y = omicron_spearman$pct_dust,
         method = "spearman")














