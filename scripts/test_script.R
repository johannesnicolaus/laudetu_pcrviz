install.packages("tidyverse")
library(tidyverse)

# read files containing relationship between primer code and gene id/ gene name
primercode <- tibble(
  Name = c("aco2 jg4641", "aco2 jg13709",
             "eno2 jg7948", "eno2 jg32633",
             "gapdh jg32651", "gapdh jg7921",
             "gpat4 jg17211", "gpat4 jg37465",
             "idh2 jg1791", "idh2 21942",
             "idh3a jg22789", "idh3a jg1859",
             "idh3g jg32257"),
  Target = c("QC1", "QF1",
           "QE1", "QK1",
           "QD1", "QD2",
           "QG1", "QM1",
           "QA1", "QH1",
           "N3", "QB2",
           "QJ1")
)



# read file containing efficiency
efficiency <- read_csv("data/efficiency.csv")

efficiency[[2]] <- (efficiency[[2]]/100)+1

colnames(efficiency) <- c("Target", "efficiency")

# TOSET: set params
control_gene_1 <- "RPL7"
control_gene_2 <- "RPL32"

# TOSET: samples to analyze
samples <- c("brain", "liver", "muscle")
# samples <- c("dmso", "lxr10-7m", "lxr10-8m", "lxr10-9m")

# TOSET: samples to be used as control
control_sample <- "brain"
control_biologicalsetname <- 2

efficiency_ctrl_gene_1 <- efficiency %>% filter_at(1, all_vars(. == control_gene_1)) %>% .[[2]]
efficiency_ctrl_gene_2 <- efficiency %>% filter_at(1, all_vars(. == control_gene_2)) %>% .[[2]]

# list files
files <- list.files("data/", full.names = T, recursive = T)

files <- files %>% str_subset("Cq")

# read files
df <- files %>% lapply(function(x){read_csv(x, name_repair = make.names)})

df <- df %>% lapply(function(x){
  x %>%
    select(Target, Content, Sample, Biological.Set.Name, Cq.Mean, Cq.Std..Dev)
})


df <- df %>% bind_rows(.id = "number")

# remove positive controls and negative controls
df <- df %>% filter(!Content %in% c("NTC", "Pos Ctrl"))

# divide by experiments (TODO: make automatic)
df <- df %>% filter(Sample %in% samples)

# extract controls
control_1 <- df %>% filter(Target == control_gene_1) %>% distinct() %>%
  select(Sample, Biological.Set.Name, mean_control_1 = Cq.Mean)
control_2 <- df %>% filter(Target == control_gene_2) %>% distinct() %>%
  select(Sample, Biological.Set.Name, mean_control_2 = Cq.Mean)

# extract cq for control sample on control genes
control_1_ctrsample <- control_1 %>% filter(Sample == control_sample, Biological.Set.Name == control_biologicalsetname) %>%
  select(Biological.Set.Name, mean_control_1_sample = mean_control_1, -Biological.Set.Name)
control_2_ctrsample <- control_2 %>% filter(Sample == control_sample, Biological.Set.Name == control_biologicalsetname) %>%
  select(Biological.Set.Name, mean_control_2_sample = mean_control_2, -Biological.Set.Name)

# extract sample to use as control
df_bio_control <- df %>% filter(Sample == control_sample, Biological.Set.Name == control_biologicalsetname) %>%
  distinct() %>%
  select(Target, Biological.Set.Name, cq_control_sample = Cq.Mean, -Biological.Set.Name)

# remove samples that are used for control
df_samples <- df %>% mutate(is.control = case_when(
  (Sample == control_sample & Biological.Set.Name == control_biologicalsetname) ~ "control",
  TRUE ~ "sample"
)) %>% filter(is.control == "sample") %>% select(-is.control) %>% distinct()

# show bad replicates with Std Dev > 5%
df <- df %>% mutate(badrep = case_when(Cq.Std..Dev > (Cq.Mean/20) ~ "above_5%",
                                 TRUE ~ "OK"))

# remove erroneous data
df_samples <- df_samples %>% filter(Cq.Std..Dev != 0)
df_bio_control <- df_bio_control %>% filter(cq_control_sample != 0)

# combine dataframes into 1
df_samples <- df_samples %>% left_join(control_1) %>% left_join(control_2) %>% left_join(df_bio_control) %>% left_join(efficiency) %>%
  mutate(mean_control_1_sample = control_1_ctrsample$mean_control_1_sample, mean_control_2_sample = control_2_ctrsample$mean_control_2_sample)

# start calculation here --------------------------------------------------
calculated_df <- df_samples %>% mutate(RE = ((efficiency)^(cq_control_sample - Cq.Mean))/
                        sqrt((efficiency_ctrl_gene_1)^(mean_control_1_sample-mean_control_1) * (efficiency_ctrl_gene_2)^(mean_control_2_sample-mean_control_2)))

# calculate sd
calculated_df_summary <- calculated_df %>% group_by(Target, Sample) %>% summarise(mean_RE = mean(RE), sd = sd(RE)) %>% ungroup()

# create graph ------------------------------------------------------------

calculated_df_summary %>% filter(!Target %in% c(control_gene_1, control_gene_2)) %>%
  inner_join(primercode) %>%
  ggplot(aes(x = Sample, y = mean_RE, fill = Sample)) +
  facet_wrap(~ Name, ncol = 4) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=mean_RE-sd, ymax=mean_RE+sd), width=.2, position=position_dodge(.9)) +
  cowplot::theme_cowplot(10) + theme(axis.text.x = element_blank()) + ylab("Relative expression (A.U.)") +
  cowplot::panel_border(color = "black")


