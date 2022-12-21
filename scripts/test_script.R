install.packages("tidyverse")

library(tidyverse)


# read file containing efficiency
efficiency <- read_csv("data/efficiency.csv")

colnames(efficiency) <- c("Target", "efficiency")

# set params
control_gene_1 <- "RPL7"
control_gene_2 <- "RPL32"

control_sample <- "brain"

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
df <- df %>% filter(Sample %in% c("brain", "liver", "muscle"))

# extract controls
control_1 <- df %>% filter(Target == control_gene_1) %>% distinct() %>%
  select(Sample, Biological.Set.Name, mean_control_1 = Cq.Mean)
control_2 <- df %>% filter(Target == control_gene_2) %>% distinct() %>%
  select(Sample, Biological.Set.Name, mean_control_2 = Cq.Mean)

# extract cq for control sample on control genes
control_1_ctrsample <- control_1 %>% filter(Sample == control_sample) %>%
  select(Biological.Set.Name, mean_control_1_sample = mean_control_1)
control_2_ctrsample <- control_2 %>% filter(Sample == control_sample) %>%
  select(Biological.Set.Name, mean_control_2_sample = mean_control_2)

# extract sample to use as control
df_bio_control <- df %>% filter(Sample == control_sample) %>% distinct() %>%
  select(Target, Biological.Set.Name, cq_control_sample = Cq.Mean)

df_samples <- df %>% filter(!Target %in% c(control_gene_1, control_gene_2)) %>%
  filter(!Sample %in% control_sample) %>%
  distinct()


# show bad replicates with Std Dev > 5%
df <- df %>% mutate(badrep = case_when(Cq.Std..Dev > (Cq.Mean/20) ~ "above_5%",
                                 TRUE ~ "OK"))

# combine dataframes into 1
df_samples <- df_samples %>% left_join(control_1) %>% left_join(control_2) %>% left_join(df_bio_control) %>% left_join(control_1_ctrsample) %>% left_join(control_2_ctrsample) %>% left_join(efficiency)


# start calculation here --------------------------------------------------

calculated_df <- df_samples %>% mutate(RE = (efficiency)^(cq_control - Cq.Mean)/
                        sqrt(efficiency_ctrl_gene_1^(mean_control_1_sample-mean_control_1) * efficiency_ctrl_gene_2^(mean_control_2_sample-mean_control_2)))



# create graph ------------------------------------------------------------
calculated_df %>% ggplot(aes(x = Target, y = log2(RE), fill = Sample)) + geom_bar(stat = "identity", position = "dodge")


# calculate SD




