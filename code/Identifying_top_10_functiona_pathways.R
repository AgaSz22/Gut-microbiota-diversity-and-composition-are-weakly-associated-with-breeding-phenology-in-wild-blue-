
library(tidyverse)

# Read the TSV file
df <- read.delim("diet_unstratified.relab.tsv", header = TRUE, sep = "\t", check.names = FALSE)

# Rename columns for clarity
colnames(df) <- c("Pathway", "EC", "LC", "EF", "LF")

# Convert abundance columns to numeric
df <- df %>%
  mutate(across(c(EC, LC, EF, LF), as.numeric))

### since these lines calculate the absolute difference, you need to check out the idrection in the original file

df_diff <- df %>%
  mutate(
    diff_chicks  = abs(EC - LC),   # Early vs Late Chicks
    diff_females = abs(EF - LF)    # Early vs Late Females
  )
