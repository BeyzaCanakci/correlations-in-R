# required libraries
library(reshape2)   # For reshaping data
library(dplyr)      # For data manipulation
library(ggplot2)    # For plotting
library(stringr)    # For string operations
library(tidyr)      # For data tidying

# Load theoretical FoldX data (already log-transformed and -1 subtracted)
main_dir2 <- "/Users/b/Desktop/wide_foldx.csv"
df <- read.csv(main_dir2)

# Load experimental expression data
main_path_ex <- "/Users/be/Desktop/results/expression.csv"
df_ex <- read.csv(main_path_ex)

# Initialize empty data.frames to store combined data
combined <- data.frame()
combined_ex <- data.frame()
combined <- rbind(combined, df)        # Append FoldX data
combined_ex <- rbind(combined_ex, df_ex)  # Append experimental data

# Rename first column to 'Position' for consistency
names(combined)[1] <- "Position"
names(combined_ex)[1] <- "Position"

# Convert FoldX data from wide to long format
long_df_1 <- pivot_longer(
  data = combined,
  cols = -Position,
  names_to = "Mutation",
  values_to = "Bind_Avg_Energy"
)

# Merge 'Position' and 'Mutation' into a new Position identifier
long_df_1 <- long_df_1 %>%
  mutate(Position = paste0(Position, Mutation)) %>% 
  select(-Mutation)

# Convert experimental data from wide to long format
long_df_2 <- pivot_longer(
  data = combined_ex,
  cols = -Position,
  names_to = "Mutation",
  values_to = "Bind_Avg_Energy"
)

# Merge 'Position' and 'Mutation' into a new Position identifier
long_df_2 <- long_df_2 %>%
  mutate(Position = paste0(Position, Mutation)) %>% 
  select(-Mutation)

# Ensure column names are aligned
names(combined_ex)[2] <- "Bind_Avg_Energy"
names(long_df_1)[2] <- "Bind_Avg_Energy"

# Round energy values to 2 decimal places
combined_ex$Bind_Avg_Energy <- round(combined_ex$Bind_Avg_Energy, 2)
long_df_1$Bind_Avg_Energy <- round(long_df_1$Bind_Avg_Energy, 2)


# Merge theoretical and experimental data on Position
combined_data_log <- merge(long_df_1, combined_ex, by = "Position", suffixes = c("_konum1", "_konum2"))
combined_data_log <- na.omit(combined_data_log)  # Remove NA rows

# Keep only the first entry per position (handle duplicates if any)
df_unique_log <- combined_data_log %>%
  group_by(Position) %>%
  slice(1) %>%
  ungroup()

# Define interface residues to filter
interface_residue <- c("G446", "Y449", "Y453", "L455", "F456", "Y473", "A475", 
                       "G476", "F486", "S477", "N487", "Y489", "Q493", "G496", 
                       "Q498", "T500", "N501", "G502", "Y505")

# Filter only interface residues from merged data
filtered_df_log <- df_unique_log[grepl(paste0("^(", paste(interface_residue, collapse = "|"), ")"), df_unique_log$Position), ]

# Calculate Pearson correlation
correlation <- cor(filtered_df_log$Bind_Avg_Energy_konum1, filtered_df_log$Bind_Avg_Energy_konum2)
correlation_text <- paste("Correlation: ", round(correlation, 3))

# Create scatter plot with diagonal line
ggplot(filtered_df_log, aes(x = Bind_Avg_Energy_konum1, y = Bind_Avg_Energy_konum2, label = Position)) +
  geom_point() +  # Scatter plot
  geom_abline(intercept = 0, slope = ifelse(correlation < 0, -1, 1), color = "red", linetype = "dashed") +  # Diagonal line showing trend
  labs(title = "Correlation of energy of interface residue between foldx vs LLR",
       x = "dG of interface residue of complex from foldx",
       y = "dG of interface residue") +
  theme_minimal() +
  annotate("text", x = -1.9, y = 1.2, label = correlation_text, size = 4, color = "purple")  # Show correlation value

# Print correlation coefficient
print(paste("Correlation coefficient: ", correlation_text))
