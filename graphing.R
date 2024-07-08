ase <- data.table::fread("temp/debug/7_all_splicing_events.tsv")
samplefile <- data.table::fread("temp/debug/0_samplefile.tsv")

canonical <- ase[gene == "CAPN3" & annotated == "canonical" & SJ_IR == "SJ"]

tests <- samplefile[sampletype == "test", sampleID]
exclusion_criteria <- samplefile[sampleID %in% tests, c('familyID', 'genes', 'transcript')]
controls <- samplefile[!(sampleID %in% tests) & !(familyID %in% exclusion_criteria$familyID) & !(genes %in% exclusion_criteria$genes) & !(transcript %in% exclusion_criteria$transcript), sampleID]

grep("count_", test_count_cols)

test_count_cols <- names(canonical)[names(canonical) %in% sapply(tests, function(x) paste0("count_", x))]
control_count_cols <- names(canonical)[names(canonical) %in% sapply(controls, function(x) paste0("count_", x))]

test_percentage_cols <- names(canonical)[names(canonical) %in% sapply(tests, function(x) paste0("pct_", x))]
control_percentage_cols <- names(canonical)[names(canonical) %in% sapply(controls, function(x) paste0("pct_", x))]

test_data <- canonical[, c(..test_percentage_cols)]
control_data <- canonical[, c(..control_percentage_cols)]


library(ggplot2)

long_data <- data.table::as.data.table(reshape2::melt(t(control_data)))
mean_data <- long_data[, .(mean_value = mean(value)/2), by = Var2]
line_data <- data.table::as.data.table(reshape2::melt(t(test_data)))
long_data$Var2 <- as.factor(long_data$Var2)
line_data$Var2 <- as.factor(line_data$Var2)
mean_data$Var2 <- as.factor(mean_data$Var2)
line_data$Var1 <- "proband"
mean_data$Var1 <- "estimated single control allele"
point_data <- line_data
line_data <- rbind(mean_data[,c(3,1,2)], line_data, use.names = F)

names(line_data)

ggplot(long_data, aes(x = Var2, y = value, fill = Var2)) +
  geom_violin(trim = FALSE, fill = "grey", color = "grey") +
  #geom_line(data = mean_data, aes(x = Var2, y = mean_value, group = 1), color = "red") +
  #geom_line(data = line_data, aes(x = Var2, y = value, group = 1), color = "black") +
  geom_point(data = point_data, aes(x = Var2, y = value, group = 1), color = "black") +
  geom_line(data = line_data, aes(x = Var2, y = mean_value, group = Var1, color = Var1)) +
  scale_color_manual(values = c("proband" = "black", "estimated single control allele" = "red")) +
  guides(fill = FALSE) +
  labs(x = "Intron Number",
       y = "Canonical Splicing / All Splicing at Intron") +    
  theme_minimal()
  

