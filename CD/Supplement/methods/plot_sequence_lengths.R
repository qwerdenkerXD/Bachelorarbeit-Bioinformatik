args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    message("Usage: Rscript plot_sequence_lengths.R <fasta-file> <out-file.png>")
    q()
}

# Load necessary libraries
library(seqinr)
library(ggplot2)
library(ggdist)
library(tibble)

# Read the FASTA files
fasta_files <- args[-length(args)]

# get array with array of sequence lengths for each fasta file
sequence_lengths <- lapply(fasta_files, function(f) {
  sequences <- read.fasta(file = f, seqtype = "AA", as.string = TRUE, seqonly=TRUE)
  sapply(sequences, nchar)
})

# create dataframe with fasta files as column names and its sequence lengths as values
groups <- lapply(seq_along(sequence_lengths), function(i) {tibble(value = sequence_lengths[[i]], group = fasta_files[i])})

# prepare dataframe for raincloud plot
data <- dplyr::bind_rows(groups)

# create raincloud plot
plt <- ggplot(data, aes(x = group, y = value)) +
  geom_jitter(position = position_jitter(width = 0.1, seed = 1), size = 0.1, colour = "blue") +
  stat_halfeye(fill = "darkgray", width = 2, position = position_nudge(x = 0.11, y=0)) +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5), axis.ticks.y = element_blank()) +
  scale_y_continuous(n.breaks = 10) +
  labs(y = "Sequenzlängen", x = "") +
  ggtitle("Raincloud Plot der Trainings-Protein Sequenzlängen")

ggsave(args[length(args)], plot = plt, width = 8, height = 3, dpi = 300)
