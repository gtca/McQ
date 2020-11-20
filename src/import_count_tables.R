# 
# Combine read/UMI counts for all the samples in a batch
# 
# Examples:
#   # combine for all the batches
#   Rscript src/import_count_tables.R     # combine for all the batches
#   # combine for all the batches that start with J2.*
#   Rscript src/import_count_tables.R J2  
#   # combine custom_counts.tsv.gz files for all the batches that start with J2.*
#   Rscript src/import_count_tables.R J2 custom_counts.tsv.gz
# 

library(tidyverse)
library(glue)

input_dir <- "data/interim/seq/"
output_dir <- "data/processed/counts/"

args <- commandArgs(trailingOnly = TRUE)
batch_pattern <- NULL
if (length(args) > 0) {
  batch_pattern <- args[1]
  message(glue("Using only samples starting with {batch_pattern}"))
  if (!endsWith(batch_pattern, ".*")) {
    batch_pattern <- glue("^{batch_pattern}.*")
  }

  if (length(args) > 1) {
    filename <- args[2]
  } else {
    filename <- "read_umi_counts_clean.tsv.gz"
  }
}

# Get all the samples names
# Do not include temporary or legacy folders starting with _
samples_names <- input_dir %>% list.files(pattern = batch_pattern) %>% basename %>% .[!startsWith(., "_")]

# Get all the batches names
get_batch_from_sample <- function(e) {
	# Account for batch names for newer plates with shorter names
	if (length(strsplit(e, "-")[[1]]) < 3) {	
		e %>% str_split("-") %>% map_chr(~ .[1])
	} else {
		e %>% str_split("-", n = 3) %>% map_chr(~ paste(.[1], .[2], sep="-"))
	}
}
batches_names <- samples_names %>% get_batch_from_sample %>% unique


# Load metadata
meta <- read_tsv("data/raw/tables/sample_bcs.txt", col_types = cols())

# Load all the samples
counts <- do.call(rbind, 
		lapply(samples_names, function(sample_name) {
			file.path(input_dir, sample_name, filename) %>%
				read_tsv(col_types = cols()) %>%
				mutate(sample_name = sample_name,
	   				   batch = sample_name %>% get_batch_from_sample) %>%
				rename(barcode = cell) %>%
				left_join(meta %>% select(barcode = `Barcode Sequence`, barcode_id = BarcodeID) %>% unique, by = "barcode") %>%
				mutate(sample_id = paste(sample_name, barcode_id, sep = "-"),
                )
		})
)


# Create an output folder for the batch
for (batch_name in batches_names) {
	dir.create(file.path(output_dir, batch_name), showWarnings = FALSE, recursive = TRUE)
	output_file <- file.path(output_dir, batch_name, filename)

	# Save counts table
	counts %>% 
		dplyr::filter(batch == batch_name) %>%
		write_tsv(output_file)

	message(paste0("Done writing ", output_file))
}
