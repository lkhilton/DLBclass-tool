library(tidyverse)

samps <- read_tsv("DLBCL_GenomeCanada/DLBCL_GenomeCanada_sampleset.tsv", col_names = "sample_id")

svs <- read_tsv(
    "/projects/rmorin/projects/gambl-repos/gambl-lhilton/results/all_the_things/lymphgen-1.0/create_inputs/complete_BCL2_BCL6_fusions_for_lymphgen.tsv",
    col_names = c("sample_id", "SV.BCL6", "SV.BCL2", "SV.MYC"),
    skip = 1
)

svs <- svs %>%
    filter(sample_id %in% samps$sample_id) %>%
    mutate(across(-matches("sample_id"), ~ case_when(
        .x == "POS" ~ 3,
        .x == "NEG" ~ 0
    ))) %>%
    column_to_rownames("sample_id")

svs <- rownames_to_column(data.frame(t(svs[samps$sample_id, ]), check.names = FALSE), "classifier_name")

write_tsv(svs, "DLBCL_GenomeCanada/DLBCL_GenomeCanada.15Jan2025.SV.GSM.tsv")
