# jhu-covid-pipeline

Analysis code and viral genomics pipeline associated with *Genomic Diversity of SARS-CoV-2 During Early Introduction into the United States National Capital Region* by Thielen et al. 2020.


Please note that line list metadata needed for Figures 1-3 (`diagnostic_metadata.csv` and `sample_metadata.csv`) are not provided due to patient privacy restrictions, and sequencing data is not included per GISAID agreements (`all_submitted_jhu_sequences.fasta` and `nextstrain_ncov_local_metadata.tsv`). Associated sample metadata is available on NCBI GenBank under accessions MT509452-MT509493 and MT646049-MT646120, and all sequence data (accessions listed in `tree_ids_authors.csv` and `global_geo_accessions.csv`) can be downloaded from [GISAID](gisaid.org).

`jhhs_aggregate_march.csv` is also provided as Table S1 in the associated manuscript; `variants_jhhs_samples.csv` is provided as Table S2; `samples_zip3_clade.csv` is provided as Table S3; `tree_ids_authors.csv` is provided as Table S4; and `global_geo_accessions.csv` is provided as Table S5.

All analyses use [Wuhan-Hu-1/2019](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) as the reference genome (`reference.fasta`).

Data in `time_series_covid19_confirmed_US.csv` provided by: Dong E, Du H, Gardner L. An interactive web-based dashboard to track COVID-19 in real time. Lancet Infect Dis 2020; 20: 533–4.
