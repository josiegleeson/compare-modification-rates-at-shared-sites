# Compare m6A modification rates between shared sites in different transcript isoforms encoded by the same gene

Applies a two-proportions z-test on the annotated output from m6anet to test for differences in the modification rates within multiple transcript isoforms that encode the same genomic position.

Requires the output file "annotated_modification_sites.csv" from running: https://github.com/josiegleeson/annotate-m6anet-output

Download the Rscript and run as follows:
```
Rscript compare_modification_rates_at_shared_genomic_positions.R annotated_modification_sites.csv group_id
Rscript compare_modification_rates_at_shared_genomic_positions.R test_data_annotated_modification_sites.csv caud
```
