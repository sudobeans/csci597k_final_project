data_files="chea3Results data graphs htseqModResults pyFilteredResults unfilteredResults results sorted_bam Homo_sapiens.GRCh38.113.gtf"
code_files="STAR randall_j_eck_r_files_unused transcriptionFactors.R microglia_tfna.ipynb installPkgs.R deseq2.R README.md"

zip -1 -r data_files.zip $data_files
zip -1 -r code_files.zip $code_files
