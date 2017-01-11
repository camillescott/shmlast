# quick riff of the JOSS pandoc command for testing purposes
pandoc --template latex.template --filter pandoc-citeproc  -V geometry:margin=1in -r markdown+yaml_metadata_block paper.md -o paper.pdf
