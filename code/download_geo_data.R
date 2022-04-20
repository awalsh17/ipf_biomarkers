# Download the data

# 20 Gb file - need to have space for this

dir.create("data")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893_IPF_metadata.csv.gz",
              "data/GSE135893_IPF_metadata.csv.gz")

download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893_ILD_annotated_fullsize.rds.gz",
              "data/GSE135893_ILD_annotated_fullsize.rds.gz")
# or use wget from terminal

# unzip the rds
unzip("data/GSE135893_ILD_annotated_fullsize.rds.gz")
# gzip -d file.gz