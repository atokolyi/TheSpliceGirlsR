# TheSpliceGirls

## Installation
Clone the repo
```bash
git clone https://github.com/atokolyi/TheSpliceGirls.git
cd TheSpliceGirls
```
Open R, install devtools, bioconductor, & TSG
```R
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install()
devtools::install()
```
## Usage
Download/update the annotation cache (with desired gencode version, default is 46).
```R
tsg_update_cache(gencode=46)
```
Annotate your vector of splice events (e.g. in format "10:112426859:112427218:clu_12345_+"). Demo data is supplied.
```R
splices = readLines(system.file("extdata", "sample.txt", package = "tsg"))
annots = tsg_annotate(splices)
```
