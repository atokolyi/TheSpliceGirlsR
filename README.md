# TheSpliceGirls
An annotation pipeline for transcript splicing events (i.e. from [leafcutter](https://davidaknowles.github.io/leafcutter/))
## Installation
Open R, install devtools, bioconductor, & TSG
```R
install.packages("devtools") # Or with conda: conda create -n tsg conda-forge::r-devtools; conda activate tsg;
install.packages("BiocManager")
devtools::install_github("atokolyi/TheSpliceGirls",upgrade="never")
```
## Usage
Load the library
```R
library(tsg)
```
Download/update the annotation cache (with desired gencode version, default is 46).
```R
tsg_update_cache(gencode=46)
```
Annotate your vector of splice events (e.g. in format "10:112426859:112427218:clu_12345_+"). Demo data is supplied.
```R
splices = readLines(system.file("extdata", "sample.txt", package = "tsg"))
annots = tsg_annotate(splices)
```
## Description of output
The above function (`tsg_annotate`) will return a dataframe with the following columns:
Name | Description
--- | ---
ID | The original identifier of the splice event
chr | Chromosome ID of the splice event, will start with "chr"
start | Splice start location (unstranded, i.e. always less than `end`)
end | Splice end location (unstranded, i.e. always greater than `start`)
clu | Cluster ID
strand | Strand (if unstranded, will be `*`)
start_stranded | Splice start location (stranded, i.e. if negative strand, start_stranded>end_stranded)
end_stranded | Splice end location (stranded, i.e. if negative strand, end_stranded<start_stranded)
gene_full_match | Boolean, if the splice event is fully contained within a gene body
gene_id | Gene ID of gene(s), (column-separated if multiple) the splice event falls within
gene_name | Gene symbol as above
gene_type | Gene type as above
splice_antisense | Boolean, true if the only matching gene is on the opposite strand of the splice event
range | Splice event coordinates in GRanges format
overlaps_exon | Boolean, if the transcript splicing event excises (fully or in part) an exon
overlaps_cds | 	Boolean, if the transcript splicing event excises (fully or in part) a CDS (a subset of above)
p5_exon_overlap	| Boolean, if the 5' splice donor matches a known exon boundary
p3_exon_overlap	| Boolean, if the 3' splice acceptor matches a known exon boundary
p5_and_p3_exon_overlap | Boolean, if both p3 and p5 are true
p5_or_p3_exon_overlap | Boolean, if either p3 or p5 are true (for a more lenient 'known' splice event definition)
unipLocTransMemb | Boolean, if the splice event excises a known transmembrane domain (from uniprot)
unipLocCytopl | Boolean, if the splice event excises a known cytoplasmic domain (from uniprot)
unipLocExtra | Boolean, if the splice event excises a known extracellular domain (from uniprot)
unipLocSignal | Boolean, if the splice event excises a known signal peptide domain (from uniprot)
unipDomain | Name of uniprot domains (comma-separated) that are excised by the splice event, blank if none
pfamDomain | Name of pfam domains (comma-separated) that are excised by the splice event, blank if none