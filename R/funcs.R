#' Update TheSpliceGirls cache
#'
#' This function updates the cache for TheSpliceGirls.
#'
#' @param gencode Desired gencode version (default: 46)
#' @export
tsg_update_cache = function(gencode="46") {
    options(timeout=6000)
    cache_dir <- tools::R_user_dir("thesplicegirls", "cache")
    if(!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive=TRUE)
    } else {
        message(paste0("Last updated: ",file.mtime(file.path(cache_dir,"ucscGenePfam.txt.gz"))),appendLF=F)
    }
    host = "http://hgdownload.soe.ucsc.edu/"
    files = paste0(host,c("goldenPath/hg38/database/ucscGenePfam.txt.gz","gbdb/hg38/uniprot/unipDomain.bb","gbdb/hg38/uniprot/unipLocCytopl.bb","gbdb/hg38/uniprot/unipLocExtra.bb","gbdb/hg38/uniprot/unipLocSignal.bb",
                         "gbdb/hg38/uniprot/unipLocTransMemb.bb"))
    message("Downloading:",appendLF=F)
    for (f in files) {
        message(paste0("\t",basename(f)),appendLF=F)
        download.file(f, destfile = file.path(cache_dir,basename(f)), mode="wb")
    }    
    options(timeout=1000)
    if (file.exists(file.path(cache_dir,paste0("gencode.v",gencode,".annotation.gff3.gz")))) {
        message(paste0("\tGencode is up to date."),appendLF=F)
    } else {
        gff = paste0("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",gencode,"/gencode.v",gencode,".annotation.gff3.gz")
        message(paste0("\tGencode (v",gencode,")..(~3 mins).."),appendLF=F)
        download.file(gff, destfile = file.path(cache_dir,basename(gff)), mode="wb")
    }
    # Download and run bigbedtobed
    #https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
    message("Converting:",appendLF=F)
    bin = "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed"
    download.file(bin, destfile = file.path(cache_dir,basename(bin)), mode="wb")
    system(paste0("chmod +x ",file.path(cache_dir,basename(bin))),intern=TRUE)
    for (f in files[2:length(files)]) {
        f = stringr::str_replace_all(basename(f),".bb","")
        message(paste0("\t",f,".bb to .bed"),appendLF=F)
        system(paste0(file.path(cache_dir,basename(bin))," ",file.path(cache_dir,paste0(f,".bb"))," ",file.path(cache_dir,paste0(f,".bed"))),intern=TRUE)
    }
}

#' Run TheSpliceGirls
#'
#' This function runs TheSpliceGirls annotation.
#'
#' @param splices A vector of leafcutter transcript splicing IDs, i.e. "11:2384682:2390412:clu_12345_+"
#' @return A dataframe containing the annotations
#' @export
tsg_annotate = function(splices) {
    cache_dir <- tools::R_user_dir("thesplicegirls", "cache")
    message(paste0("Annotating with cache version: ",file.mtime(file.path(cache_dir,"ucscGenePfam.txt.gz"))),appendLF=F)
    valid_chr = c("chrX","chrY",paste0("chr",1:22))

    message(paste0("\tAnnotating splice site information.."),appendLF=F)
    splices = cbind.data.frame(id=splices)
    splices = cbind.data.frame(splices,t(as.data.frame(stringr::str_split(splices$id,":"))))
    rownames(splices) = NULL
    colnames(splices) = c("id","chr","start","end","clu")
    splices$strand = t(as.data.frame(stringr::str_split(splices$clu,"_")))[,3]
    if (!all(grepl("chr",splices$chr))) {
        splices$chr = paste0("chr",splices$chr)
    }
    stranded = all(splices$strand %in% c("+","-"))
    if (!stranded) {
        splices[which(!splices$strand %in% c("+","-")),]$strand = "*"
    }
    splices$start_stranded = splices$start
    splices$end_stranded = splices$end
    if (stranded) {
        splices[which(splices$strand=="-"),]$start_stranded = splices[which(splices$strand=="-"),]$end
        splices[which(splices$strand=="-"),]$end_stranded = splices[which(splices$strand=="-"),]$start
    }

    message(paste0("\tLoading Gencode GFF.."),appendLF=F)
    gencode = 46
    gff = paste0("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",gencode,"/gencode.v",gencode,".annotation.gff3.gz")
    gff = rtracklayer::readGFF(file.path(cache_dir,basename(gff)))
    gff = as.data.frame(gff)
    gff = gff[which(gff$seqid %in% valid_chr),]
    gff = gff[order(gff$gene_type=="protein_coding",decreasing = T),]
    gff_genes = gff[which(gff$type=="gene"),]
    gff_exons = gff[which(gff$type=="exon"),]
    gff_cds = gff[which(gff$type=="CDS"),]
    rownames(gff_genes) = gff_genes$gene_id
    genes_granges_us = as(paste0(gff_genes$seqid,":",gff_genes$start,"-",gff_genes$end,":","*"),"GRanges")
    GenomicRanges::values(genes_granges_us) = data.frame(gene_id=gff_genes$gene_id,gene_name=gff_genes$gene_name)

    message(paste0("\tAnnotating with overlapping genes.."),appendLF=F)
    # Match splice to genes fully containing the excised region
    start_ranges = as(paste0(splices$chr,":",splices$start,":*"),"GRanges")
    ov_start=as.data.frame(GenomicRanges::findOverlaps(start_ranges,genes_granges_us))
    splices$start_match = FALSE
    for (i in 1:nrow(splices)) {
        ge = gff_genes[ov_start[which(ov_start$queryHits==i),]$subjectHits,]$gene_id
        if (length(ge)>0) {
            start_matches = c()
            for (g in ge) {
                start_matches=c(start_matches,g)
            }
            splices$start_matches[i] = paste0(start_matches,collapse=",")
            splices$start_match[i] = TRUE
        }
    }
    end_ranges = as(paste0(splices$chr,":",splices$end,":*"),"GRanges")
    ov_end=as.data.frame(GenomicRanges::findOverlaps(end_ranges,genes_granges_us))
    splices$end_match = FALSE
    for (i in 1:nrow(splices)) {
        ge = gff_genes[ov_end[which(ov_end$queryHits==i),]$subjectHits,]$gene_id   
        if (length(ge)>0) { 
            end_matches=c()
            for (g in ge) {
                end_matches=c(end_matches,g)
            }
            splices$end_matches[i] = paste0(end_matches,collapse=",")
            splices$end_match[i] = TRUE
        }
       
    }
    splices$full_match = FALSE
    splices$both_genes = ""
    for (i in 1:nrow(splices)) {
        # If both have a value
        if (splices$start_match[i] & splices$end_match[i]) {
            start_matches = unlist(stringr::str_split(splices$start_matches[i],","))
            end_matches = unlist(stringr::str_split(splices$end_matches[i],","))
            splices$full_match[i] = any(start_matches %in% end_matches)
            splices$both_genes[i] = paste0(start_matches[which(start_matches %in% end_matches)],collapse=",")
        }
    }
    # Merge above to create final gene annotation for splice events, defaulting to sense-strand gene if there's a hit
    splices$merge_id = ""
    splices$merge_name = ""
    splices$merge_type = ""
    splices$merge_OS = FALSE
    for (i in 1:nrow(splices)) {
        # $full_match if both start and end within, else partial  
        ## Order mergeids and supp by protein coding first, then others? (as other usually lncRNA)
        if (splices$full_match[i]) {
            mergeids = unlist(stringr::str_split(splices[i,]$both_genes,","))
            types = gff_genes[mergeids,]$gene_type
            
            #mergeids = mergeids[order(types=="protein_coding",decreasing = TRUE)]
            #types = types[order(types=="protein_coding",decreasing = TRUE)]
            
            samestrand = c()
            for (m in mergeids) {
                samestrand = c(samestrand,gff_genes[m,]$strand==splices[i,]$strand)
            }
            if (sum(samestrand)>0) {
                splices$merge_id[i] = paste0(mergeids[samestrand],collapse=",")
                splices$merge_name[i] = paste0(gff_genes[mergeids[samestrand],]$gene_name,collapse=",")
                splices$merge_type[i] = paste0(gff_genes[mergeids[samestrand],]$gene_type,collapse=",")
            } else {
                splices$merge_id[i] = paste0(mergeids,collapse=",")
                splices$merge_name[i] = paste0(gff_genes[mergeids,]$gene_name,collapse=",")
                splices$merge_type[i] = paste0(gff_genes[mergeids,]$gene_type,collapse=",")
                if (splices$strand[i] %in% c("+","-")) {
                    splices$merge_OS[i] = TRUE
                }
            }
        } else {
            if (splices$start_match[i] | splices$end_match[i]) {
                mergeids = unlist(stringr::str_split(c(splices[i,]$start_matches,splices[i,]$end_matches),","))
                mergeids = mergeids[which(mergeids!="")]
                types = gff_genes[mergeids,]$gene_type
                
                #mergeids = mergeids[order(types=="protein_coding",decreasing = TRUE)]
                #types = types[order(types=="protein_coding",decreasing = TRUE)]
                samestrand = c()
                for (m in mergeids) {
                    samestrand = c(samestrand,gff_genes[m,]$strand==splices[i,]$strand)
                }
                if (sum(samestrand)>0) {
                    splices$merge_id[i] = paste0(mergeids[samestrand],collapse=",")
                    splices$merge_name[i] = paste0(gff_genes[mergeids[samestrand],]$gene_name,collapse=",")
                    splices$merge_type[i] = paste0(gff_genes[mergeids[samestrand],]$gene_type,collapse=",")
                } else {
                    splices$merge_id[i] = paste0(mergeids,collapse=",")
                    splices$merge_name[i] = paste0(gff_genes[mergeids,]$gene_name,collapse=",")
                    splices$merge_type[i] = paste0(gff_genes[mergeids,]$gene_type,collapse=",")
                    if (splices$strand[i] %in% c("+","-")) {
                        splices$merge_OS[i] = TRUE
                    }
                }   
            } else {
                    splices$merge_id[i] = ""
                    splices$merge_name[i] = ""
                    splices$merge_type[i] = ""
                }
        }
        
    }

    message(paste0("\tAnnotating splice site overlaps.."),appendLF=F)
    splices = splices[,c(colnames(splices)[1:8],"full_match","merge_id","merge_name","merge_type","merge_OS")]
    colnames(splices)[9:13] = c("gene_full_match","gene_id","gene_name","gene_type","splice_antisense")
    exons_granges = as(paste0(gff_exons$seqid,":",gff_exons$start,"-",gff_exons$end,":",gff_exons$strand),"GRanges")
    splices$range = paste0(splices$chr,":",splices$start,"-",splices$end,":",splices$strand)
    splices$overlaps_exon = GenomicRanges::countOverlaps(as(splices$range,"GRanges"),exons_granges,minoverlap=2,maxgap=-1)>0
    cds_granges = as(paste0(gff_cds$seqid,":",gff_cds$start,"-",gff_cds$end,":",gff_cds$strand),"GRanges")
    splices$overlaps_cds = GenomicRanges::countOverlaps(as(splices$range,"GRanges"),cds_granges,minoverlap=2,maxgap=-1)>0
    
    gff_exons$start_stranded = gff_exons$start
    gff_exons[which(gff_exons$strand=="-"),]$start_stranded = gff_exons[which(gff_exons$strand=="-"),]$end
    gff_exons$end_stranded = gff_exons$end
    gff_exons[which(gff_exons$strand=="-"),]$end_stranded = gff_exons[which(gff_exons$strand=="-"),]$start
    gff_exons$chrstart = paste0(gff_exons$seqid,":",gff_exons$start_stranded)
    gff_exons$chrend = paste0(gff_exons$seqid,":",gff_exons$end_stranded)
    
    splices$start_stranded = as.numeric(splices$start_stranded)
    splices$end_stranded = as.numeric(splices$end_stranded)
    splices$p5_exon_overlap = paste0(splices$chr,":",splices$start_stranded) %in% gff_exons$chrend
    splices$p3_exon_overlap = paste0(splices$chr,":",splices$end_stranded) %in% gff_exons$chrstart
    splices$p5_and_p3_exon_overlap = splices$p5_exon & splices$p3_exon
    splices$p5_or_p3_exon_overlap = splices$p5_exon | splices$p3_exon

    message(paste0("\tLoading protein domain localizations.."),appendLF=F)
    unipLocTransMemb = read.table(file.path(cache_dir,'unipLocTransMemb.bed'), sep="\t", quote="")
    unipLocTransMemb = unipLocTransMemb[which(unipLocTransMemb$V1 %in% valid_chr),]
    unipLocTransMemb = as(paste0(unipLocTransMemb$V1,":",unipLocTransMemb$V2,"-",unipLocTransMemb$V3,":",unipLocTransMemb$V6),"GRanges")
    unipLocTransMemb = GenomicRanges::reduce(unipLocTransMemb)
    #values(unipLocTransMemb) = data.frame(type="Transmembrane")
    unipLocCytopl = read.table(file.path(cache_dir,'unipLocCytopl.bed'), sep="\t", quote="")
    unipLocCytopl = unipLocCytopl[which(unipLocCytopl$V1 %in% valid_chr),]
    unipLocCytopl = as(paste0(unipLocCytopl$V1,":",unipLocCytopl$V2,"-",unipLocCytopl$V3,":",unipLocCytopl$V6),"GRanges")
    unipLocCytopl = GenomicRanges::reduce(unipLocCytopl)
    #values(unipLocCytopl) = data.frame(type="Cytoplasm")
    unipLocExtra = read.table(file.path(cache_dir,'unipLocExtra.bed'), sep="\t", quote="")
    unipLocExtra = unipLocExtra[which(unipLocExtra$V1 %in% valid_chr),]
    unipLocExtra = as(paste0(unipLocExtra$V1,":",unipLocExtra$V2,"-",unipLocExtra$V3,":",unipLocExtra$V6),"GRanges")
    unipLocExtra = GenomicRanges::reduce(unipLocExtra)
    #values(unipLocExtra) = data.frame(type="Extracellular")
    unipLocSignal = read.table(file.path(cache_dir,'unipLocSignal.bed'), sep="\t", quote="")
    unipLocSignal = unipLocSignal[which(unipLocSignal$V1 %in% valid_chr),]
    unipLocSignal = as(paste0(unipLocSignal$V1,":",unipLocSignal$V2,"-",unipLocSignal$V3,":",unipLocSignal$V6),"GRanges")
    unipLocSignal = GenomicRanges::reduce(unipLocSignal)
    #values(unipLocSignal) = data.frame(type="SignalP")

    message(paste0("\tAnnotating protein domain localizations.."),appendLF=F)
    splices$unipLocTransMemb = GenomicRanges::countOverlaps(as(splices$range,"GRanges"),unipLocTransMemb,minoverlap=2,maxgap=-1)>0
    splices[which(!splices$overlaps_cds),]$unipLocTransMemb=FALSE
    splices$unipLocCytopl = GenomicRanges::countOverlaps(as(splices$range,"GRanges"),unipLocCytopl,minoverlap=2,maxgap=-1)>0
    splices[which(!splices$overlaps_cds),]$unipLocCytopl=FALSE
    splices$unipLocExtra = GenomicRanges::countOverlaps(as(splices$range,"GRanges"),unipLocExtra,minoverlap=2,maxgap=-1)>0
    splices[which(!splices$overlaps_cds),]$unipLocExtra=FALSE
    splices$unipLocSignal = GenomicRanges::countOverlaps(as(splices$range,"GRanges"),unipLocSignal,minoverlap=2,maxgap=-1)>0
    splices[which(!splices$overlaps_cds),]$unipLocSignal=FALSE

    message(paste0("\tLoading protein domains"),appendLF=F)
    unipDomainT = read.table(file.path(cache_dir,'unipDomain.bed'), sep="\t", quote="")
    unipDomainT = unipDomainT[which(unipDomainT$V1 %in% valid_chr),]
    unipDomain = as(paste0(unipDomainT$V1,":",unipDomainT$V2,"-",unipDomainT$V3,":",unipDomainT$V6),"GRanges")
    GenomicRanges::values(unipDomain) = data.frame(type=unipDomainT$V27)
    ucscGenePfamT = read.table(file.path(cache_dir,'ucscGenePfam.txt.gz'), sep="\t", quote="")
    ucscGenePfamT = ucscGenePfamT[which(ucscGenePfamT$V2 %in% valid_chr),]
    ucscGenePfam = as(paste0(ucscGenePfamT$V2,":",ucscGenePfamT$V3,"-",ucscGenePfamT$V4,":",ucscGenePfamT$V7),"GRanges")
    GenomicRanges::values(ucscGenePfam) = data.frame(type=ucscGenePfamT$V5)

    message(paste0("\tAnnotating protein domains"),appendLF=F)
    hits=as.data.frame(GenomicRanges::findOverlaps(as(splices$range,"GRanges"),unipDomain,minoverlap=2,maxgap=-1))
    doms = c()
    for (i in 1:length(splices$range)) {
        sub = hits[which(hits$queryHits==i),]
        if (dim(sub)[1]>0) {
            doms = c(doms,paste0(sort(unique(unipDomain$type[sub$subjectHits])),collapse = ","))
        } else {
            doms = c(doms,"")
        }
    }
    splices$unipDomain = doms
    splices[which(!splices$overlaps_cds),]$unipDomain=""
    
    hits=as.data.frame(GenomicRanges::findOverlaps(as(splices$range,"GRanges"),ucscGenePfam,minoverlap=2,maxgap=-1))
    doms = c()
    for (i in 1:length(splices$range)) {
        sub = hits[which(hits$queryHits==i),]
        if (dim(sub)[1]>0) {
            doms = c(doms,paste0(sort(unique(ucscGenePfam$type[sub$subjectHits])),collapse = ","))
        } else {
            doms = c(doms,"")
        }
    }
    splices$pfamDomain = doms
    splices[which(!splices$overlaps_cds),]$pfamDomain=""

    # Change types to numeric
    splices$start = as.integer(splices$start)
    splices$end = as.integer(splices$end)
    splices$start_stranded = as.integer(splices$start_stranded)
    splices$end_stranded = as.integer(splices$end_stranded)
    
    message(paste0("Done!"),appendLF=F)
    return(splices)
}