library(panelcn.mops)
library(data.table)
library(R.utils)
args=cmdArgs()
bams <- Sys.glob(file.path(args$bamdir,"*bam")) # list of files
print("LIST OF BAMS:")
print(bams)
ref_genome <- args$ref_genome # "/srv/data2/pfeijao/ref/hg19.fa"  # path of reference; has to be a config/parameter.
out <- args$out # "CNV_calls.txt"
bed.file <- args$bedfile

# # if (file.exists(out)) {
# #     print(paste("==================== Output exists, skipping", out))
# #     quit(status=0, save='no')
# # }
# print(paste("================== Running on folder", out))

# Filter BLANK:
`%like%` <- function (x, pattern) {
  grepl(pattern, x, ignore.case=TRUE)
}
bams <- bams[!(bams %like% "blank")]

# MOPS: count windows
countWindows <- getWindows(bed.file)
# CNV genes in our assay: (update if necessary!)
selectedGenes <- c("ERBB2","MET","KRAS","EGFR","CCNE1","KIT","PTEN")

cnv.calls = NULL
refs<-list()

sample.names <- basename(bams)

# read counts for all samples:
all_samples <- countBamListInGRanges(countWindows = countWindows, bam.files = bams, read.width = 150)

for(i in 1:length(sample.names)){ #for each sample:
    print(paste("Processing sample: ",sample.names[i]," ", i,"/",length(sample.names),sep=""))

    # select test (clinical) sample from all samples:
    clinical <- all_samples
    clinical@elementMetadata@listData <- list(clinical@elementMetadata@listData[[i]])
    names(clinical@elementMetadata@listData) <- names(all_samples@elementMetadata@listData)[i]
    # others are control, remove i-th sample (the test sample)
    control <- all_samples
    control@elementMetadata@listData[[i]] <- NULL

    elementMetadata(clinical) <- cbind(elementMetadata(clinical), elementMetadata(control))

    resultlist <- runPanelcnMops(clinical, countWindows = countWindows, selectedGenes = selectedGenes)

    sampleNames <- colnames(elementMetadata(clinical))
    resulttable <- createResultTable(resultlist = resultlist, XandCB = clinical, countWindows = countWindows, selectedGenes = selectedGenes, sampleNames =  sampleNames)

    cnv.calls <- rbind(cnv.calls, resulttable[[1]])

    }

# Write output:
write.table(cnv.calls,file=out,sep="\t",quote=F,row.names=F,col.names=T)
