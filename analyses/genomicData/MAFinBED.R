# libraries
library(synapseClient)
library(VariantAnnotation)

# SAGE login
synapseLogin(username=,password=) # set user and password

# read aggregated clinical data
genieClinData = read.delim(synGet("syn7392892")@filePath,header=TRUE)
# read aggregated BED file data
genieBedData = read.delim(synGet("syn7444851")@filePath,header=TRUE)
# read aggregated MAF file
genieMutData = read.delim(synGet("syn5571527")@filePath,header=TRUE,sep="\t",quote="",comment="",skip = 1)

# records with count data that preclude a VAF estimate - set VAF to 100% (1/1, alt/depth)
noVAF.idx = which((genieMutData$t_depth==0)|is.na(genieMutData$t_depth))
genieMutData$t_alt_count_num = as.numeric(levels(genieMutData$t_alt_count))[genieMutData$t_alt_count]
genieMutData$t_alt_count_num[noVAF.idx] = 1
genieMutData$t_depth[noVAF.idx] = 1

# get VRanges for all variants called in the MAF
mafVR = VRanges(seqnames=Rle(paste0("chr",genieMutData$Chromosome)),ranges=IRanges(start=genieMutData$Start_Position,end=genieMutData$End_Position),ref=genieMutData$Reference_Allele,alt=genieMutData$Tumor_Seq_Allele2,altDepth=genieMutData$t_alt_count_num,totalDepth=genieMutData$t_depth,sampleNames=genieMutData$Tumor_Sample_Barcode)
seqlevels(mafVR) = sort(seqlevels(mafVR))

# get GRanges for each bed
bedGR = split(genieBedData, factor(genieBedData$SEQ_ASSAY_ID))
bedGR = lapply(bedGR, function(x) {
  GR = GRanges(seqnames=Rle(paste0("chr",x$Chromosome)),ranges=IRanges(start=x$Start_Position,end=x$End_Position))
  seqlevels(GR) = sort(seqlevels(GR))
  return(GR)
})

# all inclusive list of samples should be from the clinical data file
# therefore factor levels for Tumor_Sample_Barcode in the MAF should be set to that of SAMPLE_ID of clinical data table 
# check that no samples are listed in the MAF that are not listed in the clinical data file
# reversing the order of the inputs would tell you which samples are submitted that have no entries (no mutations) in the MAF
if (length(setdiff(levels(genieMutData$Tumor_Sample_Barcode),levels(genieClinData$SAMPLE_ID)))==0) {
  genieMutData$Tumor_Sample_Barcode = factor(genieMutData$Tumor_Sample_Barcode, levels=levels(genieClinData$SAMPLE_ID))
}

# add factor to MAF and set to FALSE, forcing the variant to bed match for TRUE below to clear filter
genieMutData$inBED = FALSE

# collect some stats in panStats and set inBED to TRUE if variant is in corresponding BED
panStats = data.frame()
for (panelName in levels(genieClinData$SEQ_ASSAY_ID)) {
  print(panelName)
  samples.idx = which(genieClinData$SEQ_ASSAY_ID==panelName)
  panStats[panelName,"samples"] = length(samples.idx)
  samples.idx = which(genieMutData$Tumor_Sample_Barcode %in% genieClinData$SAMPLE_ID[samples.idx])
  panStats[panelName,"samples with calls"] = length(unique(genieMutData$Tumor_Sample_Barcode[samples.idx]))
  panStats[panelName,"variants calls"] = length(samples.idx)
  panStats[panelName,"unique variants calls"] = length(unique(mafVR[samples.idx]))
  genieMutData$inBED[samples.idx] = mafVR[samples.idx] %over% bedGR[[panelName]]
  panStats[panelName,"variants called out of BED"] = length(samples.idx[!genieMutData$inBED[samples.idx]])
  panStats[panelName,"unique variants called out of BED"] = length(unique(mafVR[samples.idx[!genieMutData$inBED[samples.idx]]]))
}


