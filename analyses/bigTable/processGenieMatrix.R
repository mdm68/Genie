# libraries
library(synapseClient)
library(VariantAnnotation)

# read aggregated clinical data
genieClinData = read.delim(synGet("syn7392892")@filePath)
# read aggregated BED file data
genieBedData = read.delim(synGet("syn7444851")@filePath,header=FALSE)
colnames(genieBedData) = c("Chromosome","Start_Position","End_Position","Hugo_Symbol","SEQ_ASSAY_ID","Feature_Type")
# read aggregated MAF file
genieMutData = read.delim(synGet("syn5571527")@filePath,header=TRUE,sep="\t",quote="",comment="",skip = 1)

# records with count data that preclude a VAF estimate - set VAF to 100% (1/1, alt/depth)
noVAF.idx = which((genieMutData$t_depth==0)|is.na(genieMutData$t_depth))
genieMutData$t_alt_count_num = as.numeric(levels(genieMutData$t_alt_count))[genieMutData$t_alt_count]
genieMutData$t_alt_count_num[noVAF.idx] = 1
genieMutData$t_depth[noVAF.idx] = 1

# these are the set of samples with defined associated bed files, it is assumed the each sample only has one record in the clinical data
# this is what is expected, we could check for this if needed 
# this filter goes away once all samples have a SEQ_ASSAY_ID with a bed file
samples.idx = which(genieClinData$SEQ_ASSAY_ID %in% levels(genieBedData$SEQ_ASSAY_ID))

# all inclusive list of samples should be from the clinical data file
# therefore factor levels for Tumor_Sample_Barcode in the MAF should be set to that of SAMPLE_ID of clinical data table 
# check that no samples are listed in the MAF that are not listed in the clinical data file
# reversing the order of the inputs would tell you which samples are submitted that have no entries (no mutations) in the MAF
if (length(setdiff(levels(genieMutData$Tumor_Sample_Barcode),levels(genieClinData$SAMPLE_ID)))==0) {
  genieMutData$Tumor_Sample_Barcode = factor(genieMutData$Tumor_Sample_Barcode, levels=levels(genieClinData$SAMPLE_ID))
}

# get VRanges for the unique variants
uniqVar = unique(genieMutData[genieMutData$Tumor_Sample_Barcode %in% genieClinData$SAMPLE_ID[samples.idx],c("Chromosome","Hugo_Symbol","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2")])
uniqVar$mutVR = VRanges(seqnames=Rle(paste0("chr",uniqVar$Chromosome)),ranges=IRanges(start=uniqVar$Start_Position,end=uniqVar$End_Position),ref=uniqVar$Reference_Allele,alt=uniqVar$Tumor_Seq_Allele2)
seqlevels(uniqVar$mutVR) = sort(seqlevels(uniqVar$mutVR))

# get VRanges for all variants called in the MAF
mafVR = VRanges(seqnames=Rle(paste0("chr",genieMutData$Chromosome)),ranges=IRanges(start=genieMutData$Start_Position,end=genieMutData$End_Position),ref=genieMutData$Reference_Allele,alt=genieMutData$Tumor_Seq_Allele2,altDepth=genieMutData$t_alt_count_num,totalDepth=genieMutData$t_depth,sampleNames=genieMutData$Tumor_Sample_Barcode)
seqlevels(mafVR) = sort(seqlevels(mafVR))

# get GRanges for each bed
bedGR = split(genieBedData, factor(genieBedData$SEQ_ASSAY_ID))
bedGR = lapply(bedGR, function(x){GRanges(seqnames=Rle(paste0("chr",x$chr)),ranges=IRanges(start=x$start,end=x$end))})

# get boolean vectors for unique variants across each panel
uniqVar[,names(bedGR)] = lapply(bedGR, function(x){ uniqVar$mutVR %over% x })

# some precomputations for optimizatin
sidxs <- match(mafVR, uniqVar$mutVR)
vaf = altDepth(mafVR)/totalDepth(mafVR)

print(length(samples.idx))
a = Sys.time()
for (i in 1:length(samples.idx)) {
  # get sample name
  sampleName = as.character(genieClinData$SAMPLE_ID[samples.idx[i]])
  
  # set all values to NA to start for sample, this will end up being the variants outside of SEQ_ASSAY_ID of this sample
  uniqVar[,sampleName] = NA
  
  # get called variants for sampleName in MAF
  sampleMAF.idx = which(sampleNames(mafVR)==sampleName)
  sampleCalledVariants = sidxs[sampleMAF.idx]
  
  # if sample has some variants called
  if (length(sampleCalledVariants)>0) {
    uniqVar[sampleCalledVariants,sampleName] = vaf[sampleMAF.idx]
  }
  
  # not called variants from the possibly callable variants based on sample's SEQ_ASSAY_ID - assign these to true 0
  uniqVar[setdiff(which(uniqVar[,as.character(genieClinData$SEQ_ASSAY_ID[genieClinData$SAMPLE_ID==sampleName])]),sampleCalledVariants),sampleName] = 0
  
  # explicitly assign variants outside of the sample's SEQ_ASSAY_ID to NA, not needed as sample vector was started with NAs
  # uniqVar[which(!uniqVar[,as.character(genieClinData$SEQ_ASSAY_ID[genieClinData$SAMPLE_ID==sampleName])]),sampleName] = NA
  
  # ticker per 100 cases
  if ((i %% 100)==0) {
    print(i)
    a[2] = Sys.time()
    print(a[2]-a[1])
    a[1] = a[2]
  }
}
