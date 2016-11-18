# libraries
library(synapseClient)
library(VariantAnnotation)

# login to synapse
synapseLogin()

# read aggregated clinical data
genieClinData = read.delim(synGet("syn7392892")@filePath)
# read aggregated MAF file
genieMutData = read.delim(synGet("syn5571527")@filePath,header=TRUE,sep="\t",quote="",comment="",skip = 1)

# all inclusive list of samples should be from the clinical data file
# therefore factor levels for Tumor_Sample_Barcode in the MAF should be set to that of SAMPLE_ID of clinical data table 
# check that no samples are listed in the MAF that are not listed in the clinical data file
# reversing the order of the inputs would tell you which samples are submitted that have no entries (no mutations) in the MAF
if (length(setdiff(levels(genieMutData$Tumor_Sample_Barcode),levels(genieClinData$SAMPLE_ID)))==0) {
  genieMutData$Tumor_Sample_Barcode = factor(genieMutData$Tumor_Sample_Barcode, levels=levels(genieClinData$SAMPLE_ID))
}

# records with count data that preclude a VAF estimate - set VAF to 100% (1/1, alt/depth)
noVAF.idx = which((genieMutData$t_depth==0)|is.na(genieMutData$t_depth))
genieMutData$t_alt_count_num = as.numeric(levels(genieMutData$t_alt_count))[genieMutData$t_alt_count]
genieMutData$t_alt_count_num[noVAF.idx] = 1
genieMutData$t_depth[noVAF.idx] = 1

# get VRanges for all variants called in the MAF
mafVR = VRanges(seqnames=Rle(paste0("chr",genieMutData$Chromosome)),ranges=IRanges(start=genieMutData$Start_Position,end=genieMutData$End_Position),ref=genieMutData$Reference_Allele,alt=genieMutData$Tumor_Seq_Allele2,altDepth=genieMutData$t_alt_count_num,totalDepth=genieMutData$t_depth,sampleNames=genieMutData$Tumor_Sample_Barcode)
seqlevels(mafVR) = sort(seqlevels(mafVR))

# precompute
vaf = altDepth(mafVR)/totalDepth(mafVR)
ord = order(mafVR)

# start with empty table
tbl = genieMutData[1,c("Center","Tumor_Sample_Barcode","Hugo_Symbol","HGVSp_Short","Variant_Classification","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","t_alt_count_num","t_depth")]
tbl = tbl[-1,]

# check for potential variants that may need to be evaluated for merge (cis/trans)
t = Sys.time()
for (i in 1:length(genieClinData$SAMPLE_ID)) {
  # get sample indices (in order from pre sort above)
  idx = ord[which(genieMutData$Tumor_Sample_Barcode[ord]==genieClinData$SAMPLE_ID[i])]
  # get length of idx
  l = length(idx)

  # if sample has more than one variant
  if (l>1) {
    # get differences in BPs of variant sites
    dBP = distance(mafVR[idx[1:(l-1)]],mafVR[idx[2:(l)]])
    # get difference in VAFs of variants
    dVAF = abs(diff(vaf[idx]))
    
    # potential matches - criteria of difference in BPs between of > 0 & < 6 bps difference, < 5% VAF difference
    pm = which((dBP>0) & (dBP<6) & (dVAF<.05))
    for (m in pm) {
      # calc difference in codon number
      codonDiff = abs(diff(as.numeric(sapply(strsplit(as.character(genieMutData$Protein_position[c(idx[m],idx[m+1])]),split="/"),"[",1))))
      if (is.na(codonDiff)|(codonDiff==1)) {
        tbl = rbind(tbl,genieMutData[c(idx[m],idx[m+1]),c("Center","Tumor_Sample_Barcode","Hugo_Symbol","HGVSp_Short","Variant_Classification","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","t_alt_count_num","t_depth")])
      }
    }
  }

  # time ticker per 100 cases
  if ((i %% 100)==0) {
    print(i)
    t[2] = Sys.time()
    print(t[2]-t[1])
    t[1] = t[2]
  }
}

write.table(tbl,file="~/Desktop/mergeCheck.tsv",sep="\t",row.names=FALSE)
