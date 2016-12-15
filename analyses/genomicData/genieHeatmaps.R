# libraries
library(synapseClient)
library(VariantAnnotation)

# SAGE login
synapseLogin(username=,password=) # put your name and password

# read aggregated clinical data
genieClinData = read.delim(synGet("syn7392892")@filePath,header=TRUE)
# read aggregated BED file data
genieBedData = read.delim(synGet("syn7444851")@filePath,header=TRUE)
# read aggregated MAF file
genieMutData = read.delim(synGet("syn5571527")@filePath,header=TRUE,sep="\t",quote="",comment="",skip = 1)

# load JHH Ion blacklist
chad.blacklist = read.csv("~/Desktop/chad.blacklist.csv")
chad.blacklist.VR = VRanges(seqnames=Rle(paste0("chr",chad.blacklist$Chromosome)),ranges=IRanges(start=chad.blacklist$Start_Position,end=chad.blacklist$End_Position),ref=chad.blacklist$ref,alt=chad.blacklist$alt)

# records with count data that preclude a VAF estimate - set VAF to 100% (1/1, alt/depth)
noVAF.idx = which((genieMutData$t_depth==0)|is.na(genieMutData$t_depth))
genieMutData$t_alt_count_num = as.numeric(levels(genieMutData$t_alt_count))[genieMutData$t_alt_count]
genieMutData$t_alt_count_num[noVAF.idx] = 1
genieMutData$t_depth[noVAF.idx] = 1

# all inclusive list of samples should be from the clinical data file
# therefore factor levels for Tumor_Sample_Barcode in the MAF should be set to that of SAMPLE_ID of clinical data table 
# check that no samples are listed in the MAF that are not listed in the clinical data file
# reversing the order of the inputs would tell you which samples are submitted that have no entries (no mutations) in the MAF
if (length(setdiff(levels(genieMutData$Tumor_Sample_Barcode),levels(genieClinData$SAMPLE_ID)))==0) {
  genieMutData$Tumor_Sample_Barcode = factor(genieMutData$Tumor_Sample_Barcode, levels=levels(genieClinData$SAMPLE_ID))
}

# get VRanges for all variants called in the MAF
mafVR = VRanges(seqnames=Rle(paste0("chr",genieMutData$Chromosome)),ranges=IRanges(start=genieMutData$Start_Position,end=genieMutData$End_Position),ref=genieMutData$Reference_Allele,alt=genieMutData$Tumor_Seq_Allele2,altDepth=genieMutData$t_alt_count_num,totalDepth=genieMutData$t_depth,sampleNames=genieMutData$Tumor_Sample_Barcode)
seqlevels(mafVR) = sort(seqlevels(mafVR))

# remove JHH Ion blacklist and remake mafVR
genieMutData = genieMutData[!((genieMutData$Center=="JHH")&(mafVR %in% chad.blacklist.VR)),]
mafVR = VRanges(seqnames=Rle(paste0("chr",genieMutData$Chromosome)),ranges=IRanges(start=genieMutData$Start_Position,end=genieMutData$End_Position),ref=genieMutData$Reference_Allele,alt=genieMutData$Tumor_Seq_Allele2,altDepth=genieMutData$t_alt_count_num,totalDepth=genieMutData$t_depth,sampleNames=genieMutData$Tumor_Sample_Barcode)
seqlevels(mafVR) = sort(seqlevels(mafVR))

# restrict to not "commmon" germline filter and not silent
mafFilter = (genieMutData$FILTER!="common_variant") & (genieMutData$Variant_Classification!="Silent")

# apply filter to mafVR
mafVR = mafVR[mafFilter]

# sample based on tumor type filter
samples.idx = which(genieClinData$CANCER_TYPE=="Melanoma")

# MAF indices - logical index against aggregated MAF, based on tumor type samples + filters
maf.idx = (genieMutData$Tumor_Sample_Barcode %in% genieClinData$SAMPLE_ID[samples.idx]) & mafFilter

# get VRanges for the unique variants
uniqVar = unique(genieMutData[maf.idx,c("Chromosome","Hugo_Symbol","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","HGVSp_Short","HGVSc")])
uniqVar$mutVR = VRanges(seqnames=Rle(paste0("chr",uniqVar$Chromosome)),ranges=IRanges(start=uniqVar$Start_Position,end=uniqVar$End_Position),ref=uniqVar$Reference_Allele,alt=uniqVar$Tumor_Seq_Allele2)
seqlevels(uniqVar$mutVR) = sort(seqlevels(uniqVar$mutVR))

# get GRanges for each bed
bedGR = split(genieBedData, factor(genieBedData$SEQ_ASSAY_ID))
bedGR = lapply(bedGR, function(x) {
  GR = GRanges(seqnames=Rle(paste0("chr",x$Chromosome)),ranges=IRanges(start=x$Start_Position,end=x$End_Position))
  seqlevels(GR) = sort(seqlevels(GR))
  return(GR)
})

# get boolean vectors for unique variants across each panel
uniqVar[,names(bedGR)] = lapply(bedGR, function(x){ uniqVar$mutVR %over% x })

# some precomputations for optimizatin
sidxs <- match(mafVR, uniqVar$mutVR)
vaf = altDepth(mafVR)/totalDepth(mafVR)

# construct initial variant X sample data.frame
print(length(samples.idx))
a = Sys.time()
for (i in 1:length(samples.idx)) {
  # get sample name
  sampleName = as.character(genieClinData$SAMPLE_ID[samples.idx[i]])
  
  # set all variants to NA to start for sample, this will end up being the variants outside of SEQ_ASSAY_ID of this sample
  uniqVar[,sampleName] = NA
  
  # then set all the variants that this sample's SEQ_ASSAY_ID covers to 0 (ie not called)
  uniqVar[uniqVar[,as.character(genieClinData$SEQ_ASSAY_ID[samples.idx[i]])],sampleName] = 0
  
  # then for any called variants in the MAF for this sample insert the variant allele frequency
  sampleMAF.idx = which(sampleNames(mafVR)==sampleName)
  sampleCalledVariants = sidxs[sampleMAF.idx]
  if (length(sampleCalledVariants)>0) {
    uniqVar[sampleCalledVariants,sampleName] = vaf[sampleMAF.idx]
  }
  
  # ticker per 100 cases
  if ((i %% 100)==0) {
    print(i)
    a[2] = Sys.time()
    print(a[2]-a[1])
    
    a[1] = a[2]
  }
}

# identify where sample has a call for a given unique variant but the corrresponding bed file for that sample says the variant is not covered
a = rowsum(t(1*((uniqVar[,as.character(genieClinData$SAMPLE_ID[samples.idx])]>0)&(!is.na(uniqVar[,as.character(genieClinData$SAMPLE_ID[samples.idx])])))),group=genieClinData$SEQ_ASSAY_ID[samples.idx])
# change value for the variant to TRUE (ie detectable) for variants seen in the data
for (x in as.character(unique(genieClinData$SEQ_ASSAY_ID[samples.idx]))) { uniqVar[!uniqVar[,x]&(a[x,]>0),x] = TRUE }

# repopulate the unique var X sample data.frame based on updated covering matrix
print(length(samples.idx))
a = Sys.time()
for (i in 1:length(samples.idx)) {
  # get sample name
  sampleName = as.character(genieClinData$SAMPLE_ID[samples.idx[i]])
  
  # set all variants to NA to start for sample, this will end up being the variants outside of SEQ_ASSAY_ID of this sample
  uniqVar[,sampleName] = NA
  
  # then set all the variants that this sample's SEQ_ASSAY_ID covers to 0 (ie not called)
  uniqVar[uniqVar[,as.character(genieClinData$SEQ_ASSAY_ID[samples.idx[i]])],sampleName] = 0
  
  # then for any called variants in the MAF for this sample insert the variant allele frequency
  sampleMAF.idx = which(sampleNames(mafVR)==sampleName)
  sampleCalledVariants = sidxs[sampleMAF.idx]
  if (length(sampleCalledVariants)>0) {
    uniqVar[sampleCalledVariants,sampleName] = vaf[sampleMAF.idx]
  }
  
  # ticker per 100 cases
  if ((i %% 100)==0) {
    print(i)
    a[2] = Sys.time()
    print(a[2]-a[1])
    
    a[1] = a[2]
  }
}

# count of calls (rows are site and colums are variants)
a = rowsum(t(1*((uniqVar[,as.character(genieClinData$SAMPLE_ID[samples.idx])]>0)&(!is.na(uniqVar[,as.character(genieClinData$SAMPLE_ID[samples.idx])])))),group=genieClinData$CENTER[samples.idx])
# count of true non-detects (ie coverage is there but variant is not reported)
b = rowsum(t(1*((uniqVar[,as.character(genieClinData$SAMPLE_ID[samples.idx])]==0)&(!is.na(uniqVar[,as.character(genieClinData$SAMPLE_ID[samples.idx])])))),group=genieClinData$CENTER[samples.idx])

# seen in at least 1 sample in at least 1 site
k = which(colSums(a>=1)>=1)
# site has at least 50 samples & in least 1 variant
j = which(rowSums((a[,k]+b[,k])>=50)>=1)

# get detection rate at variant level as data.frame with filters applied
d = as.data.frame(a[j,k]/(a[j,k]+b[j,k]))
colnames(d) = k
d$Center = rownames(d)
# reshape to long format for geom_tile of ggplot
d = reshape(d,direction="long",idvar="Center",varying=as.character(k),v.names="test",times=as.character(k))
colnames(d) = c("Center","idx","Detect.Rate")
d$idx = as.numeric(d$idx)
# add in gene symbol and variant descriptor -- you can customized here
d$Hugo_Symbol = uniqVar$Hugo_Symbol[d$idx]
d$Variant = paste0(uniqVar$Hugo_Symbol[d$idx]," ",uniqVar$HGVSp_Short[d$idx]," (",uniqVar$Chromosome[d$idx],":",uniqVar$Start_Position[d$idx],":",uniqVar$Reference_Allele[d$idx],":",uniqVar$Tumor_Seq_Allele2[d$idx],")")
# sort variants by coverage across site and detection rate -- this is done ultimately by sorting the levels of the factor
l = aggregate(Detect.Rate ~ (Variant + Hugo_Symbol),d,function(x) {c(min(x,na.rm=TRUE),mean(is.na(x)))},na.action=NULL)
l = l[order(-l$Detect.Rate[,2],l$Detect.Rate[,1],decreasing=TRUE),]
d$Variant = factor(d$Variant,levels=l$Variant)
# sort sites by number of variants with calls made -- this is done ultimately by sorting the levels of the factor
l = aggregate(Detect.Rate ~ Center,d,function(x) {mean((x>0)&(!is.na(x)))},na.action=NULL)
l = l[order(l$Detect.Rate,decreasing=FALSE),]
d$Center = factor(d$Center,levels=l$Center)
# plot
ggplot(d, aes(Hugo_Symbol,Center)) + geom_tile(aes(fill=Detect.Rate),colour="white") + scale_fill_gradientn(colours=c("grey95","lightblue1","red"),values=c(0,0.0000000001,1),na.value="grey50") + labs(x ="",y="") + theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5,size=8),axis.text.y=element_text(size=8)) + scale_x_discrete(expand=c(0, 0),position="top") + scale_y_discrete(expand=c(0, 0))

# get detection rate aggregated to gene level as data.frame with filters applied
d = as.data.frame(t(apply(a[j,k],1,function(x) {tapply(x,as.character(uniqVar$Hugo_Symbol[k]),sum)})/apply(a[j,k]+b[j,k],1,function(x) {tapply(x,as.character(uniqVar$Hugo_Symbol[k]),max)})))
# filter to detect rate of at least 1% in at least 1 site
d = d[,which(colSums(d>=.01)>=1)]
l = colnames(d)
d$Center = rownames(d)
# reshape to long format for geom_tile of ggplot
d = reshape(d,direction="long",idvar="Center",varying=l,v.names="test",times=l)
colnames(d) = c("Center","Gene","Detect.Rate")
# sort genes by coverage across site and detection rate -- this is done ultimately by sorting the levels of the factor
l = aggregate(Detect.Rate ~ Gene,d,function(x) {c(min(x,na.rm=TRUE),mean(is.na(x)))},na.action=NULL)
l = l[order(-l$Detect.Rate[,2],l$Detect.Rate[,1],decreasing=TRUE),]
d$Gene = factor(d$Gene,levels=l$Gene)
# sort sites by number of genes with calls made -- this is done ultimately by sorting the levels of the factor
l = aggregate(Detect.Rate ~ Center,d,function(x) {mean((x>0)&(!is.na(x)))},na.action=NULL)
l = l[order(l$Detect.Rate,decreasing=FALSE),]
d$Center = factor(d$Center,levels=l$Center)
# plot
ggplot(d, aes(Gene,Center)) + geom_tile(aes(fill=Detect.Rate),colour="white") + scale_fill_gradientn(colours=c("grey95","lightblue1","red"),values=c(0,0.000000000001,1),na.value="grey50") + labs(x ="",y="") + theme(axis.text.x=element_text(angle=45,hjust=0,vjust=0.5,size=8),axis.text.y=element_text(size=8)) + scale_x_discrete(expand=c(0, 0),position="top") + scale_y_discrete(expand=c(0, 0))
