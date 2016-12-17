# libraries
library(synapseClient)
library(VariantAnnotation)
library(ComplexHeatmap)
library(circlize)

# SAGE login
synapseLogin(username=,password=) # you could do this with prompt or coded, but will need to enter something here

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

# get GRanges for each bed
bedGR = split(genieBedData, factor(genieBedData$SEQ_ASSAY_ID))
bedGR = lapply(bedGR, function(x) {
  GR = GRanges(seqnames=Rle(paste0("chr",x$Chromosome)),ranges=IRanges(start=x$Start_Position,end=x$End_Position))
  seqlevels(GR) = sort(seqlevels(GR))
  return(GR)
})

# restrict to not "commmon" germline filter and not silent
mafFilter = (genieMutData$FILTER!="common_variant") & (genieMutData$Variant_Classification!="Silent")

# apply filter to mafVR
mafVR = mafVR[mafFilter]

# sample based on tumor type filter
samples.idx = which(genieClinData$CANCER_TYPE=="Colorectal Cancer")

# MAF indices - logical index against aggregated MAF, based on tumor type samples + filters
maf.idx = (genieMutData$Tumor_Sample_Barcode %in% genieClinData$SAMPLE_ID[samples.idx]) & mafFilter

# get VRanges for the unique variants
uniqVar = unique(genieMutData[maf.idx,c("Chromosome","Hugo_Symbol","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","HGVSp_Short","HGVSc")])
uniqVar$mutVR = VRanges(seqnames=Rle(paste0("chr",uniqVar$Chromosome)),ranges=IRanges(start=uniqVar$Start_Position,end=uniqVar$End_Position),ref=uniqVar$Reference_Allele,alt=uniqVar$Tumor_Seq_Allele2)
seqlevels(uniqVar$mutVR) = sort(seqlevels(uniqVar$mutVR))

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

# VARIANT FILTER - seen in at least 1 sample in at least 1 site
k = which(colSums(a>=5)>=1)
# SITE FILTER - site has at least 50 samples & in least 1 variant
j = which(rowSums((a[,k]+b[,k])>=50)>=1)

# get detection rate at VARIANT level as data.frame with filters applied
d = as.data.frame(a[j,k]/(a[j,k]+b[j,k]))
var.idx = k[colnames(d)]
# label HUGO + HGSVp + genomic
# var.labels = paste0(uniqVar$Hugo_Symbol[var.idx]," ",uniqVar$HGVSp_Short[var.idx]," (",uniqVar$Chromosome[var.idx],":",uniqVar$Start_Position[var.idx],":",uniqVar$Reference_Allele[var.idx],":",uniqVar$Tumor_Seq_Allele2[var.idx],")")
# label HUGO + HGSVp + genomic
var.labels = paste0(uniqVar$Hugo_Symbol[var.idx]," ",uniqVar$HGVSp_Short[var.idx])
colnames(d) = var.labels

# column order
m = apply(d,2,function(x) {c(min(x,na.rm=TRUE),mean(is.na(x)))})
gOrd = order(-m[2,],m[1,],decreasing=TRUE)
# row order
n = data.frame(withHits=as.numeric(rowSums((d>0)&!is.na(d))))
n$hitRate = n$withHits/rowSums(!is.na(d))
sOrd = order(n$withHits,decreasing=TRUE)
# heatmap annotations
h.r = rowAnnotation(siteHitCount=row_anno_barplot(as.numeric(n$withHits[sOrd]),axis=TRUE,axis_side="top",border=FALSE,baseline=0),width=unit(4,"cm"))
h.c = columnAnnotation(minDetectRate=column_anno_barplot(colSums(d[,gOrd],na.rm=TRUE)/colSums(!is.na(d[,gOrd])),axis=TRUE,border=FALSE,baseline=0),height=unit(4,"cm"))
# heatmap
h.m = Heatmap(d[sOrd,gOrd],col=colorRamp2(c(0,.000000000001,.2),c("grey95","lightblue1","red")),na_col="grey50",top_annotation=h.c,row_names_side="left",column_names_side="bottom",cluster_rows=FALSE,cluster_columns=FALSE,heatmap_legend_param=list(title="Detect Rate",color_bar="continuous",legend_direction="vertical",title_gp=gpar(fontsize=12),title_position="leftcenter"),column_title="Gene",column_title_side="bottom",row_names_max_width=unit(4,"cm"),column_names_max_height=unit(8,"cm"),row_names_gp=gpar(fontsize=10),column_names_gp=gpar(fontsize=10))
#draw
draw(h.m+h.r,heatmap_legend_side="left")
  
# get detection rate aggregated to GENE level as data.frame with filters applied
d = as.data.frame(t(apply(a[j,k],1,function(x) {tapply(x,as.character(uniqVar$Hugo_Symbol[k]),sum)})/apply(a[j,k]+b[j,k],1,function(x) {tapply(x,as.character(uniqVar$Hugo_Symbol[k]),max)})))
# GENE FILTER - detected in at least X % for at least 2 site
d = d[,which(colSums(((d>=0.05)&!is.na(d)))>=2)]

# column order
m = apply(d,2,function(x) {c(min(x,na.rm=TRUE),mean(is.na(x)))})
gOrd = order(-m[2,],m[1,],decreasing=TRUE)
# row order
n = data.frame(withHits=as.numeric(rowSums((d>0)&!is.na(d))))
n$hitRate = n$withHits/rowSums(!is.na(d))
sOrd = order(n$withHits,decreasing=TRUE)
# heatmap annotations
h.r = rowAnnotation(siteHitCount=row_anno_barplot(as.numeric(n$withHits[sOrd]),axis=TRUE,axis_side="top",border=FALSE,baseline=0),width=unit(4,"cm"))
h.c = columnAnnotation(minDetectRate=column_anno_barplot(colSums(d[,gOrd],na.rm=TRUE)/colSums(!is.na(d[,gOrd])),axis=TRUE,border=FALSE,baseline=0),height=unit(4,"cm"))
# heatmap
h.m = Heatmap(d[sOrd,gOrd],col=colorRamp2(c(0,.000000000001,.2),c("grey95","lightblue1","red")),na_col="grey50",top_annotation=h.c,row_names_side="left",column_names_side="bottom",cluster_rows=FALSE,cluster_columns=FALSE,heatmap_legend_param=list(title="Detect Rate",color_bar="continuous",legend_direction="vertical",title_gp=gpar(fontsize=12),title_position="leftcenter"),column_title="Gene",column_title_side="bottom",row_names_max_width=unit(4,"cm"),column_names_max_height=unit(8,"cm"),row_names_gp=gpar(fontsize=10),column_names_gp=gpar(fontsize=10))
#draw
draw(h.m+h.r,heatmap_legend_side="left")
