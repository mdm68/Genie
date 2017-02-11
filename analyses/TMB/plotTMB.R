# libraries
library(synapseClient)
library(tidyr)
library(plyr)
library(VariantAnnotation)
library(ggplot2)

# SAGE login
synapseLogin(username="",password="") 

# read aggregated clinical data
genieClinData = read.delim(synGet("syn7871792")@filePath,header=TRUE)
# read aggregated BED file data
genieBedData = read.delim(synGet("syn7871816")@filePath,header=TRUE)
# get GRanges for each bed
bedGR = split(genieBedData, factor(genieBedData$SEQ_ASSAY_ID))
bedGR = lapply(bedGR, function(x) {
  GR = GRanges(seqnames=Rle(paste0("chr",x$Chromosome)),ranges=IRanges(start=x$Start_Position,end=x$End_Position))
  seqlevels(GR) = sort(seqlevels(GR))
  return(GR)
})
# read aggregated MAF file
genieMutData = read.delim(synGet("syn7871815")@filePath,header=TRUE,sep="\t",quote="",comment="",skip=1)
# read CNA
genieCNAData = read.delim(synGet("syn7871784")@filePath,header=TRUE,sep="\t",check.names=FALSE)
# read fusion data
genieFusData = read.delim(synGet("syn7871793")@filePath,header=TRUE,sep="\t",check.names=FALSE)

# contrust basic panel statistics
panelStats = aggregate(SAMPLE_ID ~ (SEQ_ASSAY_ID+CENTER),genieClinData,length)
rownames(panelStats) = panelStats$SEQ_ASSAY_ID
panelStats = panelStats[,-which(colnames(panelStats)=="SEQ_ASSAY_ID")]
colnames(panelStats)[colnames(panelStats)=="CENTER"] = "Center"
colnames(panelStats)[colnames(panelStats)=="SAMPLE_ID"] = "Samples"
panelStats[,c("mut","cnv","fus")] = FALSE
panelStats$mut[unique(genieClinData$SEQ_ASSAY_ID[match(unique(genieMutData$Tumor_Sample_Barcode),genieClinData$SAMPLE_ID)])] = TRUE
panelStats$cnv[unique(genieClinData$SEQ_ASSAY_ID[match(unique(setdiff(colnames(genieCNAData),"HUGO_SYMBOL")),genieClinData$SAMPLE_ID)])] = TRUE
panelStats$fus[unique(genieClinData$SEQ_ASSAY_ID[match(unique(genieFusData$TUMOR_SAMPLE_BARCODE),genieClinData$SAMPLE_ID)])] = TRUE
panelStats$bps = NA
for (panelName in rownames(panelStats)) {panelStats[panelName,"bps"] = sum(width(reduce(bedGR[[panelName]])))}

# contruct the cancer type factor for graphing and sort this correctly
t = table(genieClinData$CANCER_TYPE)
genieClinData$CANCER_TYPE_GRAPH = as.character(genieClinData$CANCER_TYPE)
genieClinData$CANCER_TYPE_GRAPH[genieClinData$CANCER_TYPE %in% rownames(t)[which(t<100)]] = "Other"
genieClinData$CANCER_TYPE_GRAPH = factor(genieClinData$CANCER_TYPE_GRAPH,levels=rev(c("Other",rownames(sort(t[t>=100],decreasing=FALSE)))))
t = table(genieClinData$CENTER)
genieClinData$CENTER = factor(genieClinData$CENTER,levels=rownames(sort(t,decreasing=TRUE)))

# get all variant from MAF
genieVR = VRanges(seqnames=Rle(paste0("chr",genieMutData$Chromosome)),ranges=IRanges(start=genieMutData$Start_Position,end=genieMutData$End_Position),ref=genieMutData$Reference_Allele,alt=genieMutData$Tumor_Seq_Allele2)

# panel lists
largePanels = rownames(panelStats)[panelStats$cnv]
smallPanels = rownames(panelStats)[!panelStats$cnv]

# all non-silent mutations
t = aggregate(t_depth ~ Tumor_Sample_Barcode,genieMutData[genieMutData$Variant_Classification!="Silent",],length)
genieClinData$numberMutations = 0
genieClinData$numberMutations[match(t$Tumor_Sample_Barcode,genieClinData$SAMPLE_ID)] = t$t_depth
genieClinData$tmb = genieClinData$numberMutations/panelStats[as.character(genieClinData$SEQ_ASSAY_ID),"bps"]

# contruct the cancer type factor for graphing and sort this correctly
t = table(genieClinData$CANCER_TYPE)
genieClinData$CANCER_TYPE_GRAPH = as.character(genieClinData$CANCER_TYPE)
genieClinData$CANCER_TYPE_GRAPH[genieClinData$CANCER_TYPE %in% rownames(t)[which(t<100)]] = "Other"
genieClinData$CANCER_TYPE_GRAPH = factor(genieClinData$CANCER_TYPE_GRAPH,levels=rev(c("Other",rownames(sort(t[t>=100],decreasing=FALSE)))))
t = table(genieClinData$CENTER)
genieClinData$CENTER = factor(genieClinData$CENTER,levels=rownames(sort(t,decreasing=TRUE)))

# set SEQ_ASSAY_ID
genieClinData$SEQ_ASSAY_ID = factor(genieClinData$SEQ_ASSAY_ID,levels=rownames(panelStats)[order(panelStats$bps)])
panelStats$SEQ_ASSAY_ID = factor(panelStats$SEQ_ASSAY_ID,levels=rownames(panelStats)[order(panelStats$bps)])

# set TMB ticks
l = c(0,1,2,5,10,20,50,100,200,500)
# set df for graphing
df = genieClinData[genieClinData$SEQ_ASSAY_ID %in% largePanels,]

# sort by median tmb in the df
t = aggregate(tmb ~ CANCER_TYPE_GRAPH,df,function(x) {quantile(x,0.5)})

# get ecdf values
var = "CANCER_TYPE_GRAPH"
for (val in unique(df[,var])) {
  idx = val==df[,var]
  fn = ecdf(log10((df$tmb[idx]*1e6)+1))
  df$ecdf[idx] = fn(log10((df$tmb[idx]*1e6)+1))
}

# make plot
ggplot(df,aes(x=log10((tmb*1e6)+1),y=1-ecdf))+
  geom_vline(aes(xintercept=0),alpha=0.2)+
  geom_hline(aes(yintercept=0),alpha=0.2)+
  geom_ribbon(aes(ymin=0,ymax=1-ecdf),fill="lightblue",alpha=0.3,linetype=0)+
  geom_vline(mapping=aes(xintercept=`log10((tmb * 1e+06) + 1)`),data=aggregate(log10((tmb*1e6)+1)~CANCER_TYPE_GRAPH,df,median),color="steelblue")+
  geom_label(mapping=aes(log10(500+1),0.5,label=paste0("n=",tmb)),data=aggregate(tmb ~ CANCER_TYPE_GRAPH,df,length),size=2.5,fill="white",label.padding=unit(0.2,"lines"))+
  geom_point(mapping=aes(color=log10((tmb*1e6)+1)),alpha=1,size=0.5)+
  scale_color_gradient2(low="black",mid="grey20",high="red",midpoint=log10(2+1))+
  facet_wrap(facets=~factor(CANCER_TYPE_GRAPH,levels=t$CANCER_TYPE_GRAPH[order(t$tmb)]),nrow=1,strip.position="bottom")+
  scale_x_continuous(breaks=log10(l+1),labels=l,limits=log10(l[c(1,length(l))]+1))+
  scale_y_reverse(breaks=c(1,.5,0),limits=c(1.2,-.2))+
  labs(x="Variants per Mbs")+coord_flip()+
  theme_minimal()+theme(axis.title.y=element_text(size=12),axis.text.x=element_blank(),axis.title.x=element_blank(),strip.text=element_text(angle=90,vjust=1,size=12),legend.position="none")
