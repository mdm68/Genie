# libraries
library(synapseClient)

# SAGE login
synapseLogin(username=,password=) # you could do this with prompt or coded, but will need to enter something here

# get bed and clinical data from SAGE
genieClinData = read.delim(synGet("syn7871792")@filePath,header=TRUE)
genieBedData = read.delim(synGet("syn7871816")@filePath,header=TRUE)

# get unique SEQ_ASSAY_ID from clinical file and unique genes covered
panels = aggregate(Hugo_Symbol ~ SEQ_ASSAY_ID,genieBedData[genieBedData$SEQ_ASSAY_ID %in% unique(genieClinData$SEQ_ASSAY_ID),],function(x) {length(unique(x))})
panels = panels[order(panels$Hugo_Symbol),]
# preallocate matrix
panelGenes = matrix(NA,nrow(panels),nrow(panels),dimnames=list(panels$SEQ_ASSAY_ID,panels$SEQ_ASSAY_ID))

# find intersect of unique genes 
for (i in 1:nrow(panels)) {
  for (j in 1:nrow(panels)) {
    panelGenes[i,j] = length(intersect(unique(as.character(genieBedData$Hugo_Symbol[genieBedData$SEQ_ASSAY_ID==panels$SEQ_ASSAY_ID[i]])),unique(as.character(genieBedData$Hugo_Symbol[genieBedData$SEQ_ASSAY_ID==panels$SEQ_ASSAY_ID[j]]))))
  }
}

# write to file
write.table(panelGenes,"~/Desktop/out.csv")
