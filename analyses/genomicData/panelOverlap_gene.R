# libraries
library(synapseClient)

# SAGE login
synapseLogin(username=,password=) # you could do this with prompt or coded, but will need to enter something here

# get bed and clinical data from SAGE
genieClinData = read.delim(synGet("syn7871792")@filePath,header=TRUE)
genieBedData = read.delim(synGet("syn7871816")@filePath,header=TRUE)

# get unique SEQ_ASSAY_ID from clinical file
panels = as.character(unique(genieClinData$SEQ_ASSAY_ID))
# preallocate matrix
panelGenes = matrix(NA,length(panels),length(panels),dimnames=list(panels,panels))

# find intersect of unique genes 
for (i in 1:length(panels)) {
  for (j in 1:length(panels)) {
    panelGenes[i,j] = length(intersect(unique(as.character(genieBedData$Hugo_Symbol[genieBedData$SEQ_ASSAY_ID==panels[i]])),unique(as.character(genieBedData$Hugo_Symbol[genieBedData$SEQ_ASSAY_ID==panels[j]]))))
  }
}
