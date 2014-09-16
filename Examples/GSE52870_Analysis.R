#Analysis of GSE52870 PAX5 dataset with PhenomeExpress

#takes around 5-10mins to run
require("Matrix")
require("igraph")
require("data.table")
require("DESeq2") # for the RNA-seq analysis
require("BioNet") # for comparison purposes - not needed by PhenomeExpress
require("VennDiagram") # for making the Venn diagram figures
require("RCytoscape") # also requires cytoscape v2.8 to be open with the Cytoscape RPC plugin active

setwd("~/PhenomeExpress")

#source the methods
source("./src/HeterogeneousNetwork.R")
source("./src/RHWN.R")
source("./src/runGIGA.R")
source("./src/runPhenoExpress.R")

#calculate the FPKM using the effective gene length and the counts per gene
GSE52870_Pax5Restoration.GenewiseCounts <- read.delim("./GSE52870/GSE52870_Pax5Restoration-GenewiseCounts.txt")
countmatrix=GSE52870_Pax5Restoration.GenewiseCounts[,3:8]
rownames(countmatrix)=GSE52870_Pax5Restoration.GenewiseCounts$EntrezID
genelength=GSE52870_Pax5Restoration.GenewiseCounts$GeneLength
FPKMtable=(countmatrix * 10^9) /(colSums(countmatrix) * genelength)
FPKMtable=ifelse(FPKMtable>1,1,0)
countmatrix=countmatrix[Matrix::rowSums(FPKMtable)>2,]


#use DESeq2 to analyse the raw data
colData=data.frame(colnames=colnames(countmatrix),condition=c(rep("PAX5KD",3),rep("PAX5Rescue",3)))
dds=DESeqDataSetFromMatrix(countData=countmatrix,colData=colData,design=~condition)
dds$condition=factor(dds$condition, levels =c ( "PAX5KD","PAX5Rescue" ))
dds2=DESeq(dds)


#get the expression table with the fold changes and p values
res=results(dds2)
dt=as.data.frame(res[order (res$log2FoldChange),])
dt$EntrezID=rownames(dt)


#Anotate the genes with SwissProt names to match the network node names
Young_EnteztoSwiss_via_Uniprot <- read.delim("./GSE52870/GenenamesEntreztoUniprot_via_UniProt.txt")
Young_EnteztoSwiss_via_David <- read.delim("./GSE52870/GenenamesEntreztoUniprot_via_David.txt", dec=",")
Young_EnteztoSwiss_via_David=Young_EnteztoSwiss_via_David[,1:2]
Young_EnteztoSwiss=rbind(Young_EnteztoSwiss_via_David,Young_EnteztoSwiss_via_Uniprot)
Young_EnteztoSwiss=Young_EnteztoSwiss[!duplicated(Young_EnteztoSwiss),]

#note 1 entrez gene maps to more than one protein
dt=merge(dt,Young_EnteztoSwiss,by.x="EntrezID",by.y="From")
dt=na.omit(dt)
colnames(dt)[8]="name"

#load the high confidence mouse PPI network from STRING
load("./Networks/HCString_Mouse_Graph.RData")
presentList=na.omit(match(dt$name,V(HCString_Mouse)$name))

#Use pre-existing networks filter based on genes found in the transcriptomics experiment
pax5.network=induced.subgraph(HCString_Mouse,presentList)
pax5.network=decompose.graph(pax5.network)[[1]]
presentList=na.omit(match(V(pax5.network)$name,dt$name))

#filter the expression data based on proteins present in the network
dt=dt[presentList,]
dt=na.omit(dt)

#calculate the Pi value for use in the node scoring stage
dt$Pi=abs(dt$log2FoldChange)*-log10(dt$padj)
dt$absFC=abs(dt$log2FoldChange)

#select the phenotypes from the UberPheno ontology - the Phenomiser tool and manual searching of the ontolgy by relevent keywords is helpful for this
Phenotypes=c("HP:0004812","MP:0012431","HP:0012191","MP:0008211","MP:0008189")

#run Phenome Express
LeukResults=runPhenoExpress(pax5.network,dt,Phenotypes,"Mouse")

#retrieve the significant sub-networks
subnetworks=LeukResults[[1]]

#retrieve the table of p-values
sigTable=LeukResults[[2]]

#collapse all the nodes in the subnetworks from PhenomeExpress
nodes=c()
for(i in 1:length(subnetworks)) {
  tempGraph=subnetworks[[i]]
  nodes=c(nodes,V(tempGraph)$name)
}

#load the results from JActiveModules and GIGA - run externally, subnetworks >= 5 nodes kept
leukJAM <-read.table("./JActiveModules/leukJM2107", quote="\"")
leukJAM=leukJAM[!duplicated(leukJAM$V1),]
GIGA <- read.delim("./GIGA/leukGIGA.txt", header=F)

#run BioNet for comparison
pval=dt$pvalue
names(pval)=dt$name
b <- fitBumModel(pval, plot = FALSE)
scores <- scoreNodes(network = pax5.network, fb = b,fdr = 1e-25) #FDR produces similar sized module to max sized PhenomeExpress sub-network
module <- runFastHeinz(pax5.network, scores)

#count the number of seed Phenotype annotated proteins present in all the sub-networks for each tool

#First get the gene to phenotype associations for labelling seed nodes
z=getHeterogeneousNetwork(pax5.network,"Mouse")[["genePheno"]] # note contains all proteins - including ones not present in network
phenoAnnotated=z[rownames(z) %in% Phenotypes,]
phenoAnnotated=phenoAnnotated[,colSums(phenoAnnotated)>0]
phenoAnnotated=colnames(phenoAnnotated)

#calculate the number of seed phenotype annotated genes for each tool
no.Seeds.PhenomeExpress=table(ifelse(nodes %in% phenoAnnotated,1,0))
no.Seeds.leukJAM=table(ifelse(leukJAM %in% phenoAnnotated,1,0))
no.Seeds.GIGA=table(ifelse(GIGA$V2 %in% phenoAnnotated,1,0))
no.Seeds.BioNet=table(ifelse(V(module)$name %in% phenoAnnotated,1,0))

#make a Venn diagram of protein in subnetworks from each tool
nodeList=list(PhenomeExpress=nodes,JActivemodules=leukJAM,GIGA=GIGA$V2,BioNet=V(module)$name)
venn.diag=venn.diagram(nodeList,fill = c("red", "green","blue","purple"),alpha = c(0.5, 0.5,0.5,0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,  filename=NULL )
grid.draw(venn.diag)


#send all the sub-networks from PhenomeExpress to cytoscape
#colours the nodes according to the fold change
#black border if directly annotated to seed phenotype 

#useful to assign the node with the entrez ID as well - for downstream analysis in cytoscape i.e mapping to genenames or functional annotation
V(pax5.network)$EntrezID=as.character(dt$EntrezID)


for(i in 1:length(subnetworks)) {
  presentList=na.omit(match(V(subnetworks[[i]])$name,V(pax5.network)$name))
  tempGraph=induced.subgraph(pax5.network,presentList)
  FC=dt[na.omit(match(V(tempGraph)$name,dt$name)),]
  V(tempGraph)$logFC=FC$log2FoldChange
  
  seedAnnotatedGenes=ifelse(V(tempGraph)$name %in% phenoAnnotated,1,0)
  V(tempGraph)$Seed=seedAnnotatedGenes
  
  
  #do the network creation stuff
  
  #convert the igraph object to a graphNEL object and intialise the attributes
  tempGraph.NEL=igraph.to.graphNEL(tempGraph)
  tempGraph.NEL=initEdgeAttribute(tempGraph.NEL,"Confidence","numeric",0)
  tempGraph.NEL=initEdgeAttribute(tempGraph.NEL,"weight","numeric",0)
  tempGraph.NEL=initNodeAttribute(tempGraph.NEL,"logFC","numeric",0)
  tempGraph.NEL=initNodeAttribute(tempGraph.NEL,"Seed","numeric",0)
  tempGraph.NEL=initNodeAttribute(tempGraph.NEL,"EntrezID","char",0)
  
  nodeDataDefaults(tempGraph.NEL, "label") <- "name"
  nodeData(tempGraph.NEL,V(tempGraph)$name,"label") = V(tempGraph)$name
  tempGraph.NEL=initNodeAttribute(tempGraph.NEL,"label","char","name")
  
  
  #Open the cytoscape window and send the graph
  cw1 <- new.CytoscapeWindow (paste("PhenoExpress",as.character(i),sep=""), graph=tempGraph.NEL)
  #display the graph
  displayGraph (cw1)
  #select the layout
  layoutNetwork (cw1, layout.name='force-directed')
  
  #colour according to the logFC
  control.points <- c(-5,0,5)
  node.colors <- c ("#00AA00", "#00FF00", "#FFFFFF", "#FF0000", "#AA0000")
  setNodeColorRule (cw1, node.attribute.name='logFC', control.points, node.colors, mode='interpolate')
  
  
  setDefaultBackgroundColor (cw1, '#FFFFFF')
  
  #set the nodeborder to correspond to the seed phenotype annotated genes
  data.values <- c ("1", "0")
  line.widths = c ("15","1")
  setNodeBorderWidthRule (cw1, 'Seed', data.values, line.widths)
}
