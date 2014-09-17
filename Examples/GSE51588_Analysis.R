#Analysis of GSE51588 OA bone dataset with PhenomeExpress

#takes around 15 mins to run entire script

require("Matrix")
require("igraph")
require("data.table")
require("hgug4112a.db") # for the microarray analysis
require("org.Hs.eg.db")
require("limma")
require("annotate")
require("BioNet") # for comparison purposes - not needed by PhenomeExpress
require("VennDiagram") # for making the Venn diagram figures
require("RCytoscape") # also requires cytoscape v2.8 to be open with the Cytoscape RPC plugin active


setwd("~/PhenomeExpress")

#source the methods
source("./src/HeterogeneousNetwork.R")
source("./src/RHWN.R")
source("./src/runGIGA.R")
source("./src/runPhenoExpress.R")

setwd("~/PhenomeExpress/GSE51588/")

#analyse the data
targets <- read.delim("./targets.txt")
images <- read.maimages(targets,source="agilent",green.only=TRUE)
images.processed <- backgroundCorrect(images, method="normexp", offset=16)
images.processed <- normalizeBetweenArrays(images.processed, method="quantile")
images.processed <- avereps(images.processed, ID=images.processed$genes$ProbeName)
dataMatrix=as.data.frame(images.processed$E)
dataMatrix2=ifelse(dataMatrix>=5,1,0)
dataMatrix2=dataMatrix[rowSums(dataMatrix2)>4,]
detected=dataMatrix[!is.na(getEG(rownames(dataMatrix),"hgug4112a.db")),]
dataMatrix3=dataMatrix[rownames(dataMatrix) %in% rownames(detected),]
Treatment=targets$Treatment
design=model.matrix(~factor(Treatment))
colnames(design)=c("Normal","OA")

#fit the linear model
fit <- lmFit(dataMatrix3, design)
#calculate pvalues
fit <- eBayes(fit)
#calculate BH correction p values and store results table
results=topTable(fit,coef="OA",number=Inf)

#Annotate the probes with Entrez gene IDs
genes=as.data.frame(getEG(rownames(results), "hgug4112a.db" ))
colnames(genes)=c("EntrezID")
results=merge(results,genes,by="row.names")


#collapse duplicated genes
require(data.table)
dt <- data.table(results)
dt=dt[, .SD[which.max(abs(logFC)),], by=EntrezID]
dt=as.data.frame(dt)

setwd("~/PhenomeExpress/")

#Use the David and Uniprot ID maps to match the EntrezID to Swissprot for the PPI network
Young_EnteztoSwiss_via_Uniprot <- read.delim("./GSE51588/GenenamesEntreztoUniprot_via_UniProt.txt")
Young_EnteztoSwiss_via_David <- read.delim("./GSE51588/GenenamesEntreztoUniprot_via_David.txt", dec=",")
Young_EnteztoSwiss_via_David=Young_EnteztoSwiss_via_David[,1:2]
Young_EnteztoSwiss=rbind(Young_EnteztoSwiss_via_David,Young_EnteztoSwiss_via_Uniprot)
Young_EnteztoSwiss=Young_EnteztoSwiss[!duplicated(Young_EnteztoSwiss),]

dt=merge(results,Young_EnteztoSwiss,by.x="EntrezID",by.y="From")
colnames(dt)[9]="name"

#calculate the Pi value for scoring the nodes
dt$absFC=abs(dt$logFC)
dt$logAdjPval=-log10(dt$adj.P.Val)
dt$Pi=(dt$absFC*dt$logAdjPval)

#load the HumanConsensusDB PPI network
load("./Networks/ConsensusDB_graph.RData")

#filter the netwrok based on expressed genes
presentList=na.omit(match(dt$name,V(ConsensusDB_graph)$name))
OA.network=induced.subgraph(ConsensusDB_graph,presentList)
OA.network=decompose.graph(OA.network)[[1]]

#filter expression table based on genes in the network
presentList=na.omit(match(V(OA.network)$name,dt$name))
dt=dt[presentList,]

#select the phenotypes from the UberPheno ontology - the Phenomiser tool and manual searching of the ontolgy by relevent keywords is helpful for this
Phenotypes=c("HP:0005086","MP:0003724","HP:0002829","HP:0100777","MP:0004983","ZP:0006539","MP:0002896","MP:0005006")

#run Phenome Express - set inital subnetwork number to 15 rather than default 20 to give reasonably sized consensus sub-networks
OAResults=runPhenomeExpress(OA.network,dt,Phenotypes,"Human",max_number=15)

#retrieve the significant sub-networks
subnetworks=OAResults[[1]]

#retrieve the table of p-values
sigTable=OAResults[[2]]

#collapse all the nodes in the subnetworks from PhenomeExpress
nodes=c()
for(i in 1:length(subnetworks)) {
  tempGraph=subnetworks[[i]]
  nodes=c(nodes,V(tempGraph)$name)
}

#load the results from JActiveModules and GIGA - run externally, subnetworks >= 5 nodes kept
boneJAM <-read.table("./JActiveModules/boneJMSA2107", quote="\"")
boneJAM=boneJAM[!duplicated(boneJAM$V1),]
GIGA <- read.delim("./GIGA/boneGIGA.txt", header=F)

#run BioNet for comparison
pval=dt$P.Value
names(pval)=dt$name
b <- fitBumModel(pval, plot = FALSE)
scores <- scoreNodes(network = OA.network, fb = b,fdr = 1e-4) #FDR produces similar sized module to max sized PhenomeExpress sub-network
module <- runFastHeinz(OA.network, scores)

#count the number of seed Phenotype annotated proteins present in all the sub-networks for each tool

#First get the gene to phenotype associations for labelling seed nodes
z=getHeterogeneousNetwork(OA.network,"Human")[["genePheno"]] # note contains all proteins - including ones not present in network
phenoAnnotated=z[rownames(z) %in% Phenotypes,]
phenoAnnotated=phenoAnnotated[,colSums(phenoAnnotated)>0]
phenoAnnotated=colnames(phenoAnnotated)

#calculate the number of seed phenotype annotated genes for each tool
no.Seeds.PhenomeExpress=table(ifelse(nodes %in% phenoAnnotated,1,0))
no.Seeds.boneJAM=table(ifelse(boneJAM %in% phenoAnnotated,1,0))
no.Seeds.GIGA=table(ifelse(GIGA$V2 %in% phenoAnnotated,1,0))
no.Seeds.BioNet=table(ifelse(V(module)$name %in% phenoAnnotated,1,0))

#make a Venn diagram of protein in subnetworks from each tool
nodeList=list(PhenomeExpress=nodes,JActivemodules=boneJAM,GIGA=GIGA$V2,BioNet=V(module)$name)
venn.diag=venn.diagram(nodeList,fill = c("red", "green","blue","purple"),alpha = c(0.5, 0.5,0.5,0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,  filename=NULL )
grid.draw(venn.diag)


#send all the sub-networks from PhenomeExpress to cytoscape
#colours the nodes according to the fold change
#black border if directly annotated to seed phenotype 

#useful to assign the node with the entrez ID as well - for downstream analysis in cytoscape i.e mapping to genenames or functional annotation
V(OA.network)$EntrezID=as.character(dt$EntrezID)


for(i in 1:length(subnetworks)) {
  presentList=na.omit(match(V(subnetworks[[i]])$name,V(OA.network)$name))
  tempGraph=induced.subgraph(OA.network,presentList)
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
