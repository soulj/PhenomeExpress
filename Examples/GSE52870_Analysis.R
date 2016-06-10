#Analysis of GSE52870 PAX5 dataset with PhenomeExpress

#takes around 15 mins to run

#Requires expression data with fold changes and pvalues - annotated with official gene symbols for mouse or human.




require("Matrix")
require("igraph")
require("data.table")

setwd("data/PhenomeExpress/")

#source the methods
source("./src/HeterogeneousNetwork.R")
source("./src/RHWN.R")
source("./src/runGIGA.R")
source("./src/runPhenoExpress.R")


#load the expressionData
load("./Networks/MouseExpressionData.RData")

#load the high confidence mouse PPI network from STRING
load("./Networks/MouseNetwork.RData")

if (packageVersion("igraph") != "0.7.1") {
  MouseNetwork<-upgrade_graph(MouseNetwork)
  
}

expressionData=na.omit(expressionData)
expressionData$name<-as.character(expressionData$name)
presentList=na.omit(match(expressionData$name,V(MouseNetwork)$name))


#Use pre-existing networks filter based on genes found in the transcriptomics experiment
pax5.network=induced.subgraph(MouseNetwork,presentList)
pax5.network=decompose.graph(pax5.network)[[1]]
presentList=na.omit(match(V(pax5.network)$name,expressionData$name))

#filter the expression data based on proteins present in the network
expressionData=expressionData[presentList,]


#calculate the Pi value for use in the node scoring stage
expressionData$Pi=abs(expressionData$log2FoldChange)*-log10(expressionData$padj)
expressionData$absFC=abs(expressionData$log2FoldChange)

#select the phenotypes from the UberPheno ontology - the Phenomiser tool and manual searching of the ontolgy by relevent keywords is helpful for this
Phenotypes=c("HP:0004812","MP:0012431","HP:0012191","MP:0008211","MP:0008189")

#run Phenome Express
LeukResults=runPhenomeExpress(pax5.network,expressionData,Phenotypes,"Mouse",max_number=10,sampleSize=1000)
