#Analysis of E-MTAB-4304 Osteoarthritis cartilage dataset with PhenomeExpress

#takes around 15 mins to run entire script

#Requires expression data with fold changes and pvalues - annotated with official gene symbols for mouse or human.

require("igraph")
require("data.table")


setwd("data/PhenomeExpress")

#source the methods
source("./src/HeterogeneousNetwork.R")
source("./src/RHWN.R")
source("./src/runGIGA.R")
source("./src/runPhenoExpress.R")


load("./Networks/HumanExpressionData.RData")


#calculate the Pi value for scoring the nodes
expressionData$absFC=abs(expressionData$log2FoldChange)
expressionData$logAdjPval=-log10(expressionData$padj)
expressionData$Pi=(expressionData$absFC*expressionData$logAdjPval)

#load the HumanConsensusDB PPI network
load("./Networks/HumanNetwork.RData")

if (packageVersion("igraph") != "0.7.1") {
  HumanNetwork<-upgrade_graph(HumanNetwork)
  
}



#filter the netwrok based on expressed genes
presentList=na.omit(match(expressionData$name,V(HumanNetwork)$name))
OA.network=induced.subgraph(HumanNetwork,presentList)
OA.network=decompose.graph(OA.network)[[1]]

#filter expression table based on genes in the network
presentList=na.omit(match(V(OA.network)$name,expressionData$name))
expressionData=expressionData[presentList,]

#select the phenotypes from the UberPheno ontology - the Phenomiser tool and manual searching of the ontolgy by relevent keywords is helpful for this
Phenotypes=c("HP:0005086","HP:0001387","MP:0003724","MP:0003436")

#run Phenome Express - set inital subnetwork number to 15 rather than default 20 to give reasonably sized consensus sub-networks
OAResults=runPhenomeExpress(OA.network,expressionData,Phenotypes,"Human",max_number=10,sampleSize=1000)

#retrieve the significant sub-networks
subnetworks=OAResults[[1]]

#retrieve the table of p-values
sigTable=OAResults[[2]]
