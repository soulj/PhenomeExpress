# Method to run PhenomeExpress and to return the significant sub-netwroks

#Input: a igraph object (PPI network) with the node name matching the name column in the expression table,the expression data with a column named "Pi" containing the Expression Score (Pi-value), a vector of Phenotypes, the species ("Human" or "Mouse")
#optional - max size of subnetworks and FDR threshold of significant sub-networks


#takes approximately 5 mins to run

runPhenoExpress=function(inputGraph,expressionTable,Phenotypes,species,max_number=20,FDR=0.05,sampleSize=10000){
  
  require(data.table)
  require(Matrix)
  require(igraph)
  
    
  #overlay the expression values on the graph
  V(inputGraph)$Pi=expressionTable$Pi
  
  #set up heterogeneous network
  heterogeneousNetwork=getHeterogeneousNetwork(inputGraph,species)
  
  #set up the parameters for the Random walk
  alpha=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)
  transitionBetween=c(0.5,0.6,0.7,0.8)
  relativeWeighting=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)
  parameters=expand.grid(relativeWeighting,alpha)
  
  #set up the transition matrices
  transitionMatrix=lapply(transitionBetween,RWHN,heterogeneousNetwork)
  
  #run RWHN using all the matrices to score the nodes
  Scores=sapply(transitionMatrix,getRWHN,parameters,expressionTable,Phenotypes)
  scoreNames=lapply(Scores,function(k) match(k$name,V(inputGraph)$name))
  Scores=Map(cbind, Scores, NodeID = scoreNames)
  
  #run GIGA to find high scoring sub-networks
  Clusters=sapply(Scores,runGIGA,inputGraph,max_number)
  Clusters=unlist(Clusters,recursive=F)
  Clusters=lapply(Clusters,function(x) CJ(x,x))
  Clusters=do.call("rbind", Clusters)
  B=as.data.frame(Clusters)
  
  #calulate the consensus sub-networks
  C=table(B[,2])
  B <- crossprod(table(B[1:2]))
  diag(B) <- 0
  B=sweep(B,MARGIN=2,C,`/`)
  B[B<1]=0
  B=graph.adjacency(B,mode="undirected",weighted=TRUE)
  Clusters=decompose.graph(B,min.vertices=5)
  
  Clusters=lapply(seq_along(Clusters),getNames,inputGraph,Clusters)
  
  #sample random sub-networks to determine signifcance
  pvalueres=c()
  for (i in 1:length(Clusters)){
    presentList=na.omit(match(V(Clusters[[i]])$name,V(inputGraph)$name))
    graph=induced.subgraph(inputGraph,presentList)
    randomGraphs=replicate(sampleSize,randomSample(inputGraph,vcount(graph)))
    count=length(randomGraphs[randomGraphs>=sum(V(graph)$Pi)])
    pvalue=(count+1)/sampleSize
    pvalueres=c(pvalueres,pvalue)
    
  }
  
  #filter the sub-networks by the user selected FDR
  pvaluetable=data.frame(Number=c(1:length(pvalueres)),pvalue=pvalueres)
  pvaluetablesig=pvaluetable[pvaluetable$pvalue<=FDR,]
  sigNetworks=pvaluetablesig$Number
  subnetworks=list()
  for(i in 1:length(sigNetworks)) {
    subnetworks[[i]]=Clusters[[sigNetworks[[i]]]]
  }
  
  return(list(subnetworks,pvaluetablesig))
  
}

#utility function to get match names of nodes
getNames=function(x,inputGraph,Clusters) {
  presentList=na.omit(match(V(Clusters[[x]])$name,V(inputGraph)$name))
  return(induced.subgraph(inputGraph,presentList))}


#utility function to randomly find sub-network of set size and return the total score
randomSample<-function(graph,ncount) {
  
  node <- sample.int(vcount(graph), 1)
  selected <- rep(NA,ncount)
  selected[[1]]<-node
  i<-2 
  
  while(i<=ncount) {
    neigh<-neighbors(graph,node)
    
    node <- sample(neigh,1)
    
    if (sum(node==selected,na.rm=TRUE)==0) {
      selected[[i]]<-node
      i<-i+1
      #print(paste0("We now have ",i," nodes."))
    }
  }
  return(sum(V(induced.subgraph(graph, selected))$Pi))
}