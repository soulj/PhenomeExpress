#R implementation of the the GIGA perl script - much faster
#finds sub-networks from a ranked list

runGIGA=function(Scores,inputGraph,max_number=20,neighboursAdj,neighboursAdjSelf) {
  
   
  Scores=Scores[order(Scores[,3]),]
  #get the seeds
  seedCluster=initaliseCluster(Scores,inputGraph,neighboursAdj,neighboursAdjSelf)
  #expand the from the seeds
  expandedCluster=clusterExpand(Scores,seedCluster,completedCluster=c(),inputGraph,max_number,neighboursAdjSelf)
  
  #keep adding nodes until max size or pvalue doesn't improve
  while(length(expandedCluster[["cluster"]])>0) {
    seedCluster=addNewMin(Scores,expandedCluster,inputGraph,neighboursAdjSelf)
    expandedCluster=clusterExpand(Scores,seedCluster,seedCluster[["completedCluster"]],inputGraph,max_number,neighboursAdjSelf)
    
  }
  #return the final sub-networks
  finalCluster=getFinalClusters(Scores,expandedCluster,inputGraph)
  finalCluster=lapply(finalCluster,function(x) V(inputGraph)$name[x])
  return(finalCluster)
  
}


getNeighboursAdj<-function(i,adj){
  neighbours<-which(adj[i,]>0)
  return(neighbours)
  
}


#function to find the local min in the PPI network
initaliseCluster = function(Scores,inputGraph,neighboursAdj,neighboursAdjSelf) {
  
  
  #identify the neihbours of each node and find the local minima - those with no connections to a node with a higher rank.
  
  size=length(neighboursAdjSelf)
  NodeDataRankNamed=Scores[,3]
  names(NodeDataRankNamed)=Scores$geneRank
  neighbours=lapply(neighboursAdjSelf,function(x) NodeDataRankNamed[x])
  localMinIndex=lapply(neighbours, function(x) which.min(as.numeric(names(x))))
  localMin=which(localMinIndex==1)
  
  #expand the initial cluster
  
  cluster=as.list(localMin)
  names(cluster)=V(inputGraph)$name[localMin]
  oldcluster=cluster
  cluster=lapply(cluster,function(x) NodeDataRankNamed[x])
  #oldpvalue=lapply(cluster,getPvalue,size)
  oldpvalue=rep(1,length(cluster))
  names(oldpvalue)<-names(cluster)
  clustermax=localMin
  
  neighbours=neighboursAdj[localMin]
  neighbours=sapply(neighbours,function(x) NodeDataRankNamed[x])
  neighbours=sapply(neighbours,function(x) x[sort.list(as.numeric(names(x)))])
  neighbours=lapply(neighbours,'[',1)
  cluster<-mapply(c, cluster, neighbours, SIMPLIFY=FALSE)
  clustermax=sapply(cluster,function(x) names(x[2]))
  
  return(list(clustermax=clustermax,oldpvalue=oldpvalue,cluster=cluster,oldcluster=oldcluster))
  
}


#function to expand the existing sub-network
clusterExpand = function(Scores,seedCluster,completedCluster,inputGraph,max_number,neighboursAdjSelf) {
  NodeDataRankNamed=Scores[,3]
  names(NodeDataRankNamed)=Scores$geneRank
  clustermax=seedCluster[["clustermax"]]
  oldpvalue=seedCluster[["oldpvalue"]]
  cluster=seedCluster[["cluster"]]
  oldcluster=seedCluster[["oldcluster"]]
  
  cluster_new=cluster
  tobePValued=c()
  
  size=vcount(inputGraph)
  
  #keep expanding while you can (neighbours less than current max rank)
  while (length(cluster_new)>0) {
    
    neigbours=getNeigbours(cluster_new,neighboursAdjSelf)
    lastexpansion=cluster_new
    
    #assign the expanded cluster with the rank
    neigbours=lapply(neigbours,function(x) NodeDataRankNamed[x])
    #extract from 1 to the max rank in each cluster
    clustermax=as.numeric(clustermax)
    expansion=lapply(seq_along(neigbours),function(x)  neigbours[[x]][(as.numeric(names( neigbours[[x]]))<=clustermax[x])])
    names(expansion)=names(cluster_new)
    #remove any clusters that are too big
    expansion_tooBig=lapply(expansion,function(x) length(x)<max_number)
    cluster_new=expansion[unlist(expansion_tooBig)]
    completedCluster=c(completedCluster,oldcluster[!unlist(expansion_tooBig)])
    oldcluster=oldcluster[unlist( expansion_tooBig)]
    clustermax=clustermax[unlist(expansion_tooBig)]
    #seperate the clusters that have changed in size
    clusterlength=lapply(cluster_new,length)
    oldclusterlength=lapply(lastexpansion,length)
    continueCluster=lapply(names(cluster_new),function(x) clusterlength[[x]]>oldclusterlength[[x]])
    if (length(continueCluster>1)) {
      tobePValued=c(tobePValued,cluster_new[!unlist(continueCluster)])
    }
    cluster_new=cluster_new[unlist(continueCluster)]
    clustermax=clustermax[unlist(continueCluster)]
    oldcluster=oldcluster[unlist(continueCluster)]
  }
  
  if (length(tobePValued)>0) {
    oldcluster=seedCluster[["oldcluster"]]  
    oldcluster=oldcluster[match(names(tobePValued),names(oldcluster))]
    oldpvalue=oldpvalue[names(oldpvalue) %in% names(tobePValued)]
    currentPValue=lapply(tobePValued,getPvalue,size)
    expansion_improved=lapply(names(currentPValue),function(x) currentPValue[[x]]<oldpvalue[[x]])
    names(expansion_improved)<-names(currentPValue)
    cluster=tobePValued[unlist(expansion_improved)]
    worseClusters=oldcluster[!unlist(expansion_improved)]
    completedCluster=c(completedCluster,worseClusters)
    oldpvalue=currentPValue[unlist(expansion_improved)]
    #oldcluster=oldcluster[unlist(testing11)]
  }
  else {cluster=NULL}
  
  return(list(oldpvalue=oldpvalue,cluster=cluster,completedCluster=completedCluster))
}


addNewMin= function(Scores,expandedCluster,inputGraph,neighboursAdjSelf) {
  
  NodeDataRankNamed=Scores[,3]
  names(NodeDataRankNamed)=Scores$geneRank
  cluster=expandedCluster[["cluster"]]
  oldpvalue=expandedCluster[["oldpvalue"]]
  completedCluster=expandedCluster[["completedCluster"]]
  
  oldcluster=cluster
  newNeighbours=getNeigbours(cluster,neighboursAdjSelf)
  
  
  clustermax=c()
  newMin=lapply(newNeighbours,function(x) NodeDataRankNamed[x])
  names(newMin)=names(cluster)
  newMin=lapply(newMin,function(x) x[sort.list(as.numeric(names(x)))])
  addition=lapply(seq_along(newMin),function(x) newMin[[x]][!(newMin[[x]] %in% cluster[[x]])][1]) 
  cluster=mapply(c, cluster, addition, SIMPLIFY=FALSE)
  clustermax=sapply(cluster,function(x) max(as.numeric(names(x))))
  
  return(list(clustermax=clustermax,oldpvalue=oldpvalue,cluster=cluster,oldcluster=oldcluster,completedCluster=completedCluster))
  
}

getPvalue=function(x,size) {
  x=as.numeric(names(x))
  maxRank=max(x)
  pvalue=maxRank/size
  if (length(x)==1) return(pvalue)
  for (k in 1:(length(x)-1)) {
    pvalue=pvalue*((maxRank-k)/(size-k))
  }
  return(pvalue)
}

getNeigbours = function(cluster,neighboursAdj) {
  
  newNeighbours=lapply(cluster,function(x) unlist(neighboursAdj[ unlist(x)]))
  newNeighbours=lapply(newNeighbours,function(x) x[!duplicated(x)])
  
  return(newNeighbours)
}

getFinalClusters=function(Scores,expandedCluster,inputGraph){
  size=vcount(inputGraph)  
  NodeDataRankNamed=Scores[,3]
  names(NodeDataRankNamed)=Scores$geneRank
  completedCluster=expandedCluster[["completedCluster"]]
  completedCluster=lapply(completedCluster,function(x) NodeDataRankNamed[x])
  completedClusterpValue=lapply(completedCluster, getPvalue,size)
  completedCluster=completedCluster[completedClusterpValue<1/size]
  completedClusterpValue=lapply(completedCluster, getPvalue,size)
  completedClusterpValue=data.frame(name=names(completedCluster),pvalue=unlist(completedClusterpValue))
  completedClusterpValue=completedClusterpValue[order(completedClusterpValue$pvalue),]
  completedClusterpValue$localmin=match(completedClusterpValue$name,V(inputGraph)$name)
  
  localmin=completedClusterpValue[,3]
  finalCluster=c()
  localminused=c()
  alreadydone=c()
  
  
  
  for (i in 1:nrow(completedClusterpValue)) {
    if ((completedClusterpValue[i,3] %in% localminused)==FALSE) {
      candidateCluster=completedCluster[[as.character(completedClusterpValue[i,1])]]
      localminpresent=candidateCluster[candidateCluster %in% localmin ]
      localminused=c(localminpresent,localminused)
      finalCluster=c(finalCluster,list(candidateCluster))
    }
  }
  
  
  return(finalCluster)
  
}
