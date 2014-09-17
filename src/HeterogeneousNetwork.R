# Function to load the pre-made PhenomeNetwork and Bipartate gene to phenotype network.
#INPUT::igraph object containing PPI network 
# species- "Mouse" or "Human"

#Returns a hetrogenous network matrix, a PPI and Phenome ajacency matrix and a gene to phenotype bipartiate network


getHeterogeneousNetwork=function(inputGraph,species="Human") {
  require(Matrix)
  require(igraph)
  
  
  
  if (species == "Mouse") {
    
    load("~/GenePhenoGraph_Mouse.RData")
    load("~/UberPhenoGraph_Mouse.RData")
    
    
  }
  
  else {
    
    load("~/GenePhenoGraph.RData")
    load("~/UberPhenoGraph.RData")
    
  }
  
  
  #Make an adjency matrix from the graphs weighted by confidence or semantic similarity
  
  x=get.adjacency(inputGraph,attr="Confidence")
  y=get.adjacency(UberPhenoGraph,attr="V3")
  z=get.adjacency(GenePhenoGraph)
  
  #Set the phenotype-gene matrix to contain all the genes and phenotypes in the phenotype and gene networks.
  genephenomap=bipartite.mapping(GenePhenoGraph)
  z=get.incidence(GenePhenoGraph, types=genephenomap$type, attr=NULL, names=TRUE, sparse=T)
  present=colnames(z) %in% rownames(x)

  #Remove any proteins that are associated with a phenotype but are not in the PPI network as they can't be used (the associtions are pre-filtered by the PhenomeNetwork)
  z_filt=z[,present]         
  
  #find the proteins that and phenotypes that don't have any heterogenous associations, but are present in the PPI network or the Phenome network to add to the matrix
  y_filt=(!(rownames(y) %in% rownames(z)))
  y_filt2=rownames(y)[y_filt]
  
  x_filt=(!(rownames(x) %in% colnames(z)))
  x_filt2=rownames(x)[x_filt]
  
  #create the gene to phenotype association based matrix
  cols=c(colnames(z_filt),x_filt2)
  rows=c(rownames(z_filt),y_filt2)
  m1 <- Matrix(data=0,nrow=(nrow(z_filt)+length(y_filt2)), ncol=(ncol(z_filt)+length(x_filt2)), sparse=TRUE,dimname=list(rows,cols))
  #fill in the matrix with the association data
  m1[1:nrow(z_filt),1:ncol(z_filt)]=z_filt[,]
  
  #sort the matrix to match the protein and phenotype matrices
  m2=m1[rownames(y),rownames(x)]
  
  
  return(list(matrix=m2,geneMatrix=x,phenoMatrix=y,genePheno=z))
  
}