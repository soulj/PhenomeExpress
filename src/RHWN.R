#Function to perform the random walk based scoring on the heterogenous network

RWHN = function(GenePhenoTransition,heterogeneousNetwork) {
  
  #retrieve and assign the previously created matrices
  z2=heterogeneousNetwork[["matrix"]]
  x=heterogeneousNetwork[["geneMatrix"]]
  y=heterogeneousNetwork[["phenoMatrix"]]
  z=heterogeneousNetwork[["genePheno"]]
  m2=z2
  
  #calculate the transition probabilities for the hetrogenous gene-phenotype matrix using the heterogenous transition probability parameter
  #Phenotype to Gene
  present=z2[Matrix::rowSums(z2)>0,]
  present=(GenePhenoTransition*present)/Matrix::rowSums(present)
  z2=rBind(present,z2[Matrix::rowSums(z2)==0,])
  z2=z2[rownames(m2),]
  
  #do the same for the transpose
  #Gene to Phenotype
  z2t=t(z2)
  present=z2t[Matrix::rowSums(z2t)>0,]
  present=(GenePhenoTransition*present)/Matrix::rowSums(present)
  z2t=rBind(present,z2t[Matrix::rowSums(z2t)==0,])
  z2t=z2t[colnames(m2),]
  
  #caculate for the PPI matrix
  x_old=x
  #find the proteins which have a gene to phenotype interaction
  present.names=na.omit(match(colnames(z),rownames(x)))
  present.names=present.names[!duplicated(present.names)]
  present=x[present.names,]
  present=((1-GenePhenoTransition)*present)/Matrix::rowSums(present)
  #find the proteins which don't have a phenotype to gene interaction
  notPresent=x[-present.names,]
  notPresent=(notPresent)/Matrix::rowSums(notPresent)
  x=rBind(present,notPresent)
  x=x[rownames(x_old),]
  
  y_old=y
  #find the phenotypes which have a gene to phenotype interaction
  present.names=na.omit(match(rownames(z),rownames(y)))
  present.names=present.names[!duplicated(present.names)]
  present=y[present.names,]
  present=((1-GenePhenoTransition)*present)/Matrix::rowSums(present)
  #find the phenotypes which don't have a gene to phenotype interaction
  notPresent=y[-present.names,]
  notPresent=(notPresent)/Matrix::rowSums(notPresent)
  y=rBind(present,notPresent)
  y=y[rownames(y_old),]
  
  
  #create the full transition matrix
  #|geneMatrix  (x)   Gene2Pheno (z2t) |
  #|Pheno2Gene (z2)   PhenoMatrix (y)  |
    
  cols=c(colnames(x),colnames(y))
  m1 <- Matrix(data=0,nrow=(nrow(x)+nrow(y)), ncol=(nrow(x)+nrow(y)), sparse=TRUE,dimname=list(cols,cols))
  
  temp1=cBind(x,z2t)
  temp2=cBind(y,z2)
  m1=rBind(temp1,temp2)
  
  
  return(list(matrix=m1,geneMatrix=x,phenoMatrix=y))
  
}




calculate_RWHN=function(parameters,matrices,expressionTable,Phenotypes) {
  
  m1=matrices[["matrix"]]
  x=matrices[["geneMatrix"]]
  y=matrices[["phenoMatrix"]]
  
  ExpressionScores=as.data.frame(rownames(x))
  
  PersonalVectorWeighting=parameters[[1]]
  alpha=parameters[[2]]
  
  
  #set up the personalised vector for restart probabilities
  
  ExpressionScores$Score=expressionTable$Pi
  ExpressionScores$Norm=ExpressionScores$Score/sum(ExpressionScores$Score)
  ExpressionScores$Norm=ExpressionScores$Norm*PersonalVectorWeighting
  
  
  PhenotypeScores=as.data.frame(rownames(y))
  PhenotypeScores$Score=ifelse(PhenotypeScores[,1] %in% Phenotypes,(1-PersonalVectorWeighting)/length(Phenotypes),0)
  
  probabilityVector=c(ExpressionScores$Norm,PhenotypeScores$Score)
  probabilityVector=as.vector(t(probabilityVector))
  
  matrix=RWHN_compute(t(m1),probabilityVector,alpha)
  matrix=matrix[matrix$name %in% rownames(x),]
  matrix$geneRank=1:nrow(matrix)
  
  return(matrix)
  
}


getRWHN=function(matrices,parameters,expressionTable,Phenotypes) {
  #calculate of the matrix with each parameter combination
  result=apply(parameters,1,calculate_RWHN,matrices,expressionTable,Phenotypes)
  
  return(result)
  
}

#calcuate via iterative power method
RWHN_compute <- function(transitionMatrix,probabilityVector,alpha) {
  eps = 1/10^6
  iter = 0
  pi0=probabilityVector
  pi1=c(rep(0,length(probabilityVector)))

  while (sum(abs(pi0 - pi1)) > eps) {
    pi0 = pi1
    pi1 = ((1- alpha) * pi1 %*% transitionMatrix ) + ( alpha * probabilityVector)
    iter = iter + 1
  } 
  pi=as.matrix(t(pi1))
  geneRank=pi[order(-pi[,1]),]
  geneRank=as.data.frame(geneRank)
  geneRank$name=rownames(geneRank)
  
  return(geneRank)
}
