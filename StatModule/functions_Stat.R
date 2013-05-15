#=============================================================================#
# ArrayAnalysis - affyAnalysisStat                                            #
# a tool for statistical analysis of Affymetrix expression data               #
#                                                                             #
# Copyright 2010-2011 BiGCaT Bioinformatics                                   #
#                                                                             #
# Licensed under the Apache License, Version 2.0 (the "License");             #
# you may not use this file except in compliance with the License.            #
# You may obtain a copy of the License at                                     #
#                                                                             #
# http://www.apache.org/licenses/LICENSE-2.0                                  #
#                                                                             #
# Unless required by applicable law or agreed to in writing, software         #
# distributed under the License is distributed on an "AS IS" BASIS,           #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #
# See the License for the specific language governing permissions and         #
# limitations under the License.                                              #
#=============================================================================#

#######################
## computeStatistics ##
#######################

computeStatistics <-function(normDataTable, descriptionFile, defaultContr=TRUE,
	matfileName=NULL, keepAnnotation=FALSE) {

  # Load the data	
  if(is.null(normDataTable)) stop("normDataTable has to be provided")
  if(is.null(descriptionFile)) stop("descriptionFile has to be provided")
  if(is.null(dim(normDataTable))) { # this is not a data.frame or matrix
			extension<-strsplit(normDataTable,"\\.")
			extension<-paste(".",extension[[1]][length(extension[[1]])],sep="")
		try(ndt<-readExtFile(normDataTable, extension))
		if(!exists("ndt")) stop("normDataTable format was not recognized")
			normDataTable <- ndt
      rm(ndt)
	}

  descdata<-readDescFile(descriptionFile)
  arrayNames<-descdata[,1]

  #in R, proper names should fulfill criteria, which already have been automatically applied to the column (array) names when reading the data
  arrayNames <- make.names(arrayNames)
  #if(sum(arrayNames != descdata[,1])>0) warning("One or more array names have been adapted to be valid R names")
  #line above: warning not needed as arraynames are in non of the outcome files; if this changes, uncomment the line
  
  experimentFactor<-descdata[,2]
  
  #in R, proper names cannot start with a number, so add a prefix if they do
  experimentFactor <- make.names(experimentFactor)
  #if(sum(experimentFactor != descdata[,2])>0) warning("One or more experimental group names have been adapted to be valid R names")
  #for now the code still does not work with unvalid R names as group names
  if(sum(experimentFactor != descdata[,2])>0) {
    stop("One or more experimental group names are invalid R names\nplease don't use special characters other than . and _ and don't start names with a number")
  }
  
  experimentFactor<-as.factor(experimentFactor)
  
# Annotation  
  
  # if the first column header of the normalized data table is not in the 
  # list of arrayNames, this will be defined as first column of annotation
  # in the result tables (usually it contains the gene/probeset IDs)
  firstColumn = NULL
  if(sum(arrayNames==colnames(normDataTable)[1])==0) {
	firstColumn <- as.matrix(normDataTable[,1])
	colnames(firstColumn) <- colnames(normDataTable)[1]
  }

  # other annotation columns are kept only if keepAnnotation is TRUE, in this 
  # case, all other columns of normDataTable that were not recognized as
  # array names are set to annotation. This annotation will be added at the  
  # end of the result tables
  annotation = NULL
  if(keepAnnotation) {
	headers <- colnames(normDataTable)
	notArrayNames <- apply(as.matrix(headers[2:length(headers)]),1,function(x) {
		sum(arrayNames==x)})
	annotation <- normDataTable[,(which(notArrayNames==0)+1)]
  }

# Normalized data
  # recreate the normDataTable based on the desc file to be sure that all  
  # array groups are defined and the columns are correctly ordered:
  normDataTable <- apply(as.matrix(arrayNames),1,function(x) {
	  as.numeric(normDataTable[,which(colnames(normDataTable)==x)])})
  normDataTable <- as.data.frame(normDataTable)
	colnames(normDataTable) <- arrayNames
  
  if(is.null(dim(normDataTable))) {
	stop("could not match array names from description file to normalized data
		file")
  }

# Compute statistical analysis	
  design <- model.matrix(~ 0 + experimentFactor)
  rownames(design) <- arrayNames
  colnames(design) <- gsub("experimentFactor","",colnames(design))

  #fit model
  fit <- lmFit(normDataTable,design=design) 
  #estimate eBayes values and perform moderated t-test
  fit <- eBayes(fit)
# Use contrast matrix to generate group comparisons  
  filesNew<-NULL
  if(defaultContr) { # always FALSE when coming from the web form.  
	if(length(levels(experimentFactor))<=4){
		cont.matrix <- defaultMatrix(design)
		defFiles<-saveComparison(cont.matrix,fit,normDataTable,
				annotation,firstColumn)
		filesNew<-c(filesNew,defFiles)
	}
	else
		print ("[[----Cant compute default matrices for more than four 
		groups----]]")
  }
  if(!is.null(matfileName)) {
  #be careful, since the saved contrast matrix still uses the not corrected names, also pass the not corrected names
	cont.matrix <- enterMatrix(matfileName,descdata[,2])
	advFiles<-saveComparison(cont.matrix,fit,normDataTable,
			annotation,firstColumn)
	filesNew<-c(filesNew,advFiles)
  }
  
  return(filesNew)
}

#=======================================================================#
# computeStatistics sub-functions: contrast matrix and save comparisons #
#=======================================================================#

#===== Read norm file, depending on the file type. Used from R and web ### Use gdata library
readExtFile <- function(filename,extension,outprint=FALSE){
  library(gdata)
  dataTable = NULL;
  switch(extension,
         ".txt" = dataTable<-trim(read.delim(filename, fill = FALSE, as.is=TRUE)),
         ".csv" = dataTable<-trim(read.csv(filename, fill = FALSE, as.is=TRUE)),
         ".xls" = dataTable<-trim(read.xls(filename, as.is=TRUE)),
         ".xlsx" = dataTable<-trim(read.xls(filename, as.is=TRUE))
	)
  if(is.null(dataTable)) stop(paste("extension",extension,"not regognised"))
  if(outprint){ # outprint is used from the website
		cols<-colnames(dataTable)
		for(i in 1:length(cols)){
			print(cols[i])
		}
   }else{
	return (dataTable)
   }
}

#===== Read desc file. Return array names and groups. Used from R and web
readDescFile <- function(descfile, outprint=FALSE){ 
  extension<-strsplit(descfile,"\\.")
  extension<-paste(".",extension[[1]][length(extension[[1]])],sep="")
  descript<-readExtFile(descfile, extension)

  if(is.null(dim(descript))) {
	stop("description file format was not recognized")
  }else{ # desc file format: either identical than from QC or just 2 columns
	if(dim(descript)[2]==2) descript<-cbind(rep("#",dim(descript)[1]),descript)
  }
  
  # when coming from the web form, "@" means "annotation"
  arrayNames <- as.vector(descript[descript[,3]!="@",2])
  arrayNames[is.numeric(arrayNames)]<-paste("X",arrayNames[is.numeric(arrayNames)],sep="")
  experimentFactor <- as.vector(descript[descript[,3]!="@",3])
  
  if(extension != ".txt"){
	arrayNames <- gsub(" $","",arrayNames)
	experimentFactor <- gsub(" $","",experimentFactor)
  }
  
  if(outprint){ # outprint is used from the website
    print(paste(paste(arrayNames,"\n",sep=""),collapse=""))
	print(paste(paste(experimentFactor,"\n",sep=""),collapse=""))
  }else{
	return<-cbind(arrayNames,experimentFactor)
  }
}


#===== Default comparison (for less than 5 experimental factors)
defaultMatrix <- function(design) {
  if(is.null(design)) stop("design has to be provided")
  levelsValue <- colnames(design)
  contrast.matrix <- factor(c(1:length(levelsValue)))
  for(i in 1:(length(levelsValue)-1)){
	for(j in (i+1):length(levelsValue)){
		contrastsValue <- paste(levelsValue[i],"-",levelsValue[j],sep="")
		cont.matrix <- makeContrasts(contrasts=contrastsValue,levels = design)
		contrast.matrix <- cbind(contrast.matrix,cont.matrix)
	}
  }
  return (contrast.matrix[,-1])
}

#===== Advanced comparison: enter custom contrasts 
enterMatrix <- function(matfileName,groupNames) {
  if(is.null(matfileName)) stop("matrixFile has to be provided")
  if(is.null(groupNames)) stop("groupNames has to be provided")
  try(mat<-read.delim(matfileName, header=FALSE))
  if(!exists("mat")) mat<-as.matrix(matfileName)
  if(!exists("mat")) stop("matrixFile is neither a string containing the data
	nor an existing file name")
  contrast.matrix <- NULL
  for(i in 1:(dim(mat)[1])){
    cont.matrix <- makeContrasts(contrasts = mat[i,1],levels = levels(as.factor(groupNames)))
	  contrast.matrix <- cbind(contrast.matrix,cont.matrix)
  }
  return (contrast.matrix)
}

#===== Create a toptable and save it into a file
toptable <- function(contrast.fit,i,normDataTable,fileName,
		firstColumn,annotation) {
  fileName <- gsub("\\(|\\)| ","",fileName)
  fileName <- gsub("\\*","x",fileName)
  fileName <- gsub("/","by",fileName)
  fileName <- gsub("\\(|\\)| ","",fileName)
  cat(paste("--[[ Saving table for contrast ",fileName, " ]]--\n", sep=""))
  
  # create the table
  toptab <- topTable(contrast.fit, adjust.method="BH", coef=i, 
		number=dim(normDataTable)[1], resort.by="P")
  # add annotation		
  if(!is.null(firstColumn)) {
	toptab <- cbind(firstColumn[match(rownames(toptab),
			rownames(normDataTable))],toptab)
	colnames(toptab)[1] <- colnames(firstColumn)
  }			
  if(!is.null(annotation)) {
	toptab <- cbind(toptab,annotation[match(rownames(toptab),
			rownames(normDataTable)),])
  }
  # modify the fold-changes so they are symetric around zero
  fc <- as.matrix(2^toptab[, "logFC"])
  fc[(toptab[, "logFC"] < 0)&(!is.na(toptab[, "logFC"]))] <- 
		-1/fc[(toptab[, "logFC"] < 0)&(!is.na(toptab[, "logFC"]))]
  colnames(fc) <- "Fold Change"
  m.col <- grep("logFC",colnames(toptab))
  toptab <- cbind(toptab[,1:m.col,drop=FALSE],fc,
		toptab[,(m.col+1):dim(toptab)[2],drop=FALSE])
  # write the table into a file		
  fileName <- paste(fileName,".txt",sep="")
  write.table(toptab,file=fileName,sep="\t",row.names=FALSE,append=FALSE, quote=FALSE)
  return (fileName)
}

#===== For each contrast, save the comparison into files
saveComparison <- function(cont.matrix,fit,normDataTable,
		annotation = NULL, firstColumn=NULL){
  if(is.null(cont.matrix)) stop("cont.matrix has to be provided")
  if(is.null(fit)) stop("fit object has to be provided")
  if(is.null(normDataTable)) stop("normDataTable has to be provided")

  contrast.fit <- contrasts.fit(fit, cont.matrix)
  contrast.fit <- eBayes(contrast.fit)
  
  files <- NULL
  if(is.null(dim(cont.matrix))){
	# Only one contrast
	fileName <- paste(paste(names(cont.matrix[cont.matrix!=0]),
		"*(",cont.matrix[cont.matrix!=0],")",sep=""), collapse = '+')
	i<-1
	files<-c(files,toptable(contrast.fit,i,normDataTable,fileName,
			firstColumn,annotation))
  }else{
	# More than one contrast
	  for (i in 1:length(colnames(cont.matrix))) {
		fileName <- paste("comp_",colnames(cont.matrix)[i])
		files<-c(files,toptable(contrast.fit,i,normDataTable,fileName,
				firstColumn,annotation))
	  }
  }
  return (files)
}


#####################
## createCutOffTab ##
#####################

createCutOffTab <- function(files,cutOffPval=NULL,cutOfflogFC=NULL,
		cutOffAveExpr=NULL) {
  if(is.null(files)) stop("vector of fit-model file names has to be provided")
  if(is.null(cutOffPval)&is.null(cutOfflogFC)&is.null(cutOffAveExpr)) {
	cat("--[[no cut-off was provided: skip the gene table creation]]--\n")
  }else{
	for(i in 1:length(files)){
		tab <- read.delim (files[i],header=TRUE)
		if(!is.null(cutOffPval)){
			tab<-tab[(tab[,"P.Value"]<=cutOffPval),]
			colnames(tab)[colnames(tab)=="P.Value"] <- paste("P.Value<=",
				cutOffPval)
		}
		if(!is.null(cutOfflogFC)){
			tab<-tab[(tab[,"logFC"]>=cutOfflogFC),]
			colnames(tab)[colnames(tab)=="logFC"]<-paste("logFC>=",cutOfflogFC)
		}
		if(!is.null(cutOffAveExpr)){
			tab<-tab[(tab[,"AveExpr"]>=cutOffAveExpr),]
			colnames(tab)[colnames(tab)=="AveExpr"] <- paste("AveExpr>=",
				cutOffAveExpr)
		}
		fileName <- paste(gsub(".txt","",files[i]),"_chosen_cutOff.txt",sep="")
		cat (paste("--[[ Saving table",fileName,"]]--\n"))
		write.table(tab,file = fileName ,sep="\t",row.names=FALSE, quote=FALSE)	
	}
  }
}
	
#################################
#####    plot histograms    #####
#################################

createPvalHist <- function(files) {
  #extract data from the comparison tables
  for(i in 1:length(files)){
	tab <- read.delim (files[i],header=TRUE)
	tabPval <- tab$P.Value
	cat (paste("--[[ Saving",files[i],"P.Value_hist.png ]]--\n"))
	png(paste(gsub(".txt","",files[i]),"_pvalue_hist.png",sep=""),width=1000,
		height=1000)
	hist (tabPval,main=paste("p-value histogram for",gsub(".txt","",files[i])),
		xlab = "p-values",col="blue")
	dev.off()
  }
}		

createFCHist <- function(files) {
  #extract data from the comparison tables
  for(i in 1:length(files)){
	tab <- read.delim (files[i],header=TRUE)
	tabFC <- tab$Fold.Change
	cat (paste("--[[ Saving",files[i],"FC_hist.png ]]--\n"))
	png(paste(gsub(".txt","",files[i]),"_FC_hist.png",sep=""),width=1000,
		height=1000)
	hist (tabFC,main=paste("adapted fold change histogram for",gsub(".txt","",files[i])),
		xlab = "adapted fold changes",col="green", breaks=60)
	dev.off()
  }
}


#####################################
###Create significant pvalue table###
#####################################

createPvalTab <- function(files,pvaluelist=c(0.001,0.01,0.05, 0.1),
		adjpvaluelist=0.05,foldchangelist=c(1.1,1.5,3), html=FALSE) {

library(R2HTML)
  if(is.null(files)) stop("Files list has to be entered")
  tab <- read.delim (files[1],header=TRUE)
  lengthTab<-length(tab$P.Value)

  write.table(t(cbind(c("No_of_genes:",lengthTab),rep("",2))),
		file="Summary_tables.txt",sep="\t",	row.names=FALSE,col.names=FALSE, quote=FALSE)
  if(html) HTML(t(cbind(c("No_of_genes:",length(tab$P.Value)),rep("",2))),file="Summary_tables.p.html",row.names = FALSE,append=FALSE)		

  ### P-values
  #expected
  expPval<-as.vector(apply(as.data.frame(pvaluelist),1,function(x) { 
			c(y<-round(lengthTab*x),0.5 * y,0.5 * y)}))

  P.Values<-c("Comparisons",paste("pVal<",rep(pvaluelist,each=3),
	c("tot","up","down")))
  Expected<-c("Expected NoOfGenes",expPval)
  P.Values<-cbind(P.Values,Expected)
   if(html) P.Values.html<-P.Values
  #extract data from the comparison tables
 
  for(i in 1:length(files)){
	tab <- read.delim (files[i],header=TRUE)
	row1 <- NULL
	for (count in 1:length(pvaluelist)){
		noOfPval <- up <- down <- 0
		#compare p-value with entered p-values 
		noOfPval <- sum(tab[,"P.Value"]<pvaluelist[count])
		up <- sum((tab[,"P.Value"]<pvaluelist[count])&(tab[,"logFC"]>=0))
		down <- sum((tab[,"P.Value"]<pvaluelist[count])&(tab[,"logFC"]<0))		
		row1 <- c(row1,noOfPval,up,down)
	}
	P.Values<-cbind(P.Values,c(gsub(".txt","",files[i]),row1))	
	if(html) P.Values.html<-cbind(P.Values.html,c(substr(gsub(".txt","",files[i]),1,30),row1))	
  }

  write.table(t(P.Values),file = "Summary_tables.txt",sep="\t",
			row.names=FALSE,col.names=FALSE,append=TRUE, quote=FALSE)
  if(html) HTML(t(P.Values.html),file="Summary_tables.p.html", row.names = FALSE)	
  
  ### Adjusted p-Values
  rowrest <- paste("Adj_pVal<",rep(adjpvaluelist,each=3),c("tot","up","down"))
  Adj.pValues <- c("Comparisons",rowrest)
  lengthTab<-length(Adj.pValues)
  if(html) Adj.pValues.html<-Adj.pValues
  
  for(i in 1:length(files)){
	tab <- read.delim (files[i],header=TRUE)
	row1 <- NULL
	for (count in 1:length(adjpvaluelist)){
		adjNoOfPval <- adjUp <- adjDown <- 0
		#compare adjusted p-values with entered adjusted p-values 
		adjNoOfPval <- sum(tab[,"adj.P.Val"]<adjpvaluelist[count])	
		adjUp <- sum((tab[,"adj.P.Val"]<adjpvaluelist[count])&(tab[,"logFC"]>=0))
		adjDown <- sum((tab[,"adj.P.Val"]<adjpvaluelist[count])&(tab[,"logFC"]<0))
		row1 <- c(row1,adjNoOfPval,adjUp,adjDown)
	}
	Adj.pValues<-cbind(Adj.pValues,c(gsub(".txt","",files[i]),row1))		
	if(html) Adj.pValues.html<-cbind(Adj.pValues.html,c(substr(gsub(".txt","",files[i]),1,30),row1))	
  }		
  write.table(t(cbind(rep("",lengthTab),Adj.pValues)),file = "Summary_tables.txt",sep="\t",
			row.names=FALSE,col.names=FALSE,append=TRUE, quote=FALSE)
  if(html) HTML(t(Adj.pValues.html),file="Summary_tables.p.html", row.names = FALSE)	
  
  ### Fold-Changes
  Fold.Changes <- c("Comparison",paste("|FC|>=",
	rep(foldchangelist,each=3),c("tot","up","down")))
  lengthTab<-length(Fold.Changes)
  if(html) Fold.Changes.html<-Fold.Changes
		
  for(i in 1:length(files)){	
	tab <- read.delim (files[i],header=TRUE)	
	row1 <- NULL
	for (count in 1:length(foldchangelist)){
		FCTot <- FCUp <- FCDown <- 0
		#compare FC with entered FC 
		FCTot <- sum(abs(tab[,"Fold.Change"])>=foldchangelist[count] & (!is.na(tab[,"Fold.Change"])))	
		FCUp <- length(tab[(tab$Fold.Change>=foldchangelist[count]) & (!is.na(tab$Fold.Change)),1]) #up
		FCDown <- length(tab[(tab$Fold.Change<=(-foldchangelist[count]))& (!is.na(tab$Fold.Change)),1])#down
		row1 <- c(row1,FCTot,FCUp,FCDown)
	}

	Fold.Changes<-cbind(Fold.Changes,c(gsub(".txt","",files[i]),row1))
	  if(html) Fold.Changes.html<-cbind(Fold.Changes.html,c(substr(gsub(".txt","",files[i]),1,30),row1))	
  }
  
	write.table(t(cbind(rep("",lengthTab),Fold.Changes)),
		file = "Summary_tables.txt",sep="\t", row.names=FALSE,
		col.names=FALSE,append=TRUE, quote=FALSE)
    if(html) HTML(t(Fold.Changes.html),file="Summary_tables.p.html", row.names = FALSE)	
  
  cat (paste("--[[ Saving Summary_tables.txt ]]--\n"))
}

