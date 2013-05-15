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
version_nb <- "1.0.0"
cat("Script run using R version ",R.Version()$major,".",R.Version()$minor,
  " and affyAnalysisStat version_",version_nb,"\n",sep="")

#set memory to maximum on windows 32bit machines
if(length(grep("w32",R.Version()$os,fixed=TRUE))>0) memory.size(4095)

###############################################################################
# Load R libraries and affyAnalysisStat functions  							  #
###############################################################################
require("limma", quietly = TRUE)
require("gdata", quietly = TRUE) #trim function
print ("Libraries has been loaded")

correctDIR <- function(d) { 
  lastChar <- substr(d,nchar(d),nchar(d))
  if((lastChar != "/") && (lastChar != "/")) d <- paste(d,"/",sep="")
  return(d)
}
if(exists("SCRIPT.DIR")) SCRIPT.DIR <- correctDIR(SCRIPT.DIR)
if(exists("WORK.DIR")) WORK.DIR <- correctDIR(WORK.DIR)
reload <- function() {
source(paste(SCRIPT.DIR,"functions_Stat.R",sep=""))
print ("Functions have been loaded")
}
reload();
if(!exists("libdir")) { # libdir exists only for GenePattern usage
  setwd(SCRIPT.DIR)
  setwd(WORK.DIR)
}

###############################################################################
# Statistical analysis with limma 
###############################################################################
print ("Statistical analysis with limma")
if(defaultContr==FALSE & is.null(matfileName)) stop("Need at least to compute default contrasts or to enter a contrast matrix file")
files <- computeStatistics(normDataTable,descriptionFile,defaultContr,
		matfileName,keepAnnotation) 

###############################################################################
# Creating table based on the cut-offs
###############################################################################
if(cutOffTable){
	print ("Creating table based on the cut-offs")
	createCutOffTab(files,cutOffPval,cutOfflogFC,cutOffAveExpr)
}
###############################################################################
# Plotting p-value histograms
###############################################################################
if(plotPvalHist){
	print ( "p-value histograms" )
	createPvalHist (files)
}

###############################################################################
# Create significant genes table
###############################################################################
if(summaryTable){
	print ( "Creating table of significant genes" )
	createPvalTab(files,pvaluelist,adjpvaluelist,foldchangelist, html)
}
