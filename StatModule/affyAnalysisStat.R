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

###############################################################################
# Set input parameters                          				             
###############################################################################

#see descriptions below
normDataTable <- "" 
descriptionFile <- ""
keepAnnotation <- FALSE
defaultContr <- TRUE
matfileName <- NULL
cutOffTable <- TRUE
cutOffPval <- 0.5
cutOfflogFC <- 2
cutOffAveExpr <- 7
plotPvalHist <- TRUE
summaryTable <- TRUE
pvaluelist <- c(0.001,0.01,0.05, 0.1)
adjpvaluelist <- c(0.05)
foldchangelist <- c(1.1,1.5,3)

###############################################################################
# Set directories 		                          				              
###############################################################################

# Change these Paths if needed #
SCRIPT.DIR <- getwd()
WORK.DIR <- getwd()

###############################################################################
# Run functions 		                          				              
###############################################################################

source(paste(SCRIPT.DIR,"run_affyAnalysisStat.R",sep=""))

###############################################################################
# PARAMETER DESCRIPTION                         				              #
###############################################################################

# "automated calls" means usage from GenePattern and arrayanalysis.org
# Note that flag letters are only used in these automated calls

# n = normDataTable = file or data.frame or matrix object containing normalized  
#		data - headers = TRUE (mandatory)
# d = descript description file describing the array names and experimental 
#        groups - headers = TRUE (mandatory)	
# k = keepAnnotation = boolean for adding annotation in result table. Before 
#		using it, make sure your the arrays present in the normalized data
#		table are all fully described in the description file.
# c = defaultContr = boolean for deafult comparisons as in if 3 groups are 
#		present then grps : 3-2,2-1,3-1 is given by default. Option not 
#		available for more than four groups. Always FALSE with web usage.
# m = matfileName = file containing contrast matrix for advanced comparisons eg  
#		provided in example folder or string containing comparisons. headers = 
#		FALSE. For eg. : c("control+24hr+48hr)/3","control+(24hr+48hr)/2")
# C = cutOffTable = boolean for creating the cut-offs table
# p = cutOffPval = p-value cutoff 
# f = cutOfflogFC = logFC cutoff 
# e = cutOffAveExpr = Average expression cutoff [table contains genes with p-value 
#		<= provided value & logFC >= provided value & average expression >= 
#		provided value]
# H = plotPvalHist = boolean for plotting the p-values histograms
# S = summaryTable = boolean for computing the summary tables
# P = pvaluelist = provide list of p-values, separated by comma, for the result
#		table. You can enter a as many p-values as you want.
# A = adjpvaluelist = provide list of adjusted p-values [table contains the number 
#		of genes falling in each class from all the comparisons ]
# F = foldchangelist = provide list of fold-changes, separated by comma, for the 
#		result table. You can enter as many fold-changes as you want.
