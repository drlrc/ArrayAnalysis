#=============================================================================#
# ArrayAnalysis - ArrayAnalysisPath                                           #
# a tool for visualisation of expression set on a pathways                    #
#                                                                             #
# Copyright 2010-2013 BiGCaT Bioinformatics                                   #
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
# Set directories and server address
###############################################################################

#SCRIPT.DIR <- paste(getwd(),"",sep="/")
#PATHWAY.DIR <- "/home/x_array/arrayanalysis.dev/arrayanalysis/pathways/"
#PATHWAY.DIR <- paste(getwd(),"test/pathways/",sep="/")
PATHWAY.DIR <- "/home/anwesha/PathVisio-Data/mouse-pathways"
#WORK.DIR <- paste(getwd(),"",sep="/")

###############################################################################
# Set input parameters                          				             
###############################################################################

inputFile <- "/home/anwesha/PathVisio-Data/arrayanalysistest/data.txt"
#inputFile <- "/home/x_array/arrayanalysis.dev/arrayanalysis/temp/example_set1_2011-08-22_13-54_44_Stat/control_24h-treated_24h.txt"
#inputFile <- paste(getwd(),"myTest.txt",sep="/")
#dbFile <- paste(getwd(),"Hs_Derby_2010601.bridge",sep="/")
dbDir <- "/home/anwesha/PathVisio-Data/gene-databases"
#species<-"HomoSapiens"
resultDir <- "/home/anwesha/PathVisio-Data/arrayanalysistest"


gSample <- "log2FC"
gColor <- "blue,white,red" 
gValue <- "-2,0,2"
rSample <- "pvalue"
rColor <- "green"
rExpr <- "[pvalue]<0.05"
zExpr <- "[pvalue]<0.05"

source(paste(SCRIPT.DIR,"run_affyAnalysisPath.R",sep=""),local=TRUE)

# DESCRIPTION : 
# inputFile <- name of the file containing statistical comparison table 
#				[give full path if the file is not present in your WORK.DIR]
# dbFile 	<- database file name (.bridgedb file corresponding to your species) 
#				[give full path if the file is not present in your WORK.DIR]
# gSample <- string representing a list of unique sample names (file headers), 
#				on which a gradient is applied. Sample names are separated by ";" 
#				Attribute NULL value if no gradient is required
# gColor 	<- string representing a list of list of colors: 2 to 3 colors are 
#				applied to each gradient. The first separator ";" splits the 
#				groups of colors belonging to the same gradient and the second 
#				"," splits the colors for each given gradient.
# gValue 	<- string representing a list of list of values corresponding to 
#				gColor and is organized the same way.
# rSample <- string representing a list of unique ids (headers no present in 
#				gColorSet or any string). These color sets may contain one or 
#				several color rules. The color set separator is ";". 
#				Attribute NULL value if no color rule is required
# rColor 	<- rColor is a string representing a list of list of colors: for 
#				each rule, only one color is set but many rules may be defined 
#				for one color set. The first separator ";" splits the  groups 
#				of colors belonging to the same color set the second "," splits 
#				the colors for each given rule.
# rExpr		<- string representing a list of list of boolean expression 
#				corresponding to rColor and is organized the same way. 
# zExpr 	<- boolean expression for zscore calculation
