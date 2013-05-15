#=============================================================================#
# ArrayAnalysis - ArrayAnalysisPath                                           #
# a tool for visualisation of expression data on a pathways                   #
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

gexFile <- paste(inputFile,".pgex",sep="")

## Load library
library(methods)
library(XMLRPC)
print("Libraries are loaded")

## Set working directory
#setwd (WORK.DIR)

## Set address
server <- "http://localhost:7777"

##Make pgex
print("Importing a Data file")
pgex <- xml.rpc(server,"PathVisio.importData", inputFile, dbDir, resultDir)
print (pgex)	

## Create Visualization
print("Creating visualisation")
xml.rpc(server,"PathVisio.createVisualization", "/home/anwesha/PathVisio-Data/arrayanalysistest/data.txt.pgex", gSample, gColor, gValue, rSample, rColor, rExpr)
print (vis)

## Get results
print("Calculating Pathway Statistics")
result <- xml.rpc(host,"PathVisio.calculatePathwayStatistics", PATHWAY.DIR, gexFile, dbDir, zExpr, resultDir)
print (result)
print("End")


