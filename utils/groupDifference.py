#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:51:28 2018

@author: yusuf
"""


import argparse
import numpy as np
import scipy.stats as stt
import matplotlib.pylab as plt
from helpers import calculateZScore,fdr,calculateGroupDifference,drawHistogram2Dataset,drawBoxPlot,drawViolinPlot


# -r 2_b_strStrDistance_Desikan86_det/results/direct_none_subtract.res -o ./ -s /home/yusuf/data/TBI/resources/qa_commonList_Deterministic_Desikan86.txt -mt similarity -st dist  -pe png -pt box -al group
# -r 1_a_strFuncCoupling_Desikan86_det/results/accuracy_direct_logScaleEdgesStructure_normalizeEdges_edgesIgnoreDiag_LinAss_0.res -o ./ -s /home/yusuf/data/TBI/resources/qa_commonList_Deterministic_Desikan86.txt -st dist  -pe png -pt violin
#get parameters
parser = argparse.ArgumentParser(description='calculate correlation between the age and matching accuracies of the subjects')
parser.add_argument('-rM','--similarityRelToMale', help='file path to the results of the matching experiment relative to males', required=True)
parser.add_argument('-sF','--similarityRelToFemale', help='file path to the results of the matching experiment relative to females', required=True)
parser.add_argument('-s','--subjectsList', help='path to the file that contains path to the connectomes', required=True)
parser.add_argument('-o','--outputFolder', help='file path to save the distribution of values for the two populations', required=True)
parser.add_argument('-t','--title', help='title for plot', required=False, type=str, default='') #You can escape white space in the title from command line with "\ " as in <two\ words>
parser.add_argument('-pe','--plotExtension', help='extension of the plot file', required=False,type=str,choices=['png','svg'],default='png')
parser.add_argument('-pt','--plotType', help='draw violin or box plot', required=False,type=str,choices=['violin','box'],default='violin')

args = vars(parser.parse_args())
similarityRelToMale=args['similarityRelToMale']
similarityRelToFemale=args['similarityRelToFemale']
subjectListPath=args['subjectsList']
outputFolder=args['outputFolder']
title=args['title']
plotExtension=args['plotExtension']
plotType=args['plotType']

##############load subjet IDs for which we have a score#########################
with open(subjectListPath,"r") as f:
    cohort =  f.read().splitlines()
numScans=len(cohort)

idMale=[]
idFemale=[]
orderMale=[]
orderFemale=[]
ageMale=[]
ageFemale=[]
fullSubjectList=[]

for order,line in enumerate(cohort):
    fullSubjectList.append(line.split('\t')[0])
    if (line.split('\t')[2]=='m'):
        idMale.append(line.split('\t')[0])
        ageMale.append(line.split('\t')[1])
        orderMale.append(order)
    elif (line.split('\t')[2]=='f'):
        idFemale.append(line.split('\t')[0])
        ageFemale.append(line.split('\t')[1])
        orderFemale.append(order)
        
##############load results of the structure-function coupling experiment###################
fileContentRelToMale =  open(similarityRelToMale,"r").read().splitlines()
fileContentRelToFemale =  open(similarityRelToFemale,"r").read().splitlines()
numNodes = int(fileContentRelToMale[1].split('\t')[0])
numSubjects = int(fileContentRelToMale[1].split('\t')[1])
measureType = str(fileContentRelToMale[3]) 
scoreName = str(fileContentRelToMale[5]) 

scores=np.zeros((numSubjects,2))

for i in range(numSubjects):
    scores[i][0] = float(fileContentRelToMale[7].split('\t')[i])
    scores[i][1] = float(fileContentRelToFemale[7].split('\t')[i])




########################calculate group difference########################
colors=['aqua','darkorchid'] #dodgerblue , magenta

orderGroups=[orderMale,orderFemale]
nameGroups=["male","female"]


# comparisons: mm vs ff, mm vs mf, ff vs fm, mf vs fm
scores_comparison=[ [scores[orderMale][:,0],scores[orderFemale][:,1]], 
                             [scores[orderMale][:,0],scores[orderMale][:,1]] , 
                             [scores[orderFemale][:,1],scores[orderFemale][:,0]] , 
                             [scores[orderMale][:,1],scores[orderFemale][:,0]] ]
names_comparisons=["mmVSff: males similarity relative to males vs females similarity relative to females", 
                   "mmVSmf: males similarity relative to males vs males similarity relative to females", 
                   "ffVSfm: females similarity relative to females vs females similarity relative to males", 
                   "mfVSfm: males similarity relative to females vs females similarity relative to males"]
shortNames_comparisons=["mmVSff","mmVSmf","ffVSfm","mfVSfm"]
isPaired=[False,True,True,False]

numComparisons=len(scores_comparison) 

effectSize=np.zeros(numComparisons)
pValue=np.zeros(numComparisons)
corrected_pValue=np.zeros(numComparisons)

#calculate group difference btw healthy vs patients at three time points
for i in range(numComparisons): 
    effectSize[i], pValue[i] = calculateGroupDifference(scores_comparison[i][0],scores_comparison[i][1],parametric=True,paired=isPaired[i])
    histogramPath=outputFolder+"histograms/"+measureType+"_"+shortNames_comparisons[i]+"."+plotExtension
    drawHistogram2Dataset(scores_comparison[i][0],scores_comparison[i][1],effectSize[i],pValue[i],"ES","pVal",shortNames_comparisons[i].split("VS")[0],shortNames_comparisons[i].split("VS")[1],"",histogramPath,xLabel=measureType,yLabel="frequency",color1=colors[0],color2=colors[1])

corrected_pValue=fdr(pValue)

reportFile=open(outputFolder+scoreName+"_groupDifference.txt",'w')
reportFile.write("============== "+measureType+"  "+scoreName+" score ==============\n")
for i in range(numComparisons): 
    if(corrected_pValue[i]<=0.05):
        reportFile.write("**")
    reportFile.write("\t"+names_comparisons[i]+"\t:\t%0.2f (%f -- %f) \n" % (effectSize[i],pValue[i],corrected_pValue[i]))
reportFile.close()


##########draw boxplot or violin plot matching scores of systems###################
colors=['#2096BA','#351C4D', '#AB3E16','#F5AB99','#849974','#F7DFD4'] #shutter blue, nightfall, rust, tropical pink, fresh,  macaron, 

outputPath=outputFolder+scoreName+"_groupDifference."+plotExtension
data=[scores[orderMale][:,0],scores[orderMale][:,1],scores[orderFemale][:,0],scores[orderFemale][:,1]]
dataLabels=['m-m','m-f','f-m','f-f']


minY=min([min(l) for l in data])
maxY=max([max(l) for l in data])
offset=(maxY-minY)/5.0
if plotType=="box":
    drawBoxPlot(data,dataLabels,title,outputPath,xLabel='',yLabel=measureType,colors=colors,rotation=0,plotScatter=True,yLim=[minY-offset/2.0,maxY+offset])
elif plotType=="violin":
    drawViolinPlot(data,dataLabels,title,outputPath,xLabel='',yLabel=measureType,colors=colors,rotation=0,plotScatter=True)

