#!/usr/bin/env python
# @Author: Jun Ding
# @Date: 05/20/2020

import pdb,sys,os
from File import *
import networkx as nx
import math 
from BioUtils import BioList
from scipy.stats import mannwhitneyu
import seaborn as sns
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import random 

# get all the nodes from the input list: X
def getNodes(X):
	XProteins=[item[:2] for item in X]
	Nodes=[]
	for i in XProteins:
		for j in i:
			if j not in Nodes:
				Nodes.append(j)
	return Nodes

# get all the edges from the input list : X
def getEdges(X):
	dEG={}
	for i in X:
		key='%s,%s'%(i[0],i[1])
		dEG[key]=float(i[-1])
	return dEG
	
# get the score of path X	
def getScore(G,X):
	ts=0
	for i in range(len(X)-1):
		si=G.edges[X[i],X[i+1]]['weight']
		ts+=si
	return ts

def extractEdges(Proteins,Edges):
	outEdges={}
	for i in Edges:
		[A,B]=i.split(",")
		if A in Proteins and B in Proteins:
			outEdges[i]=round(Edges[i],2)
	return outEdges
			

#-----------------------------------------------------------------------
# The main program starts here!
def main():
	parser=argparse.ArgumentParser()
	parser.add_argument("-i","--input_genes",required=True,help="the input gene list for the SARS-CoV-2 association analysis")
	parser.add_argument("-v","--viral_host_interactions",required=True,help="The interactions between viral proteins and host proteins")
	parser.add_argument("-p","--host_protein_interactions",required=True,help="The interaction between host proteins")
	parser.add_argument("-o","--output",required=True,help="The specified output directory")
	parser.add_argument("-b","--background",required=False,default=1000, help="Integer, the number of background proteins/genes used to calculate the connectivity significance of input genes")
	args=parser.parse_args()
	
	fnCOVID19=args.viral_host_interactions
	fnPPI=args.host_protein_interactions
	fnMarker=args.input_genes
	output=args.output
	try:
		bgN=int(args.background)
	except:
		print("-b must be integer, please check your input!")
		sys.exit(0)
		
	#----------------------------
	# read in covid19 interactions
	print("processing input files ...")
	COVInteractions=TabFile(fnCOVID19).read("\t")[1:]
	SProteins=list(set([item[1] for item in COVInteractions]))
	source=list(set([item.split(" ")[0] for item in SProteins]))
	#pdb.set_trace()
	
	dCI={}
	for i in COVInteractions:
		dCI[i[1].split(" ")[0]+","+i[3]]=[i[1],i[3]]
	COVInteractions=[[item[1].split(" ")[0],item[3],item[-2]] for item in COVInteractions]
	COVProteins=[]
	for i in COVInteractions:
		COVProteins+=i[:2]
	COVProteins=list(set(COVProteins))
	# read in PPIs
	PPI=TabFile(fnPPI).read("\t")[1:]
	PPI=[[item[0],item[2],item[3]] for item in PPI]
	# read in markers
	markers=LineFile(fnMarker).read()
	markers=[item for item in markers if item!='']
	
	print("building protein networks ...")
	AllN=getNodes(COVInteractions)+getNodes(PPI)+markers
	
	AllNodes=[]
	for i in AllN:
		if i not in AllNodes:
			AllNodes.append(i)
			
	#AllNodes=list(set(AllN))
	AEdges=getEdges(COVInteractions+PPI)
	AllEdges={}
	for i in AEdges:
		iKey=i
		iScore=AEdges[i]
		if iScore>0:
			iScore=math.log(1.0/iScore,2)
			AllEdges[iKey]=iScore
	G=nx.Graph()
	for i in AllNodes:
		G.add_node(i)

	for i in AllEdges:
		[A,B]=i.split(',')
		si=AllEdges[i]
		G.add_edge(A,B,weight=si)
		
	print("learning optimal path from SARS-CoV-2 to input genes ...")
	targets=markers
	SPaths=[]
	AllProteins=[]
	AllUProteins=[]
	for i in source:
		for j in targets:
			if nx.has_path(G,i,j):
				pij_all=nx.all_shortest_paths(G,i,j)
				sij_all=[]
				for pijk in pij_all:
					sijk=getScore(G,pijk)
					sij_all.append([sijk,pijk])
				sij_all=sorted(sij_all,key=lambda x:x[0])
				[sij,pij]=sij_all[0]
				
				AllProteins+=pij
				uij=dCI[pij[0]+','+pij[1]]+pij[2:]
				AllUProteins+=uij
				SPaths.append([i,j,sij,uij])
				print(j)
		
	AllProteins=list(set(AllProteins))
	outEdges=extractEdges(AllProteins,AEdges)
	network=[]
	for i in outEdges:
		[A,B]=i.split(",")
		iScore=outEdges[i]
		if i in dCI:
			[A,B]=dCI[i]
		ni=[A,iScore,B]
		network.append(ni)
	
	NodeInfo=[['Node','Type']]
	for i in AllUProteins:
		if i in SProteins:
			NodeInfo.append([i,'virus'])
		elif i in COVProteins:
			NodeInfo.append([i,'source'])
		elif i in targets:
			NodeInfo.append([i,'target'])
		else:
			NodeInfo.append([i, 'Intermediate'])

	print("finding the optimal paths from SARS-CoV-2 proteins to all host proteins")
	APaths=[]
	ct=0
	
	# fix the random initialization 
	
	random.seed(a=10)
	bgNodes=random.sample(AllNodes[1:],bgN)
	
	
	for i in source:
		for j in bgNodes:
			if nx.has_path(G,i,j):
				pij_all=nx.all_shortest_paths(G,i,j)
				sij_all=[]
				for pijk in pij_all:
					sijk=getScore(G,pijk)
					sij_all.append([sijk,pijk])
				sij_all=sorted(sij_all,key=lambda x:x[0])
				[sij,pij]=sij_all[0]
				APaths.append([i,j,sij,pij])
			ct+=1
			print(ct)
						
	SPath=[float(item[2]) for item in SPaths]
	APath=[float(item[2]) for item in APaths]
	pv1=mannwhitneyu(SPath,APath,alternative='less')
	print("p-value: %s"%(pv1[1]))
	XX=[SPath,APath]
	#df=pd.DataFrame(data=XX)
	#df.index=['HLH Genes','All Genes']
	
	if os.path.exists(output)==False:
		os.mkdir(output)
		
	print("writing results ...")
	print("exporting inferred paths ...")
	BioList(SPaths).ex2File("%s/SPaths.txt"%(output),"\t")
	BioList(APaths).ex2File("%s/APaths.txt"%(output),"\t")
	print("exporting network file (.sif) ...")
	BioList(network).ex2File("%s/network.sif"%(output),'\t')
	print("exporting network node attribute file (.txt)")
	BioList(NodeInfo).ex2File("%s/NodeInfo.txt"%(output),'\t')
	sns.boxplot(data=XX)
	plt.xticks(range(len(XX)),['Input Genes','All Genes'])
	plt.ylabel("Connectivity Score")
	plt.savefig("%s/Connectivity.pdf"%(output))
	
if __name__=="__main__":
	main()
	
