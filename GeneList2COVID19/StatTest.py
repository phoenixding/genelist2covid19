"""
@Name: StatTest
@Author:Jun Ding
@Date: Mar.3,2011
@Version: 1.0
module 1: HyperGeometricTest
function: perform hypergeometric testing to get the p-value
usage:
python HyperGeometricTest(N,M,n,m)

module 2: pbinom
function:compute the cumulative probability densiuty function of the binomial distribution up P(X<=x)

"""
import math,pdb,sys,os
import numpy as np 
from scipy.stats import hypergeom as hyperG

def hyperGPV(N,M,n,m):
	pv=1-hyperG.cdf(m-1,N,n,M)
	return pv
	
