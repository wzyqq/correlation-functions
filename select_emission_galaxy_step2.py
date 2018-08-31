# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 09:30:58 2016
@author: wzy
"""
import numpy as np
import matplotlib.pyplot as plt
z=[]
ra=[]
dec=[]
i1=[]
i2=[]
data=np.loadtxt("select whole properties")
for j in data:
	if( (23.2<=j[4]<=24.2) and (0.05<(j[4]-j[5])<0.35) ):
		if(0<=j[6]<=24):
			z.append(j[2])
			ra.append(j[0])
			dec.append(j[1])
		elif( (j[6]==-99.0) and (0<=j[7]<=24) ):
			z.append(j[2])
			ra.append(j[0])
			dec.append(j[1])
Z=np.array(z)
RA=np.array(ra)
DEC=np.array(dec)
result=np.array([RA,DEC,Z]).T
np.savetxt("select emission",result,"%10.7f")