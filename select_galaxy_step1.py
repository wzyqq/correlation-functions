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
m=[]
g=[]
r=[]
i1=[]
i2=[]
data=np.loadtxt("/home/wzy/Desktop/data/cfht_w1_ks_mass.cut.cat.masked")
for j in data:
	if(1.5<=j[51]<=2 and (0<=j[28]<=24 or 0<=j[30]<=24)):
		z.append(j[51])
		ra.append(j[52])
		dec.append(j[53])
		m.append(j[11])
		g.append(j[26])
		r.append(j[27])
		i1.append(j[28])
		i2.append(j[30])
M=np.array(m)
Z=np.array(z)
RA=np.array(ra)
DEC=np.array(dec)
G=np.array(g)
R=np.array(r)
I1=np.array(i1)
I2=np.array(i2)
result=np.array([RA,DEC,Z,M,G,R,I1,I2]).T
np.savetxt("select whole properties",result,"%10.7f")
	
