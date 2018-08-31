# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 09:30:58 2016
@author: wzy
"""
ngalaxy	=	17471
rgalaxy	=	500000
r_all	=	1242313
import numpy as np
import matplotlib.pyplot as plt
y=np.loadtxt('select emission')
position=np.loadtxt('random_mask')
c=np.random.randint(0,ngalaxy,rgalaxy)
c_list1=[]
position1=[]
position2=[]
d=np.random.randint(0,r_all,rgalaxy)
for k in range(rgalaxy):
    c_list1.append(y[c[k]][2])
    position1.append(position[d[k]][0])
    position2.append(position[d[k]][1])
z=np.array(c_list1)
POSITION1=np.array(position1)
POSITION2=np.array(position2)
result=np.array([position1,position2,z]).T
np.savetxt("emission random sample", result,fmt='%10.7f')
#plt.hist(z,bins=np.linspace(1.5,2,6))
#plt.title('random z distrbution')
#plt.show()