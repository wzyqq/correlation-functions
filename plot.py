import numpy as np
import matplotlib.pyplot as plt
#whole		=	np.loadtxt('../whole_sample/whole_sample_42')
#random		=	np.loadtxt('../whole_random_sample/whole_random_sample_42')
whole		=	np.loadtxt('select emission')
random		=	np.loadtxt('emission random sample')
ra_whole	=	whole[:,0]
dec_whole	=	whole[:,1]
z_whole		=	whole[:,2]	
ra_random	=	random[:,0]
dec_random	=	random[:,1]
z_random	=	random[:,2]

plt.subplot(4,2,1)
plt.plot(ra_whole,dec_whole,'r.',markersize=1)
plt.subplot(4,2,2)
plt.plot(ra_random,dec_random,'b.',markersize=0.05)
plt.subplot(4,2,3)
plt.hist(z_whole,bins=np.linspace(1.5,2,501))
plt.subplot(4,2,4)
plt.hist(z_random,bins=np.linspace(1.5,2,501))
plt.subplot(4,2,5)
plt.hist(ra_whole,bins=np.linspace(30,39,100))
plt.subplot(4,2,6)
plt.hist(ra_random,bins=np.linspace(30,39,100))
plt.subplot(4,2,7)
plt.hist(dec_whole,bins=np.linspace(-7,-4,100))
plt.subplot(4,2,8)
plt.hist(dec_random,bins=np.linspace(-7,-4,100))
plt.show()