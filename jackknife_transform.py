import numpy as np
import astropy
from astropy.cosmology import FlatLambdaCDM
subsample = 100

for i in range(subsample):
	observe_file  = "emission_sample/emission_sample_%d"%(i)
	random_file = "emission_random_sample/emission_random_sample_%d"%(i)	
	observe_position  = "emission_sample_position/emission_sample_position_%d"%(i)
	random_position = "emission_random_sample_position/emission_random_sample_position_%d"%(i)	
	observe = np.loadtxt(observe_file)
	random  = np.loadtxt(random_file)
	#i_limit=24.0
	#emission_observe=1919
	observelen=len(observe)
	#emission_random=122605
	randomlen=len(random)
	observe_result=np.zeros([observelen, 3])
	random_result=np.zeros([randomlen, 3])
	cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
	observe_distance=np.zeros(observelen)
	random_distance=np.zeros(randomlen)
	for k in range(observelen):
		observe_distance[k]=cosmo.comoving_distance(observe[k,2]).value
		observe_result[k,0]=observe_distance[k]*np.cos(observe[k,1]*3.1415927/180.0)*np.cos(observe[k,0]*3.1415927/180.0)
		observe_result[k,1]=observe_distance[k]*np.cos(observe[k,1]*3.1415927/180.0)*np.sin(observe[k,0]*3.1415927/180.0)
		observe_result[k,2]=observe_distance[k]*np.sin(observe[k,1]*3.1415927/180.0)
		#observe_result[k,3]=(observe_distance[k]*10**((i_limit-observe[k,3])/5.0))**3
	for k in range(randomlen):
		random_distance[k]=cosmo.comoving_distance(random[k,2]).value
		random_result[k,0]=random_distance[k]*np.cos(random[k,1]*3.1415927/180.0)*np.cos(random[k,0]*3.1415927/180.0)
		random_result[k,1]=random_distance[k]*np.cos(random[k,1]*3.1415927/180.0)*np.sin(random[k,0]*3.1415927/180.0)
		random_result[k,2]=random_distance[k]*np.sin(random[k,1]*3.1415927/180.0)
		#random_result[k,3]=(random_distance[k]*10**((i_limit-random[k,3])/5.0))**3
	np.savetxt(observe_position,observe_result,fmt='%14.8lf\t%14.8lf\t%14.8lf')
	np.savetxt(random_position,random_result,fmt='%14.8lf\t%14.8lf\t%14.8lf')