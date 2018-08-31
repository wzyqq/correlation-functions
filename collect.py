import numpy as np
subsample = 100
observe=[]
random=[]
observeall=[]
randomall=[]
for i in range(subsample):
	observe_position  = "emission_sample_position/emission_sample_position_%d"%(i)
	random_position = "emission_random_sample_position/emission_random_sample_position_%d"%(i)	
	observe_length  = "emission_sample_position_length/emission_sample_position_length_%d"%(i)
	random_length = "emission_random_sample_position_length/emission_random_sample_position_length_%d"%(i)	
	observe = np.loadtxt(observe_position)
	random  = np.loadtxt(random_position)
	observelength = np.array([len(observe)])
	randomlength  = np.array([len(random)])
	for k in range(len(observe)):
		observeall.append(observe[k])
	for k in range(len(random)):
		randomall.append(random[k])
	np.savetxt(observe_length,observelength,fmt='%10d')
	np.savetxt(random_length,randomlength,fmt='%10d')
np.savetxt("emission_sample_position_all",observeall,fmt='%14.8lf\t%14.8lf\t%14.8lf')
np.savetxt("emission_random_sample_position_all",randomall,fmt='%14.8lf\t%14.8lf\t%14.8lf')
