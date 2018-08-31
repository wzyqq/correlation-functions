import numpy as np
import matplotlib.pyplot as plt
subsample = 100
subsample_axis = 10
position=np.loadtxt('emission random sample')
data = np.loadtxt('select emission')
plt.plot(position[:,0],position[:,1],'r.')
random_location = np.int_(np.linspace(0,len(position)-1,len(position)))
data_location   = np.int_(np.linspace(0,len(data)-1    ,len(data)    ))

x_bound = np.zeros(subsample_axis+1)
y_bound = np.zeros([subsample_axis,subsample_axis+1])
xdata   = np.zeros(subsample_axis)
data_order     = data[np.lexsort(data[:,::-1].T)]
position_order = position[np.lexsort(position[:,::-1].T)]
for i in range(subsample_axis+1):
	x_bound[i] =  np.percentile(position_order[:,0],100.0/subsample_axis*i)


random_loc0_list = []
random_loc1_list = []
random_temp      = []
random_xdata     = []
data_loc0_list = []
data_loc1_list = []
data_temp      = []
data_xdata     = []

for i in range(subsample_axis):
	random_loc0_list.append( np.where(position_order[:,0]<x_bound[i]))
	random_loc1_list.append( np.where(position_order[:,0]>=x_bound[i+1]))
	random_temp.append( np.hstack((random_loc0_list[i],random_loc1_list[i])))
	random_xdata.append( position_order[np.delete(random_location,random_temp[i])])
	data_loc0_list.append( np.where(data_order[:,0]<x_bound[i]))
	data_loc1_list.append( np.where(data_order[:,0]>=x_bound[i+1]))
	data_temp.append( np.hstack((data_loc0_list[i],data_loc1_list[i])))
	data_xdata.append( data_order[np.delete(data_location,data_temp[i])])

random_xdata[subsample_axis-1] = np.append(random_xdata[subsample_axis-1],position_order[np.where(position_order[:,0]==x_bound[subsample_axis])],axis=0)
data_xdata[subsample_axis-1]   = np.append(data_xdata[subsample_axis-1],data_order[np.where(data_order[:,0]==x_bound[subsample_axis])],axis=0)

for i in range(subsample_axis):
	for j in np.int_(np.linspace(0,subsample_axis,subsample_axis+1)):
		y_bound[i][j] =  np.percentile(random_xdata[i][:,1],100.0/subsample_axis*j)

data_loc3_list   = []
data_loc4_list   = []
random_loc3_list = []
random_loc4_list = []
middata          = []
mid_dataorder    = []
midafterdata     = []
midrandom        = []
mid_randomorder  = []
midafterrandom   = []
dataaftertemp    = []
randomaftertemp  = []
finaldata        = []
finalrandom      = []

for i in range(subsample_axis):
	middata.append(np.array([data_xdata[i][:,1],data_xdata[i][:,0],data_xdata[i][:,2]]).T)
	mid_dataorder.append(middata[i][np.lexsort(middata[i][:,::-1].T)])
	midafterdata.append(np.array([mid_dataorder[i][:,1],mid_dataorder[i][:,0],mid_dataorder[i][:,2]]).T)
	midrandom.append(np.array([random_xdata[i][:,1],random_xdata[i][:,0],random_xdata[i][:,2]]).T)
	mid_randomorder.append(midrandom[i][np.lexsort(midrandom[i][:,::-1].T)])
	midafterrandom.append(np.array([mid_randomorder[i][:,1],mid_randomorder[i][:,0],mid_randomorder[i][:,2]]).T)

dataafterlocation=[]
randomafterlocation=[]

for i in range(subsample_axis):
	dataafterlocation.append(np.int_(np.linspace(0,len(midafterdata[i])-1,len(midafterdata[i]))))
	randomafterlocation.append(np.int_(np.linspace(0,len(midafterrandom[i])-1,len(midafterrandom[i]))))

for i in range(subsample_axis):
	for j in range(subsample_axis):
		data_loc3_list.append( np.where(midafterdata[i][:,1]<y_bound[i][j]))
		data_loc4_list.append( np.where(midafterdata[i][:,1]>=y_bound[i][j+1]))
		dataaftertemp.append( np.hstack((data_loc3_list[subsample_axis*i+j],data_loc4_list[subsample_axis*i+j])))
		finaldata.append( midafterdata[i][np.delete(dataafterlocation[i],dataaftertemp[subsample_axis*i+j])])
		random_loc3_list.append( np.where(midafterrandom[i][:,1]<y_bound[i][j]))
		random_loc4_list.append( np.where(midafterrandom[i][:,1]>=y_bound[i][j+1]))
		randomaftertemp.append( np.hstack((random_loc3_list[subsample_axis*i+j],random_loc4_list[subsample_axis*i+j])))
		finalrandom.append( midafterrandom[i][np.delete(randomafterlocation[i],randomaftertemp[subsample_axis*i+j])])
	finaldata[subsample_axis*i+subsample_axis-1] = np.append(finaldata[subsample_axis*i+subsample_axis-1],
		midafterdata[i][np.where(midafterdata[i][:,1]==y_bound[i][subsample_axis])],axis=0)
	finalrandom[subsample_axis*i+subsample_axis-1] = np.append(finalrandom[subsample_axis*i+subsample_axis-1],
		midafterrandom[i][np.where(midafterrandom[i][:,1]==y_bound[i][subsample_axis])],axis=0)

for i in range(subsample):

	datafilename = "emission_sample/emission_sample_%d"%(i)
	randfilename = "emission_random_sample/emission_random_sample_%d"%(i)	
	np.savetxt(datafilename, finaldata[i],fmt='%14.8f')
	np.savetxt(randfilename, finalrandom[i],fmt='%14.8f')
