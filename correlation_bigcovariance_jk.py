import numpy as np
import matplotlib.pyplot as plt
parallel_max = 5
parallel_min = 5
perpendicular_max = 28
subsample = 100
r = np.logspace(-1,1.7,28)

x = np.loadtxt('output')
data = x.reshape((parallel_max,subsample,perpendicular_max))
w = np.zeros([subsample,perpendicular_max*parallel_min])
for i in range(subsample):
	w[i] = data[0:parallel_min:1,i,:].reshape(perpendicular_max*parallel_min)
w_mean= np.mean(w,axis=0)
covariance=np.zeros([parallel_min*perpendicular_max,parallel_min*perpendicular_max])
bar = np.zeros(parallel_min*perpendicular_max)
for i in range(parallel_min*perpendicular_max):
	for j in range(parallel_min*perpendicular_max):
		for k in range(subsample):
			covariance[i][j] = covariance[i][j] + 1.0*(subsample-1)/subsample * (w[k][i] - w_mean[i])*(w[k][j] - w_mean[j])
for i in range(parallel_min*perpendicular_max):
	bar[i] = np.sqrt(covariance[i][i])
plt.errorbar(r, w_mean[0*perpendicular_max:1*perpendicular_max], yerr=bar[0*perpendicular_max:1*perpendicular_max], capsize=2, marker='s', markersize=1, color='r', label='emission100')
plt.errorbar(r, w_mean[1*perpendicular_max:2*perpendicular_max], yerr=bar[1*perpendicular_max:2*perpendicular_max], capsize=2, marker='s', markersize=1, color='g', label='emission200')
plt.errorbar(r, w_mean[2*perpendicular_max:3*perpendicular_max], yerr=bar[2*perpendicular_max:3*perpendicular_max], capsize=2, marker='s', markersize=1, color='b', label='emission300')
plt.errorbar(r, w_mean[3*perpendicular_max:4*perpendicular_max], yerr=bar[3*perpendicular_max:4*perpendicular_max], capsize=2, marker='s', markersize=1, color='m', label='emission400')
plt.errorbar(r, w_mean[4*perpendicular_max:5*perpendicular_max], yerr=bar[4*perpendicular_max:5*perpendicular_max], capsize=2, marker='s', markersize=1, color='c', label='emission500')
'''
x1 = np.loadtxt('../whole/output')
data1 = x1.reshape((parallel_max,subsample,perpendicular_max))
w1 = np.zeros([subsample,perpendicular_max*parallel_min])
for i in range(subsample):
	w1[i] = data1[0:parallel_min:1,i,:].reshape(perpendicular_max*parallel_min)
w_mean1= np.mean(w1,axis=0)
covariance1=np.zeros([parallel_min*perpendicular_max,parallel_min*perpendicular_max])
bar1 = np.zeros(parallel_min*perpendicular_max)
for i in range(parallel_min*perpendicular_max):
	for j in range(parallel_min*perpendicular_max):
		for k in range(subsample):
			covariance1[i][j] = covariance1[i][j] + 1.0*(subsample-1)/subsample * (w1[k][i] - w_mean1[i])*(w1[k][j] - w_mean1[j])
for i in range(parallel_min*perpendicular_max):
	bar1[i] = np.sqrt(covariance1[i][i])
plt.errorbar(r, w_mean1[0*perpendicular_max:1*perpendicular_max], yerr=bar1[0*perpendicular_max:1*perpendicular_max], capsize=2, fmt='--', markersize=1, color='r', label='whole100')
plt.errorbar(r, w_mean1[1*perpendicular_max:2*perpendicular_max], yerr=bar1[1*perpendicular_max:2*perpendicular_max], capsize=2, fmt='--', markersize=1, color='g', label='whole200')
plt.errorbar(r, w_mean1[2*perpendicular_max:3*perpendicular_max], yerr=bar1[2*perpendicular_max:3*perpendicular_max], capsize=2, fmt='--', markersize=1, color='b', label='whole300')
#plt.errorbar(r, w_mean1[3*perpendicular_max:4*perpendicular_max], yerr=bar1[3*perpendicular_max:4*perpendicular_max], capsize=2, fmt='--', markersize=1, color='c', label='whole400')
#plt.errorbar(r, w_mean1[4*perpendicular_max:5*perpendicular_max], yerr=bar1[4*perpendicular_max:5*perpendicular_max], capsize=2, fmt='--', markersize=1, color='r', label='whole500')


x2=np.loadtxt('../test2000000/w_mean')
data2 = x2.reshape((parallel_max,perpendicular_max))
plt.plot(r,data2[0])
plt.plot(r,data2[1])
plt.plot(r,data2[2])
plt.plot(r,data2[3])
plt.plot(r,data2[4])
'''
plt.title('correlation functions')
plt.xlabel("$r_{p}$(Mpc/h)",size="xx-large")
plt.ylabel('$w_{p}$',size="xx-large")
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.legend()
plt.show()
#plt.savefig("compare1",format="pdf")
np.savetxt('w_mean',w_mean,'%12.6f')
np.savetxt("covariance_matrix",covariance,"%12.7f")