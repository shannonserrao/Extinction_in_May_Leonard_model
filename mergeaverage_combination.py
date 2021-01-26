#
#print('Hello, we start the averaging process of the extinction times')
#folder=input("enter the folder name")
#V=input("The Volume parameter")
#kg_ratio=input("the kappa/gamma ratio")
#er=input("Is it extinction or relaxation")
#path2="V"+str(V)
#path3="/"+folder

#for i in range(1, eof):
	
#folder=input("enter the folder name")
import pandas as pd
import os
import numpy as np
import csv
import glob
from scipy import stats
#from glob import glob
from fnmatch import fnmatch
import matplotlib	
matplotlib.use('agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 15})
cwd=os.getcwd()
#plt.rcParams.update({'font.size': 22})
#plt.rc('xlabel', labelsize=20) 
#plt.rc('ylabel', labelsize=20) 
#plt.rc('title', labelsize=20) 

path1="/home/shann87/shann87/fortran/Data/3_1/0001x0001/cal0p1/"
path2="0500_0p50"
# paths2=["1000_0p90","1000_1p00","1000_1p20","1000_1p50","1000_1p80","001000_1p90","1000_2p00"]
# paths2=["1000_1p00","1000_1p20","1000_1p50","1000_1p80","001000_1p90","1000_2p00"]
paths2=["1000_1p50","1000_1p80","1000_2p00"]
dflist=[]
sizelist=[]
weightlist=[]
for path2 in paths2:
	# path2="1000_1p20"
	# path2="020000_2p00"

	#path2="V100kap0p7.txt"
	[V, k_over_g]=path2.split('_')
	str_V=V
	str_k=k_over_g
	str_k.replace('.','p')
	V=int(V)
	k_over_g=float(k_over_g.replace('p','.'))
	#V=200
	#k_over_g=2.2
	# print(V)
	# print(k_over_g)

	##########################################################################################################

	###########################333               Read data     ##########################################3
	#########################################################################################################
	#print(path1+path2)
	os.chdir(path1+path2)

	#list=os.listdir(path1+path2)
	df=[]
	
	#for timestamp in list:
	#	os.chdir(path1+path2+timestamp)
	for infile in glob.glob("*.dat"):
		#datafiles=glob('*.dat')
		df_new = pd.read_csv(infile,sep='\s+',names = ["ID", "extinctions", "pA","pB","pC","time"])
		#print(df_new) 
		df.append(df_new)
		#df.append(df_new.loc[df_new['extinctions'] == 2])
	df=pd.concat(df, axis=0)	
		#del df_new
	#df_final=[]
	result=df.loc[df["extinctions"] == 2]
	weights = np.ones_like(result['time']) / result.shape[0]
	sizelist.append(result.shape[0])
	weightlist.append(weights)
	dflist.append(result['time'])
	# print(type(dflist))

	#df["time"].mean()
	#df.columns = ["ID", "extinctions", "pA","pB","pC","time"]
	#print(df)
	#print(df_new)
	# dataset=(result["time"]-result["time"].mean())/(result["time"].std())
	# print(dataset.std())
	# ktest=stats.kurtosistest(dataset, axis=0, nan_policy='propagate')
	# normtest=stats.normaltest(dataset, axis=0, nan_policy='propagate')
	# kurto=stats.kurtosis(dataset)

	# print('k/g is '+str_k)
	# print('V is' +str_V)
	# print('mean is ' + str(result["time"].mean()))
	# print('std is ' + str(result["time"].std()))

	# print('kurtosis is  ' + str(kurto))

	# print('Results of norm test :' + str(normtest))

	# print('Results of k test :' + str(ktest))
# print(type(dflist))
# print(len(dflist))
print(sizelist)
##################################################################################################333
######################3       COMPUTE HISTOGRAMS      3######################################
# An "interface" to matplotlib.axes.Axes.hist() method
# x_multi = [np.random.randn(n) for n in [10000, 5000, 2000, 500, 300, 200, 50]]
# print(x_multi[0].shape)

plt.hist(dflist, bins=20,range=[100,8000], histtype='bar',weights=weightlist, alpha=1.0, rwidth=0.8)# range=[0, 30000000],
# print(bins)
# plt.legend((r'$\kappa= 1.5$',r'$\kappa= 1.8$',r'$\kappa= 1.9$',r'$\kappa= 2.0$'))
plt.legend((r'$\kappa= 1.5$',r'$\kappa= 1.8$',r'$\kappa= 2.0$'))
plt.grid(axis='y', alpha=0.75)
plt.grid(axis='x', alpha=0.75)
plt.xlabel('TE ', fontsize=15)
plt.ylabel('pdf', fontsize=15)
#plt.title('Probability distribution of Mean Time Extinction')
# plt.text(0.6, 0.8, r'$\gamma=1, \kappa=%1.2f, \rho=0.1,$' '\n' '$V=%4d$' % (k_over_g, V), horizontalalignment='center', verticalalignment='center',transform=plt.gca().transAxes, fontsize=20)
# maxfreq = n.max()
# Set a clean upper y-axis limit.
# plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
os.chdir(cwd)
# plt.savefig('Aug26_hist_plot_k_%s_V_%s.png' % (str_k, str_V))
plt.savefig('Nov5_multi_hist_plot.png')

######################################################################################################3


# import sys
# sys.stdout=open('/home//shann87/fortran/extinctiontimes.dat', 'w+')
# #f=open('/home//shann87/fortran/extinctiontimes.dat', 'w+')
# print "The extinction time for V ", V, "for kap/gamma ", k_over_g, "ext time =", result["time"].mean()
# #f.write(displayvar)
# sys.stdout.close()	




################################################################################################################
#####################3				Rough				
#os.chdir(path1+path2)
#list=os.listdir(path1+path2)
#os.chdir(path1+path2+list[1])
#datafiles=glob('*.dat')
#print(datafiles)
#data = pd.read_csv(datafiles[0],delimiter=' ',names = ["Sequence", "Start", "End", "Coverage"])
#df.loc[df['column_name'] == some_value]
#df["weight"].mean()
#data.head()
#print(data.describe())



