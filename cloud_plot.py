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
import matplotlib.cm as cm
matplotlib.use('agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 10})
cwd=os.getcwd()
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

#plt.rcParams.update({'font.size': 22})
#plt.rc('xlabel', labelsize=20) 
#plt.rc('ylabel', labelsize=20) 
#plt.rc('title', labelsize=20) 

# path1="/home/shann87/shann87/fortran/Data/3_1/0001x0001/cal0p1/"
# path2="0500_0p50"
# path2="1000_1p20"
# path2="020000_2p00"

path1="/home/shann87/shann87/fortran/Data/3_1/0001x0001/cloud_folder_temp/001000_1p20"
path2=""
#path2="V100kap0p7.txt"
stringkg="001000_1p20"
[V, k_over_g]=stringkg.split('_')
str_V=V
str_k=k_over_g
str_k.replace('.','p')
V=int(V)
k_over_g=float(k_over_g.replace('p','.'))
#V=200
#k_over_g=2.2
# print(V)
# print(k_over_g)
kap=k_over_g
rho=0.1
gam=1
fp2_dat=rho*V/(gam+kap)
fp1_dat=rho*V/gam
fp6_1_dat=rho*V/gam
fp6_2_dat=rho*(gam-kap)*V/gam**2
##########################################################################################################

###########################333               Read data     ##########################################3
#########################################################################################################
#print(path1+path2)
os.chdir(path1+path2)

#Table:000500_0p50 ...0000002_00_cloud.dat  s=160
#		001000_1p20  ....0000006_00_cloud.dat	s=60 df[time]/1000 K
#		001000_1p80 .....0000008_00_cloud.dat  s=30

#list=os.listdir(path1+path2)
df=[]
#for timestamp in list:
#	os.chdir(path1+path2+timestamp)
for infile in glob.glob("*.dat"):
	#datafiles=glob('*.dat')
	if infile=='0000006_00_cloud.dat':
		df_new = pd.read_csv(infile,sep='\s+',names = ["ID", "extinctions", "pA","pB","pC","time"])
		# #print(df_new) 
		# df.append(df_new)
		# #df.append(df_new.loc[df_new['extinctions'] == 2])
	# df=pd.concat(df_new, axis=0)	
	#del df_new
# #df_final=[]
df=df_new
print(df["ID"].nunique())

n_rows=df.shape[0]
df_FP2=df.iloc[0:int(2.0/3)*n_rows]
print(df_FP2.shape)
# result=df.loc[df["extinctions"] == 2]
# #df["time"].mean()
# #df.columns = ["ID", "extinctions", "pA","pB","pC","time"]
# #print(df)
# #print(df_new)
# dataset=(result["time"]-result["time"].mean())/(result["time"].std())
# # print(dataset.std())
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



##################################################################################################333
######################3       COMPUTE CLOUD PLOT      3######################################
# An "interface" to matplotlib.axes.Axes.hist() method
# n, bins, patches = plt.hist(x=df["time"]/100000, bins=10, color='#0504aa',alpha=0.7, rwidth=0.85)# range=[0, 30000000],
fig = plt.figure(figsize=(10, 9))

ax = fig.add_subplot(111, projection='3d')
ax.scatter(fp2_dat, fp2_dat, fp2_dat, s=660 ,alpha=1.0, c='green', marker='^')
# ax.scatter(fp6_2_dat, 0, fp6_1_dat, s=660 ,alpha=1.0, c='blue', marker='^')


p=ax.scatter(df["pB"], df["pA"], df["pC"], s=40 ,alpha=0.5, c=df['time']/1000, cmap='viridis_r')
ax.scatter(0, 0, fp1_dat, s=660 ,alpha=1.0, c='red', marker='^')

ax.view_init(elev=45,azim=45)

# print(bins)
# plt.grid(axis='y', alpha=0.75)
# map1 = ax.imshow(np.stack([df['time'], df['time']]),cmap='viridis')

ax.set_xlabel(r'$N_{a}$', fontsize=25,labelpad=20)
ax.set_ylabel(r'$N_{b}$', fontsize=25,labelpad=20)
ax.zaxis.set_rotate_label(False) 
ax.set_zlabel(r'$N_{c}$', fontsize=25,labelpad=10)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='z', labelsize=20)
position=fig.add_axes([0.88,0.15,0.02,0.7]) #[left,bottom, width, height]
cbar=fig.colorbar(p,cax=position)#orientation="horizontal"
cbar.set_label('time (K)', fontsize=20)
cbar.ax.tick_params(labelsize=20)


# #plt.title('Probability distribution of Mean Time Extinction')
# plt.text(0.6, 0.8, r'$\gamma=1, \kappa=%1.2f, \rho=0.1,$' '\n' '$V=%4d$' % (k_over_g, V), horizontalalignment='center', verticalalignment='center',transform=plt.gca().transAxes, fontsize=20)
# maxfreq = n.max()
# # Set a clean upper y-axis limit.
# plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
os.chdir(cwd)
plt.savefig('Nov1_cloud_plot_k_%s_V_%s.png' % (str_k, str_V),bbox_inches='tight')

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



