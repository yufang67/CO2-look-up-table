import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors as colors

#######################LL############################################
data_LL = np.loadtxt("LL_p_T_c")
eLL=[]
vLL=[]
PLL=[]
eLL=data_LL[:,0]
vLL=data_LL[:,1]
PLL=data_LL[:,4]

col_LL = int(np.sqrt(len(eLL)))


p_LL = np.zeros((col_LL,col_LL)) 
e_LL = np.zeros((col_LL,col_LL))
v_LL = np.zeros((col_LL,col_LL))


e_LL = np.reshape(eLL,(col_LL,col_LL))
v_LL = np.reshape(vLL,(col_LL,col_LL))
p_LL = np.reshape(PLL,(col_LL,col_LL))
#p_LL = np.reshape(PLL,(col_LL,col_LL))
###############################LH####################################
data_LH = np.loadtxt("LH_p_T_c")
eLH=[]
vLH=[]
PLH=[]
eLH=data_LH[:,0]
vLH=data_LH[:,1]
PLH=data_LH[:,4]
#
col_LH = int(np.sqrt(len(eLH)))
#
#
P_LH = np.zeros((col_LH,col_LH)) 
e_LH = np.zeros((col_LH,col_LH))
v_LH = np.zeros((col_LH,col_LH))
#
e_LH = np.reshape(eLH,(col_LH,col_LH))
v_LH = np.reshape(vLH,(col_LH,col_LH))
#P_LH = np.reshape(PLH,(col_LH,col_LH))
p_LH = np.reshape(PLH,(col_LH,col_LH))
###############################HT#####################################
data_HT = np.loadtxt("HT_p_T_c")
eHT=[]
vHT=[]
PHT=[]
eHT=data_HT[:,0]
vHT=data_HT[:,1]
PHT=data_HT[:,4]

col_HT = int(np.sqrt(len(eHT)))


P_HT = np.zeros((col_HT,col_HT)) 
e_HT = np.zeros((col_HT,col_HT))
v_HT = np.zeros((col_HT,col_HT))


e_HT = np.reshape(eHT,(col_HT,col_HT))
v_HT = np.reshape(vHT,(col_HT,col_HT))
p_HT = np.reshape(PHT,(col_HT,col_HT))
########################R###########################################
data_R = np.loadtxt("R_p_T_c")
eR=[]
vR=[]
PR=[]
eR=data_R[:,0]
vR=data_R[:,1]
PR=data_R[:,4]

col_R = int(np.sqrt(len(eR)))


P_R = np.zeros((col_R,col_R)) 
e_R = np.zeros((col_R,col_R))
v_R = np.zeros((col_R,col_R))


e_R = np.reshape(eR,(col_R,col_R))
v_R = np.reshape(vR,(col_R,col_R))
p_R = np.reshape(PR,(col_R,col_R))
##########################two-phase#########################################
#data_TP = np.loadtxt("TP2_error-P-T-x.txt")
data_TP = np.loadtxt("TP2_T-x-vV-p-u-v.txt")
eTP=[]
vTP=[]
PTP=[]
eTP=data_TP[:,4]
vTP=data_TP[:,5]
PTP=data_TP[:,6]


col_TP = int(np.sqrt(len(eTP)))
p_TP = np.zeros((col_TP,col_TP)) 
e_TP = np.zeros((col_TP,col_TP))
v_TP = np.zeros((col_TP,col_TP))
e_TP = np.reshape(eTP,(col_TP,col_TP))
v_TP = np.reshape(vTP,(col_TP,col_TP))
p_TP = np.reshape(PTP,(col_TP,col_TP))

#%%
#norm=colors.PowerNorm(10)
plt.figure(figsize=(8,4))
#
plt.pcolormesh(v_LL, e_LL/1e3, p_LL, cmap='jet',shading='gouraud')
plt.clim(100,900)
plt.pcolormesh(v_LH, e_LH/1e3, p_LH, cmap='jet',shading='gouraud')
plt.clim(100,900)
plt.pcolormesh(v_R, e_R/1e3, p_R, cmap='jet',shading='gouraud')
plt.clim(100,900)
plt.pcolormesh(v_HT, e_HT/1e3, p_HT, cmap='jet',shading='gouraud')
plt.clim(100,900)
plt.pcolormesh(v_TP, e_TP/1e3, p_TP, cmap='jet',shading='gourand')
plt.clim(100,900)


####
ax = plt.gca() 
#ax.set_xlim([0.8e-3, 1.7e-1])
#ax.set_ylim([-4.5e2, 1e2])
#ax.set_xlim([1.9e-3, 2.2e-3])
#ax.set_ylim([-2.2e2, -1.8e2])
ax.tick_params(labelsize=5)
#plt.clim(1e-5,1e-1)
#plt.clim(1e-10,1e-1)
#plt.title('Pressure Reletive Errors',fontsize=17)
plt.ylabel('Internal Energy, e $[kJkg^{-1}]$',fontsize=15,position=(0,0.5),rotation = "vertical")
plt.xlabel('Specific Volume, v $[m^{-3}kg]$',fontsize=15,rotation = "horizontal")
plt.xticks(size = 12)
plt.yticks(size = 12)
ax.set_xscale('log')
plt.grid(True)
plt.colorbar(fraction=0.046,pad=0.01)
plt.text(3.3e-1,550, ' m/s', fontsize=11,fontweight='bold')

plt.tight_layout()
##,format='%.0e'
plt.savefig("sound.pdf")
plt.show()