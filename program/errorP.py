import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#######################LL############################################
data_LL = np.loadtxt("err_pBL_LL.txt")
eLL=[]
vLL=[]
PLL=[]
eLL=data_LL[:,0]
vLL=data_LL[:,1]
PLL=data_LL[:,2]

col_LL = int(np.sqrt(len(eLL)))


p_LL = np.zeros((col_LL,col_LL)) 
e_LL = np.zeros((col_LL,col_LL))
v_LL = np.zeros((col_LL,col_LL))


e_LL = np.reshape(eLL,(col_LL,col_LL))
v_LL = np.reshape(vLL,(col_LL,col_LL))
p_LL = np.reshape(PLL,(col_LL,col_LL))
#p_LL = np.reshape(PLL,(col_LL,col_LL))
###############################LH####################################
data_LH = np.loadtxt("err_pBL_LH.txt")
eLH=[]
vLH=[]
PLH=[]
eLH=data_LH[:,0]
vLH=data_LH[:,1]
PLH=data_LH[:,2]
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
data_HT = np.loadtxt("err_pBL_HT.txt")
eHT=[]
vHT=[]
PHT=[]
eHT=data_HT[:,0]
vHT=data_HT[:,1]
PHT=data_HT[:,2]

col_HT = int(np.sqrt(len(eHT)))


P_HT = np.zeros((col_HT,col_HT)) 
e_HT = np.zeros((col_HT,col_HT))
v_HT = np.zeros((col_HT,col_HT))


e_HT = np.reshape(eHT,(col_HT,col_HT))
v_HT = np.reshape(vHT,(col_HT,col_HT))
p_HT = np.reshape(PHT,(col_HT,col_HT))
########################R###########################################
data_R = np.loadtxt("err_pBL_R.txt")
eR=[]
vR=[]
PR=[]
eR=data_R[:,0]
vR=data_R[:,1]
PR=data_R[:,2]

col_R = int(np.sqrt(len(eR)))


P_R = np.zeros((col_R,col_R)) 
e_R = np.zeros((col_R,col_R))
v_R = np.zeros((col_R,col_R))


e_R = np.reshape(eR,(col_R,col_R))
v_R = np.reshape(vR,(col_R,col_R))
p_R = np.reshape(PR,(col_R,col_R))
##########################two-phase#########################################
data_TP = np.loadtxt("TP2_error-P-T-x.txt")
data_TP2 = np.loadtxt("TP2_T-x-vV-p-u-v.txt")
eTP=[]
vTP=[]
PTP=[]
eTP=data_TP[:,3]
vTP=data_TP[:,4]
PTP=data_TP[:,0]

eTP2=[]
vTP2=[]
vari2=[]
eTP2=data_TP2[:,4]
vTP2=data_TP2[:,5]
vari2=data_TP2[:,6]

col_TP = int(np.sqrt(len(eTP)))
col_TP2 = int(np.sqrt(len(eTP2)))
p_TP = np.zeros((col_TP,col_TP)) 
e_TP = np.zeros((col_TP,col_TP))
v_TP = np.zeros((col_TP,col_TP))
e_TP = np.reshape(eTP,(col_TP,col_TP))
v_TP = np.reshape(vTP,(col_TP,col_TP))
p_TP = np.reshape(PTP,(col_TP,col_TP))

e_TP2 = np.reshape(eTP2,(col_TP2,col_TP2))
v_TP2 = np.reshape(vTP2,(col_TP2,col_TP2))
vari_TP = np.reshape(vari2,(col_TP2,col_TP2))
#%%
#norm=colors.PowerNorm(10)
plt.figure(figsize=(8,4))
#
#plt.pcolormesh(v_R, e_R/1e3, P_R*1e4, cmap='jet',shading='gouraud')
#plt.pcolormesh(v_LL, e_LL/1e3, P_LL*1e4, cmap='jet',shading='gouraud')
#plt.pcolormesh(v_LH, e_LH/1e3, P_LH*1e4, cmap='jet',shading='gouraud')
#plt.pcolormesh(v_TP, e_TP/1e3, P_TP*1e4, cmap='jet',shading='gourand')
#plt.pcolormesh(v_HT, e_HT/1e3, P_HT*1e4, cmap='jet',shading='gouraud')

##log
plt.pcolormesh(v_LL, e_LL/1e3, p_LL, cmap='jet',shading='gouraud',norm=colors.LogNorm(vmin=p_LL.min(),vmax=p_LL.max()))
plt.clim(1e-7,1e-1)
plt.pcolormesh(v_R, e_R/1e3, p_R, cmap='jet',shading='gouraud',norm=colors.LogNorm())
plt.clim(1e-7,1e-1)
plt.pcolormesh(v_LH, e_LH/1e3, p_LH, cmap='jet',shading='gouraud',norm=colors.LogNorm())
plt.clim(1e-7,1e-1)
plt.pcolormesh(v_TP, e_TP/1e3, p_TP, cmap='jet',norm=colors.LogNorm(vmin=p_TP.min(),vmax=p_TP.max()))
plt.clim(1e-7,1e-2)
plt.pcolormesh(v_HT, e_HT/1e3, p_HT, cmap='jet',norm=colors.LogNorm(vmin=P_HT.min(),vmax=P_HT.max()))
plt.clim(1e-7,1e-1)

#plt.pcolormesh(v_TP2, e_TP2/1e3, vari_TP, cmap='jet')
#plt.clim(10,500)




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

plt.tight_layout()
##,format='%.0e'
plt.savefig("log_pmap.pdf")
plt.show()