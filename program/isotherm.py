import numpy as np
import matplotlib.pyplot as plt


sp_LL=np.loadtxt('spline_LL.txt')
sp_LH=np.loadtxt('spline_LH.txt')
sp_HT=np.loadtxt('spline_HT.txt')
sp_R=np.loadtxt('spline_R.txt')
##
AA=np.loadtxt('HT_p_T_c')
BB=np.loadtxt('LH_p_T_c')
CC=np.loadtxt('R_p_T_c')
DD=np.loadtxt('LL_p_T_c')
EE=np.loadtxt('TP_T-x-Perr-p-u-v.txt')

###%%% Top %%%
#Top = np.zeros((10,2))
#Top[:,0] = np.linspace(0.0088,0.1706,num=10);
#Top[:,1] = 68569;

u_HT_A = AA[:,0]
v_HT_A = AA[:,1]
#p_HT_A = AA[:,2]
T_HT_A = AA[:,3]
t_HT_A = AA[:,3]


u_LH_B = BB[:,0]
v_LH_B = BB[:,1]
#p_LH_B = BB[:,2]
T_LH_B = BB[:,3];

u_R_C = CC[:,0];
v_R_C = CC[:,1];
#p_R_C = CC(:,3);
T_R_C = CC[:,3];

u_LL_D = DD[:,0];
v_LL_D = DD[:,1];
#p_LL_D = DD(:,3);
T_LL_D = DD[:,3];

u_TP_E = EE[:,4]
v_TP_E = EE[:,5]
#p_TP_E = EE[:,3]
T_TP_E = EE[:,0]
   
###########################isobar single phase 8,5,2,1######################################
u_HT_isoT400 = []
v_HT_isoT400 = []
u_HT_isoT350 = []
v_HT_isoT350 = []
u_HT_isoT300 = []
v_HT_isoT300 = []
u_HT_isoT230 = []
v_HT_isoT230 = []
u_HT_isoT250 = []
v_HT_isoT250 = []
#
u_LH_isoT400 = []
v_LH_isoT400 = []
u_LH_isoT350 = []
v_LH_isoT350 = []
u_LH_isoT300 = []
v_LH_isoT300 = []
#
u_TP_isoT2166 = []
v_TP_isoT2166 = []
u_TP_isoT230 = []
v_TP_isoT230 = []
u_TP_isoT250 = []
v_TP_isoT250 = []
u_TP_isoT300 = []
v_TP_isoT300 = []
#
u_R_isoT230 = []
v_R_isoT230 = []
#
u_LL_isoT350 = []
v_LL_isoT350 = []
u_LL_isoT300 = []
v_LL_isoT300 = []
u_LL_isoT250 = []
v_LL_isoT250 = []
u_LL_isoT230 = []
v_LL_isoT230 = []
#
for i in range(len(u_HT_A)):
    if T_HT_A[i] <=6.03e2 and T_HT_A[i]>=5.99e2:
        u_HT_isoT400.append(u_HT_A[i]/1000)
        v_HT_isoT400.append(v_HT_A[i])
    if T_HT_A[i] <=3.51e2 and T_HT_A[i]>=3.48e2:
        u_HT_isoT350.append(u_HT_A[i]/1000)
        v_HT_isoT350.append(v_HT_A[i])
    if T_HT_A[i] <=3.01e2 and T_HT_A[i]>=2.95e2:
        u_HT_isoT300.append(u_HT_A[i]/1000)
        v_HT_isoT300.append(v_HT_A[i])
    if T_HT_A[i] <=2.31e2 and T_HT_A[i]>=2.21e2:
        u_HT_isoT230.append(u_HT_A[i]/1000)
        v_HT_isoT230.append(v_HT_A[i])
    if T_HT_A[i] <=2.51e2 and T_HT_A[i]>=2.45e2:
        u_HT_isoT250.append(u_HT_A[i]/1000)
        v_HT_isoT250.append(v_HT_A[i])
#
for i in range(len(u_LH_B)):
#    if T_LH_B[i] <=4.01e2 and T_LH_B[i]>=3.99e2:
#        u_LH_isoT400.append(u_LH_B[i]/1000)
#        v_LH_isoT400.append(v_LH_B[i])
    if T_LH_B[i] <=3.51e2 and T_LH_B[i]>=3.497e2:
        u_LH_isoT350.append(u_LH_B[i]/1000)
        v_LH_isoT350.append(v_LH_B[i])
    if T_LH_B[i] <=3.01e2 and T_LH_B[i]>=2.999e2:
        u_LH_isoT300.append(u_LH_B[i]/1000)
        v_LH_isoT300.append(v_LH_B[i])
#############################isobar two-phase####################################
for i in range(len(u_TP_E)):
    if T_TP_E[i] <=2.167e2 and T_TP_E[i]>=2.164e2:
#        print(T_TP_E[i])
        u_TP_isoT2166.append(u_TP_E[i]/1000)
        v_TP_isoT2166.append(v_TP_E[i])
    if T_TP_E[i] <=2.3001e2 and T_TP_E[i]>=2.298e2:
        u_TP_isoT230.append(u_TP_E[i]/1000)
        v_TP_isoT230.append(v_TP_E[i])
    if T_TP_E[i] <=2.5001e2 and T_TP_E[i]>=2.499e2:
        u_TP_isoT250.append(u_TP_E[i]/1000)
        v_TP_isoT250.append(v_TP_E[i])
    if T_TP_E[i] <=3.001e2 and T_TP_E[i]>=2.999e2:
        u_TP_isoT300.append(u_TP_E[i]/1000)
        v_TP_isoT300.append(v_TP_E[i])
#################################R#########################################
for i in range(len(u_R_C)):
    if T_R_C[i] <=2.301e2 and T_R_C[i]>=2.299e2:
        u_R_isoT230.append(u_R_C[i]/1000)
        v_R_isoT230.append(v_R_C[i])
#########################################################################
for i in range(len(u_LL_D)):
    if T_LL_D[i] <=2.301e2 and T_LL_D[i]>=2.299e2:
        u_LL_isoT230.append(u_LL_D[i]/1000)
        v_LL_isoT230.append(v_LL_D[i])
    if T_LL_D[i] <=2.5001e2 and T_LL_D[i]>=2.499e2:
        u_LL_isoT250.append(u_LL_D[i]/1000)
        v_LL_isoT250.append(v_LL_D[i])
    if T_LL_D[i] <=3.001e2 and T_LL_D[i]>=2.999e2:
        u_LL_isoT300.append(u_LL_D[i]/1000)
        v_LL_isoT300.append(v_LL_D[i])
    if T_LL_D[i] <=3.51e2 and T_LL_D[i]>=3.497e2:
        u_LL_isoT350.append(u_LL_D[i]/1000)
        v_LL_isoT350.append(v_LL_D[i])
########################################################################
    
    
plt.figure(figsize=(9,4)) 
#
plt.semilogx(sp_R[:,0],sp_R[:,2]/1e3,color='r',linewidth=2)
plt.semilogx(sp_R[:,1],sp_R[:,2]/1e3,color='k',linewidth=2)

plt.semilogx(sp_HT[:,0],sp_HT[:,2]/1e3,linestyle='--',color='k',linewidth=2)
plt.semilogx(sp_HT[:,1],sp_HT[:,2]/1e3,linestyle='--',color='k',linewidth=2)

plt.semilogx(sp_LH[:,0],sp_LH[:,2]/1e3,color='r',linewidth=2)
plt.semilogx(sp_LH[:,1],sp_LH[:,2]/1e3,linestyle='--',color='k',linewidth=2)

plt.semilogx(sp_LL[:,0],sp_LL[:,2]/1e3,linestyle='-',color='b',linewidth=2)
plt.semilogx(sp_LL[:,1],sp_LL[:,2]/1e3,linestyle='--',color='k',linewidth=2)
###############################################################################
plt.semilogx(v_HT_isoT400,u_HT_isoT400,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_HT_isoT350,u_HT_isoT350,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_HT_isoT300,u_HT_isoT300,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_HT_isoT250,u_HT_isoT250,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_HT_isoT230,u_HT_isoT230,color='k',linewidth=2,linestyle=':')
#
plt.semilogx(v_LH_isoT400,u_LH_isoT400,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_LH_isoT350,u_LH_isoT350,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_LH_isoT300,u_LH_isoT300,color='k',linewidth=2,linestyle=':')
#
plt.semilogx(v_R_isoT230,u_R_isoT230,color='k',linewidth=2,linestyle=':')
#
plt.semilogx(v_LL_isoT350,u_LL_isoT350,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_LL_isoT300,u_LL_isoT300,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_LL_isoT250,u_LL_isoT250,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_LL_isoT230,u_LL_isoT230,color='k',linewidth=2,linestyle=':')
#########################TWO-PHASE########################################
plt.semilogx(v_TP_isoT2166,u_TP_isoT2166,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_TP_isoT230,u_TP_isoT230,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_TP_isoT250,u_TP_isoT250,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_TP_isoT300,u_TP_isoT300,color='k',linewidth=2,linestyle=':')

#####control des labels
ax = plt.gca() 
#ax.set_xlim([0.015, 1e-1])
ax.set_ylim([-4.5e2, 3e2])
#plt.title('Isothermal curves',fontsize=17)
plt.ylabel('Internal Energy, e $[kJkg^{-1}]$',fontsize=15,position=(0,0.5),rotation = "vertical")
plt.xlabel('Specific Volume, v $[m^{-3}kg]$',fontsize=15,rotation = "horizontal")
plt.xticks(size = 12)
plt.yticks(size = 12)
plt.grid(True)
#plt.axis([-4,4,-0.3,1.5])
#plt.xlabel('X', color='C1')
#plt.ylabel('X', color='0.5')  # grayscale color
#ax.xaxis.set_label_coords(0.5,-0.05)
#ax.yaxis.set_label_coords(-0.08,0.5)
#####################################
ax.text(2.0e-2,200, '$600K$', fontsize=10)
ax.text(2.1e-2,-25, '$350K$', fontsize=10)
ax.text(2.8e-3,-200, '$300K$', fontsize=10)
ax.text(0.9e-2,-210, '$250K$', fontsize=10)
ax.text(2.0e-2,-210, '$230K$', fontsize=10)
ax.text(3.3e-2,-210, '$216.6K$', fontsize=10)
#############################################################################
plt.tight_layout()
plt.savefig("isotherm.pdf")
plt.show()
