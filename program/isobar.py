import numpy as np
import matplotlib.pyplot as plt


sp_LL=np.loadtxt('spline_LL.txt')
sp_LH=np.loadtxt('spline_LH.txt')
sp_HT=np.loadtxt('spline_HT.txt')
sp_R=np.loadtxt('spline_R.txt')
##
AA=np.loadtxt('HT_p_T_c')
BB=np.loadtxt('LH_p_T_c')
#CC=np.loadtxt('R_p_T_c')
DD=np.loadtxt('LL_p_T_c')
EE=np.loadtxt('TP_T-x-Perr-p-u-v.txt')

u_HT_A = AA[:,0]
v_HT_A = AA[:,1]
p_HT_A = AA[:,2]
#T_HT_A = AA[:,4]
#c_HT_A = AA[:,5]

u_LH_B = BB[:,0]
v_LH_B = BB[:,1]
p_LH_B = BB[:,2]
#T_LH_B = BB(:,4);
#c_LH_B = BB(:,5);

#u_R_C = CC(:,1);
#v_R_C = CC(:,2);
#p_R_C = CC(:,3);
#T_R_C = CC(:,4);
#c_R_C = CC(:,5);

u_LL_D = DD[:,0];
v_LL_D = DD[:,1];
p_LL_D = DD[:,2];
#T_LL_D = DD(:,4);
#c_LL_D = DD(:,5);
u_TP_E = EE[:,4]
v_TP_E = EE[:,5]
p_TP_E = EE[:,3]
   
###########################isobar single phase 8,5,2,1######################################
u_HT_isop8 = []
v_HT_isop8 = []
u_HT_isop30 = []
v_HT_isop30 = []
u_HT_isop5 = []
v_HT_isop5 = []
u_HT_isop2 = []
v_HT_isop2 = []
u_HT_isop1 = []
v_HT_isop1 = []
#
u_LH_isop = []
v_LH_isop = []
u_LH_isop5 = []
v_LH_isop5 = []
u_LH_isop30 = []
v_LH_isop30 = []
#
u_LL_isop30 = []
v_LL_isop30 = []
#
u_TP_isop8 = []
v_TP_isop8 = []
u_TP_isop5 = []
v_TP_isop5 = []
u_TP_isop2 = []
v_TP_isop2 = []
u_TP_isop1 = []
v_TP_isop1 = []
u_TP_isop05 = []
v_TP_isop05 = []
for i in range(len(u_HT_A)):
    if p_HT_A[i] <=30.1e6 and p_HT_A[i]>=29.8e6:
        u_HT_isop30.append(u_HT_A[i]/1000)
        v_HT_isop30.append(v_HT_A[i])
    if p_HT_A[i] <=8.1e6 and p_HT_A[i]>=7.9e6:
        u_HT_isop8.append(u_HT_A[i]/1000)
        v_HT_isop8.append(v_HT_A[i])
    if p_HT_A[i] <=5.01e6 and p_HT_A[i]>=4.9e6:
        u_HT_isop5.append(u_HT_A[i]/1000)
        v_HT_isop5.append(v_HT_A[i])
    if p_HT_A[i] <=2.01e6 and p_HT_A[i]>=1.974e6:
        u_HT_isop2.append(u_HT_A[i]/1000)
        v_HT_isop2.append(v_HT_A[i])
    if p_HT_A[i] <=1.01e6 and p_HT_A[i]>=0.976e6:
        u_HT_isop1.append(u_HT_A[i]/1000)
        v_HT_isop1.append(v_HT_A[i])
#        p_HT_isop.append(p_HT_A[i])
for i in range(len(u_LH_B)):
    if p_LH_B[i] <=8.001e6 and p_LH_B[i]>=7.9e6:
        u_LH_isop.append(u_LH_B[i]/1000)
        v_LH_isop.append(v_LH_B[i])
    if p_LH_B[i] <=5.001e6 and p_LH_B[i]>=4.95e6:
        u_LH_isop5.append(u_LH_B[i]/1000)
        v_LH_isop5.append(v_LH_B[i])
    if p_LH_B[i] <=30.1e6 and p_LH_B[i]>=29.9e6:
        u_LH_isop30.append(u_LH_B[i]/1000)
        v_LH_isop30.append(v_LH_B[i])
#
for i in range(len(u_LL_D)):
    if p_LL_D[i] <=30.00001e6 and p_LL_D[i]>=29.9e6:
        u_LL_isop30.append(u_LL_D[i]/1000)
        v_LL_isop30.append(v_LL_D[i])
#############################isobar two-phase####################################
for i in range(len(u_TP_E)):
    if p_TP_E[i] <=0.52e6 and p_TP_E[i]>=0.495e6:
        u_TP_isop05.append(u_TP_E[i]/1000)
        v_TP_isop05.append(v_TP_E[i])
    if p_TP_E[i] <=5.001e6 and p_TP_E[i]>=4.99e6:
        u_TP_isop5.append(u_TP_E[i]/1000)
        v_TP_isop5.append(v_TP_E[i])
    if p_TP_E[i] <=2.001e6 and p_TP_E[i]>=1.99e6:
        u_TP_isop2.append(u_TP_E[i]/1000)
        v_TP_isop2.append(v_TP_E[i])
    if p_TP_E[i] <=1.001e6 and p_TP_E[i]>=0.999e6:
        u_TP_isop1.append(u_TP_E[i]/1000)
        v_TP_isop1.append(v_TP_E[i])


########################################################################
    
    
plt.figure(figsize=(9,4)) 

plt.semilogx(sp_R[:,0],sp_R[:,2]/1e3,color='r',linewidth=2)
plt.semilogx(sp_R[:,1],sp_R[:,2]/1e3,color='k',linewidth=2)

plt.semilogx(sp_HT[:,0],sp_HT[:,2]/1e3,linestyle='--',color='k',linewidth=2)
plt.semilogx(sp_HT[:,1],sp_HT[:,2]/1e3,linestyle='--',color='k',linewidth=2)

plt.semilogx(sp_LH[:,0],sp_LH[:,2]/1e3,color='r',linewidth=2)
plt.semilogx(sp_LH[:,1],sp_LH[:,2]/1e3,linestyle='--',color='k',linewidth=2)

plt.semilogx(sp_LL[:,0],sp_LL[:,2]/1e3,linestyle='-',color='b',linewidth=2)
plt.semilogx(sp_LL[:,1],sp_LL[:,2]/1e3,linestyle='--',color='k',linewidth=2)

plt.semilogx(v_HT_isop30,u_HT_isop30,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_HT_isop8,u_HT_isop8,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_HT_isop5,u_HT_isop5,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_HT_isop2,u_HT_isop2,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_HT_isop1,u_HT_isop1,color='k',linewidth=2,linestyle=':')
#
plt.semilogx(v_LH_isop,u_LH_isop,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_LH_isop5,u_LH_isop5,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_LH_isop30,u_LH_isop30,color='k',linewidth=2,linestyle=':')
#
plt.semilogx(v_LL_isop30,u_LL_isop30,color='k',linewidth=2,linestyle=':')
#########################TWO-PHASE########################################
plt.semilogx(v_TP_isop5,u_TP_isop5,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_TP_isop2,u_TP_isop2,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_TP_isop1,u_TP_isop1,color='k',linewidth=2,linestyle=':')
plt.semilogx(v_TP_isop05,u_TP_isop05,color='k',linewidth=2,linestyle='--')


#plt.semilogx([3.62E-002], [-279.073504], 'bo',markersize=5)


#####control des labels
ax = plt.gca() 
#ax.set_xlim([0.015, 1e-1])
ax.set_ylim([-4.5e2, 3e2])
#plt.title('Isobaric curve',fontsize=17)
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
ax.text(1.7e-3,50, '$50MPa$', fontsize=10)
ax.text(2.7e-3,50, '$30MPa$', fontsize=10)
ax.text(8.9e-3,50, '$8MPa$', fontsize=10)
ax.text(1.8e-2,50, '$5MPa$', fontsize=10)
ax.text(2.8e-2,50, '$2MPa$', fontsize=10)
ax.text(5.5e-2,50, '$1MPa$', fontsize=10)
ax.text(1.1e-1,50, '$0.5MPa$', fontsize=10)
#############################################################################
plt.tight_layout()
plt.savefig("isobar.pdf")
plt.show()
