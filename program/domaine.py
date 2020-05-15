import numpy as np
import matplotlib.pyplot as plt


sp_LL=np.loadtxt('spline_LL.txt')
sp_LH=np.loadtxt('spline_LH.txt')
sp_HT=np.loadtxt('spline_HT.txt')
sp_R=np.loadtxt('spline_R.txt')
###%%% A corner %%%
A = np.zeros((10,2))
A[:,0] = np.linspace(1.23287e-3,2.138e-3,num=10)
A[:,1] = -190.311e3
###%%% B corner %%%
B = np.zeros((10,2))
B[:,0] = np.linspace(1.50208e-3,20.0e-3,num=10)
B[:,1] = -107.978e3
###%%% C corner %%%
C = np.zeros((10,2))
C[:,0] = np.linspace(20.0e-3,76.597e-3,num=10)
C[:,1] = -107.978e3
###%%% Top %%%
Top = np.zeros((10,2))
Top[:,0] = np.linspace(0.00430863,0.354716,num=10);
Top[:,1] = 540000;
   
########################################################################
    
    
plt.figure(figsize=(9,4)) 

#plt.semilogx(sp_R[:,0],sp_R[:,2]/1e3,color='orange',linewidth=2)
#plt.semilogx(sp_R[:,1],sp_R[:,2]/1e3,color='orange',linewidth=2)

#plt.semilogx(sp_HT[:,0],sp_HT[:,2]/1e3,linestyle='--',color='r',linewidth=2)
#plt.semilogx(sp_HT[:,1],sp_HT[:,2]/1e3,linestyle='--',color='r',linewidth=2)

#plt.semilogx(sp_LH[:,0],sp_LH[:,2]/1e3,color='m',linewidth=2)
#plt.semilogx(sp_LH[:,1],sp_LH[:,2]/1e3,linestyle='--',color='m',linewidth=2)

plt.semilogx(sp_LL[:,0],sp_LL[:,2]/1e3,linestyle='-',color='b',linewidth=2)
plt.semilogx(sp_LL[:,1],sp_LL[:,2]/1e3,linestyle='--',color='b',linewidth=2)
plt.semilogx(sp_LL[:,3],sp_LL[:,2]/1e3,linestyle='--',color='g',linewidth=2@)
###%%%%%%% corners %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#plt.semilogx(A[:,0],A[:,1]/1e3,linestyle='-',color='b',linewidth=2)
#plt.semilogx(B[:,0],B[:,1]/1e3,linestyle='-',color='m',linewidth=2)
#plt.semilogx(C[:,0],C[:,1]/1e3,linestyle='-',color='orange',linewidth=2)
#plt.semilogx(Top[:,0],Top[:,1]/1e3,linestyle='-',color='r',linewidth=2)



##############################################################################
#plt.plot(t, z,linewidth = 2,dashes = [10, 5, 25, 5],marker='o',markersize=8,color = 'blue')
#linestyles = ['-', '--', '-.', ':']
#dashes = [10, 5, 100, 5] # 10 points on, 5 off, 100 on, 5 off  

#####control des labels
ax = plt.gca() 
#ax.set_xlim([0.015, 1e-1])
#ax.set_ylim([-4.5e2, 6e2])
plt.title('Physical domain',fontsize=17)
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
ax.plot([2.138e-3], [-190.311], 'ro')
ax.plot([1.8626809245968707e-003], [-207.34331603119953], 'ro')
ax.annotate('critical point', fontweight='bold',xy=(2.338e-3,-196.311), xytext=(3.038e-3,-250.311),
            arrowprops=dict(facecolor='black', width=0.01,headwidth=5))
ax.annotate('saturation line',fontweight='bold', xy=(4.338e-3,-148.311), xytext=(8.338e-3,-159.311),
            arrowprops=dict(facecolor='black', width=0.01,headwidth=5))
ax.text(1.238e-3,-230.311 , 'LL', fontsize=11,fontweight='bold')
ax.text(1.838e-3,-154.311 , 'LH', fontsize=11,fontweight='bold')
ax.text(1.538e-2,-50.311 , 'HT', fontsize=11,fontweight='bold')
ax.text(6.538e-2,-125.311 , 'R', fontsize=11,fontweight='bold')
ax.text(1.538e-2,-300.311 , 'TP', fontsize=11,fontweight='bold')
ax.text(2.38e-3,30, '$50MPa$', fontsize=10)
ax.text(0.9e-1,20, '$0.5MPa$', fontsize=10)
#############################################################################
plt.tight_layout()
#plt.savefig("domain.pdf")
plt.show()
