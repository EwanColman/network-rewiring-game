import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf
from scipy.special import erfinv

import pickle as pk

plt.rcParams['xtick.bottom'] =plt.rcParams['xtick.labelbottom'] =True
plt.rcParams['xtick.top'] =plt.rcParams['xtick.labeltop'] = False

fs=15
fig = plt.figure(figsize=(12,5))
gs = fig.add_gridspec(3, 7)


########## LEFT SIDE OF PLOT ################

l_range=[i/100 for i in range(100)]

beta_values=[0.001,0.05,0.1]
alpha_values=[0.6,0.5,0.4]
for i in range(3):
    for j in range(3):
        a=alpha_values[i]
        b=beta_values[j]
        ax = fig.add_subplot(gs[i, j])
        #ax.set_title('$\\alpha='+str(a)+'$, $\\beta='+str(b)+'$')

        if i==0:
            ax.set_title('$\\beta='+str(b)+'$',size=fs)
#        
        up_rate=[]
        down_rate=[]
        transition_rate=[]
        energy=[]
        for l in l_range:
            
            r=(-3+((9+8*(1+l))**(1/2)))/2
            
            mu0=(1/2)*(1-l*r)/(r-2*r*l+l)
            mu1=(1/2)*(r*(1-l))/(r-2*r*l+l)
            
            var0=(((r*l*(mu0**2)+(1-r*l)*(mu1**2))*b+2*(r*l*mu0+(1-r*l)*mu1)*(1-r*b)*mu0)/(2*r-b*r*(1+r)))-mu0**2
            var1=(((r*l*(mu0**2)+(1-r*l)*(mu1**2))*b+2*(r*l*mu0+(1-r*l)*mu1)*(1-(r+((1-r)/(1-l)))*b)*mu1)/(2*(r+((1-r)/(1-l)))-b*((1+(2-l)*r-(1+l)*(r**2))/(1-l))))-mu1**2
            
            up=(1-l)*(1/2)*(1-erf((a-mu1)/((2*var1)**(1/2))))
            down=l*(1/2)*(1+erf((a-mu0)/((2*var0)**(1/2))))
            
            up_rate.append(up)
            down_rate.append(down)
        
            transition_rate.append(up-down)
            energy.append(sum(transition_rate))
        energy=[-c for c in energy]
#            
        #Energy
        plt.plot(l_range[0:100],energy[0:100],'k',linewidth=2)
        
        # transition rate
        #plt.plot(l_range[0:int(100*l_max)],energy[0:int(100*l_max)],'b',linewidth=2)
        
        #plt.plot(l_range[int(100*l_max):100],energy[int(100*l_max):100],'b:',linewidth=2)
        
        # plt.plot([0,1],[0,0],'k:',linewidth=2)
        
        plt.xticks([])
        plt.yticks([])
        plt.xlim([-0.05,1.05])
        plt.ylim([min(energy)-5,min(energy)+40])
        
#        plt.plot([l_min,l_min],[-1,1],'k:')
#        plt.plot([l_max,l_max],[-1,1],'k:')
#        
        
        if i==2:
            print(i,j,a,b)
            plt.xticks([0,1],size=fs)
            
        if j==2:   
            ax1 = ax.twinx()
            plt.yticks([])
            plt.ylabel('$\\alpha='+str('%.1f' % a)+'$',size=fs)
        
                
       
        
fig.text(0.13,0.01,'Proportion in high sickness group, $\lambda$',size=fs)
fig.text(0.09,0.51,'Energy',size=fs,rotation='vertical')

fig.text(0.52,0.88,'B',size=2*fs)
fig.text(0.08,0.88,'A',size=2*fs)


############### RIGHT SIDE OF PLOT ##################
#
trials=1
N=100# number of players
T=100

# right hand side

alpha_values=[1-0.01*i for i in range(101)]
beta_values=[0.002*i for i in range(101)]
#
M=pk.load(open('../plot_data/M_'+str(N)+'_'+str(T)+'.p','rb'))
#
Values=np.array(M)

#
#
#
        
ax2 = fig.add_subplot(gs[0:,4:])
#ax2.set_title('X')
#plt.plot([0,1],[0,1])

plt.xlabel('$\\beta$',size=fs)
plt.ylabel('$\\alpha$',size=fs,labelpad=10,rotation='horizontal')

plt.xticks([0,0.1,0.2],size=fs)
plt.yticks([0.2,0.4,0.6,0.8],size=fs)
#
im=plt.imshow(Values,aspect='auto',extent=(0,0.2,0,1),vmax=2,cmap='Greys') # displays in color ,vmin=-1, vmax=1
#cbar=plt.colorbar(fraction=0.05, pad=0.04)
#cbar.outline.set_visible(False)
#cbar.set_label('$\\log_{10}($ time in initial state $)$',size=fs)
#cbar.ax.tick_params(labelsize=fs, length=0)


#### THE RED LINE #################

alpha_values=[1-0.01*i for i in range(101)]
beta_values=[0.002*i for i in range(101)]

T=100
N=100

k=(17)**(1/2)
#z=erfinv(2*((1-(1/T))**(1/N))-1)
#for z in [2,3]:

#z=erfinv(1-(2/(N*T)))

for p_value in [1/2]:
 
    z=erfinv(2*p_value**(1/(N*T))-1)

    print(z)
    
    beta=[]
    for alpha in alpha_values:
        if alpha>1/2:
            a=4
            b=5*k-17
            c=10*k-38
        else:
            a=1
            b=1
            c=1
        
        if alpha==1/2:
            beta.append(0)
        else:
            beta.append((a/(b+c*(z/(1-2*alpha))**2)))
        
    plt.plot(beta,alpha_values,c='r',linewidth=4,alpha=0.5)

################

plt.savefig('../figures/phase_space.pdf',format='pdf',dpi=256, bbox_inches='tight')



