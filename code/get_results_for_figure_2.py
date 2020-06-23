import matplotlib.pyplot as plt
import numpy.random as rdm
import networkx as nx
from scipy.optimize import fsolve
import numpy as np
from scipy.special import erf
from scipy.special import erfinv

import pickle as pk

trials=1
N=100# number of players
T=100


alpha_values=[1-0.01*i for i in range(101)]
beta_values=[0.002*i for i in range(101)]


M=[]

for alpha in alpha_values:
    print(alpha)
    m=[]
    for beta in beta_values:

        split_list=[]
        for trial in range(trials):
        
            # this dictionary keeeps the disease status of every player
            sickness={}
            choice={}
            
            # set the initial conditions
            for i in range(N):
                # choose random players
                choice[i]=i
                while choice[i]==i:
                    choice[i]=int(N*rdm.random())
                
                # homogeneous sickness
                sickness[i]=1/2
            
            # split is true if at least one is on the pther side of the threshold
            split=False
            # count the rounds
            t=0
            
            while split==False and t<T:
                t=t+1
                # each timestep change the choice of each player
                number_of_times_chosen=[0 for i in range(N)]
                
                for i in range(N):
                    # keep track of how many times that player was chosen
                    number_of_times_chosen[choice[i]]=number_of_times_chosen[choice[i]]+1
                    
                old_sickness=sickness.copy()
                for i in range(N):
                    infection=beta*old_sickness[choice[i]]
                    
                   
                    sickness[i]=sickness[i]+infection-number_of_times_chosen[i]*sickness[i]*beta
                    if alpha>0.5 and sickness[i]>alpha: # and old_sickness[i]>alpha:
                        split=True
                        
                    elif alpha<=0.5 and sickness[i]<alpha: # and old_sickness[i]<alpha:
                        split=True


                    #update the choice of the player
                
                    if number_of_times_chosen[choice[i]]+old_sickness[choice[i]]>1+alpha:
                        previous_choice=choice[i]
                        #choice[i]=i
                        while choice[i]==i or choice[i]==previous_choice:
                            choice[i]=int(N*rdm.random())
            if t==T:
                split_list.append(1)
            else:
                split_list.append(0)
        m.append(np.mean(split_list))
        #print(alpha,beta,T)

    M.append(m)
    
# pickle it!!

pk.dump(M,open('../plot_data/M_'+str(N)+'_'+str(T)+'.p','wb'))

#Values=np.array(M)
#
#
#
#fig=plt.figure()
###plt.xticks([i*10 for i in range(9)],[m_values[i*10] for i in range(9)], fontsize=fs)
##plt.yticks([9+i*10 for i in range(9)],[int(h_values[9+i*10]) for i in range(9)], fontsize=fs)
# 
#fs=15
#plt.xlabel('$\\beta$',size=fs)
#plt.ylabel('$\\alpha$',size=fs,labelpad=10,rotation='horizontal')
#
#plt.xticks([0,0.1,0.2],size=fs)
#plt.yticks([0.2,0.4,0.6,0.8],size=fs)
##
#im=plt.imshow(Values,aspect='auto',extent=(0,0.2,0,1),cmap='Blues') # displays in color ,vmin=-1, vmax=1
##cbar=plt.colorbar(fraction=0.05, pad=0.04)
##cbar.outline.set_visible(False)
##cbar.set_label('$\\log_{10}($ time in initial state $)$',size=fs)
##cbar.ax.tick_params(labelsize=fs, length=0)
#
#
##### THE RED LINE #################
#
#k=(17)**(1/2)
##z=erfinv(2*((1-(1/T))**(1/N))-1)
##for z in [2,3]:
#
##z=erfinv(1-(2/(N*T)))
#
#for p_value in [1/2]:
# 
#    z=erfinv(2*p_value**(1/(N*T))-1)
#
#    print(z)
#    
#    beta=[]
#    for alpha in alpha_values:
#        if alpha>1/2:
#            a=4
#            b=5*k-17
#            c=10*k-38
#        else:
#            a=1
#            b=1
#            c=1
#        
#        if alpha==1/2:
#            beta.append(0)
#        else:
#            beta.append((a/(b+c*(z/(1-2*alpha))**2)))
#        
#    plt.plot(beta,alpha_values,c='r',linewidth=3,alpha=0.5)
#
#################
#
#plt.savefig('June_2nd/figures/phase_space.png',dpi=256, bbox_inches='tight')



