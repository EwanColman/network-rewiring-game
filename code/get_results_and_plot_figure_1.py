import matplotlib.pyplot as plt
from matplotlib import gridspec

import numpy.random as rdm
import networkx as nx
from scipy.optimize import fsolve
import numpy as np


alpha=0.5
beta=0.005

N=100 # number of players

T=5000 # number of time steps

threshold=-1-alpha # less than this and player chooses random


#recovery_parameter=1-infection_parameter
# this disctionary keeeps the disease status of every player
sickness={}

# this sstores the choice that the player makes for the current round
# at forst it is random
choice={}
# set the initial conditions
for i in range(N):
    # choose random players
    choice[i]=i
    while choice[i]==i:
        choice[i]=int(N*rdm.random())
    
    # random sickness
    sickness[i]=rdm.random()
    #sickness[i]=0.5

# normalise to get a mean of 1/2
total=sum(sickness.values())
for i in sickness:
    sickness[i]=(1/2)#*N*sickness[i]/total


# this dictionary tracks it over time
disease={}

# add the initial value
for i in range(N):
    disease[i]=[sickness[i]]

# keep track of these for plotting
top_values=[]
bottom_values=[]

# start the simulation
for t in range(T):
    
    # at this point (t=1000) count how many are in each group
    if t==1000:
       # here compute the individuals above or below the threshold
        top_group=[i for i in sickness if sickness[i]>alpha]
        bottom_group=[i for i in sickness if sickness[i]<=alpha]

        n0=len(top_group)
        # now calculate r the precise way
        sol=fsolve(lambda x: np.exp(x)-(((x**2)-((1+(n0/N))*x)+1)/(1-x)),0.99)
        r=sol[0]       
        
            
    # each timestep change the choice of each player
    number_of_times_chosen=[0 for i in range(N)]
    
    for i in range(N):
        # keep track of how many times that player was chosen
        number_of_times_chosen[choice[i]]=number_of_times_chosen[choice[i]]+1

    old_sickness=sickness.copy()
    for i in range(N):
        infection=beta*old_sickness[choice[i]]
        payoff=-number_of_times_chosen[choice[i]]-(1/beta)*infection
  
        sickness[i]=sickness[i]+infection-number_of_times_chosen[i]*sickness[i]*beta
        
        # track some quantities over time
        disease[i].append(sickness[i])
             
        #update the choice of the player
        if payoff<threshold:
            previous_choice=choice[i]
            #choice[i]=i
            while choice[i]==i or choice[i]==previous_choice:
                choice[i]=int(N*rdm.random())

    # at these times we record some outputs for plotting   
    if t % 100==0 and t>=1000:
        top_values=top_values+[sickness[i] for i in top_group]
        bottom_values=bottom_values+[sickness[i] for i in bottom_group]
     
    
#############
p=n0/N


mu0=(1/2)*(1-p*r)/(r-2*r*p+p)
mu1=(1/2)*(r*(1-p))/(r-2*r*p+p)

# low-sickness group
k_m1=r
k_m2=r*(1+r)

m2=((r*p*(mu0**2)+(1-r*p)*(mu1**2))*beta+2*(r*p*mu0+(1-r*p)*mu1)*(1-k_m1*beta)*mu0)/(2*k_m1-beta*k_m2)
var0=m2-mu0**2


# high-sickness group
k_m1=r+((1-r)/(1-p))
k_m2=(1+(2-p)*r-(1+p)*(r**2))/(1-p)

m2=((r*p*(mu0**2)+(1-r*p)*(mu1**2))*beta+2*(r*p*mu0+(1-r*p)*mu1)*(1-k_m1*beta)*mu1)/(2*k_m1-beta*k_m2)
var1=m2-mu1**2

#This is jsut a visual thing if you ant to make the curves smoother
resolution=1

n=10


plt.figure(figsize=(10,5))
fs=12
gs = gridspec.GridSpec(1, 2, width_ratios=[9, 3]) 

ax0 = plt.subplot(gs[0])
plt.xlabel('Time',size=fs)
plt.ylabel('Sickness',size=fs)
plt.xticks([0,1000,2000,3000,4000,5000],size=fs)
plt.yticks([0,0.2,0.4,0.6,0.8,1],size=fs)

for i in range(n):
    if i==0:
        lab='Individual player'
    else:
        lab=None
    plt.plot([disease[i][int(resolution*t)] for t in range(int(T/resolution))],c='k',alpha=0.4,label=lab)

# the means of the two groups
plt.plot([0,T],[mu0,mu0],':',c='k',alpha=1,label='$\mu_{0}$ and $\mu_{1}$')
plt.plot([0,T],[mu1,mu1],':',c='k',alpha=1)

# add the threshold
plt.plot([0,T],[alpha,alpha],'-.',c='k',alpha=0.6,label='$\\alpha$')
plt.ylim([0.2,0.8])

plt.legend(loc=4,prop={'size':fs})

################
ax1 = plt.subplot(gs[1])

# first curve
x_range=[i/100 for i in range(int(100*alpha),80)]
normal=[(1/100)*n0*np.exp(-((i-mu0)**2)/(2*var0))/np.sqrt(2*np.pi*var0) for i in x_range]
plt.plot(normal,x_range,c='r',linestyle='--',linewidth=2)

#second curve
x_range=[i/100 for i in range(20,int(100*alpha))]
normal=[(1/100)*(N-n0)*np.exp(-((i-mu1)**2)/(2*var1))/np.sqrt(2*np.pi*var1) for i in x_range]
plt.plot(normal,x_range,c='r',linestyle='--',linewidth=2,label='Analytical')


# add 0.5 to get the bar to appear in the right place i.e half way between i and i+1
x_range=[(i+0.5)/100 for i in range(100)]

#first histogram
dist=[0 for i in range(100)]
for value in top_values:
    dist[int(100*value)]=dist[int(100*value)]+1/40
plt.barh(x_range,dist,height=0.008,color='k',label='Simulated')

#second histogram
dist=[0 for i in range(100)]
for value in bottom_values:
    dist[int(100*value)]=dist[int(100*value)]+1/40
plt.barh(x_range,dist,height=0.008,color='k')

plt.ylim([0.2,0.8])
#plt.xticks([i/10 for i in range(2,8)]
plt.yticks([])
plt.xticks([5,10],size=fs)
plt.xlabel('Frequency',size=fs)

plt.legend(loc=4,prop={'size':fs})

plt.subplots_adjust(wspace=0)
#
#plt.xlabel('Sickness,$x$')
#plt.ylabel('Frequency')


plt.savefig('../figures/Time_series_and_distribution.pdf',format='pdf',dpi=256,bbox_inches='tight')

top_group=[i for i in sickness if sickness[i]>alpha]
print(n0,len(top_group))
