import matplotlib.pyplot as plt
import numpy.random as rdm
import pickle as pk
import numpy as np

from scipy.stats import pearsonr
from scipy.stats import linregress
from network_generator import small_world


trials=100

# network parameters
N=100# number of players

# game paramters
alpha=0.5
threshold=-1-alpha # less than this and player chooses random
infection_parameter=0.05
T=100 # number of time steps
   

fig=plt.figure(figsize=(10,10))

m=0
degree_matrix=[]
sickness_matrix=[]

k_values=[i for i in range(1,51)]
r_values=[i*0.01 for i in range(50)]

#k_values=[49]
#r_values=[0]



for k in k_values:
    degree_row=[]
    sickness_row=[]
    for r in r_values:
        
#    for heterogeneity in [0]:#i*0.05 for i in range(0,90,10)]:
        print('k=',k,'r=',r)
        m=m+1
        
        sickness_correlation_list=[]
        degree_correlation_list=[]
        
        for trial in range(trials):
            #print('Trial',trial)
            
 
            nodes,edges=small_world(N,k,r)
    
            neighbours={}           
            for i in nodes:
                neighbours[i]=[]
                
            for edge in edges:
                neighbours[edge[0]].append(edge[1])
                neighbours[edge[1]].append(edge[0])
     
            # this disctionary keeeps the disease status of every player
            sickness={}
            choice={}
            # set the initial conditions
            for i in nodes:
                # choose randomly from neighbours
                choice[i]=neighbours[i][int(len(neighbours[i])*rdm.random())]
            
                # random sickness
                sickness[i]=rdm.random()
                #sickness[i]=0.5
            
            total=sum(sickness.values())
            for i in sickness:
                sickness[i]=(1/2)*N*sickness[i]/total
            
            
            disease={}        
            # track some quantities over time
            for i in nodes:
                disease[i]=[sickness[i]]
      
            for t in range(T):
                # each timestep change the choice of each player
                number_of_times_chosen={}
                for node in nodes:
                    number_of_times_chosen[node]=0
                
                for i in nodes:
                    # keep track of how many times that player was chosen
                    number_of_times_chosen[choice[i]]=number_of_times_chosen[choice[i]]+1
                    
            
                old_sickness=sickness.copy()
                for i in nodes:
                    infection=infection_parameter*old_sickness[choice[i]]
    
                    payoff=-number_of_times_chosen[choice[i]]-(1/infection_parameter)*infection
          
                    sickness[i]=sickness[i]+infection-number_of_times_chosen[i]*sickness[i]*infection_parameter
                    
                    # track some quantities over time
                    disease[i].append(sickness[i])             
                    #update the choice of the player
                    if payoff<threshold:
                        previous_choice=choice[i]
    
                        choice[i]=neighbours[i][int(len(neighbours[i])*rdm.random())]
    
            #mean_sickness={}
            degree={}
            for i in nodes:
                degree[i]=len(neighbours[i])
#                mean_sickness[i]=sum(disease[i])/T 
            
            neighbour_sickness={}
            for i in nodes:
                neighbour_sickness[i]=sum([sickness[j] for j in neighbours[i]])/len(neighbours[i])
            
            x=[sickness[i] for i in nodes]
            y=[degree[i] for i in nodes]
#            plt.figure()
#            plt.title('Degree, k='+str(k)+', r='+str(r))
#            plt.scatter(x,y)
            
            degree_correlation,p=pearsonr(x,y)
            
            degree_correlation_list.append(degree_correlation)
            y=[neighbour_sickness[i] for i in nodes]
#            plt.figure()
#            plt.title('Sickness, k='+str(k)+', r='+str(r))
#            plt.scatter(x,y)
    
            difference=np.mean([abs(sickness[i]-neighbour_sickness[i]) for i in nodes])
    
            sickness_correlation,p=pearsonr(x,y)
            slope, intercept, r_value, p_value, std_err=linregress(x,y)
            sickness_correlation_list.append(slope)
#            if p<0.01:
#                sickness_correlation_list.append(1)
#            elif p<0.1:
#                sickness_correlation_list.append(0.5)



        sickness_row.append(np.mean(sickness_correlation_list))
        degree_row.append(np.mean(degree_correlation_list))
        #print(degree_correlation,sickness_correlation)
        #print(sickness_correlation)
        
    degree_matrix.append(degree_row)
    sickness_matrix.append(sickness_row)
##################################################
pk.dump(sickness_matrix,open('../plot_data/sickness_matrix.p','wb'))

#
#Values=np.array(sickness_matrix)
#
#plt.figure()
##plt.xticks([i*10 for i in range(9)],[m_values[i*10] for i in range(9)], fontsize=fs)
##plt.yticks([9+i*10 for i in range(9)],[int(h_values[9+i*10]) for i in range(9)], fontsize=fs)
# 
#fs=15
#plt.xlabel('Rewire probability',size=fs)
#plt.ylabel('Mean degree',size=fs)
#
#plt.xticks([0,2,4,6,8],[0,0.2,0.4,0.6,0.8])
#plt.yticks([0,2,4,6,8],[1,10,20,30,40])
#
#plt.imshow(Values,cmap='Greys') # displays in color
#cbar=plt.colorbar(fraction=0.05, pad=0.04)
#cbar.outline.set_visible(False)
#cbar.set_label('Sickness correlation',size=fs)
#cbar.ax.tick_params(labelsize=fs, length=0)
# 
#plt.savefig('June_2nd/figures/Sickness_correlation.png',dpi=256,bbox_inches='tight')

