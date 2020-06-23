import matplotlib.pyplot as plt
import numpy.random as rdm
import networkx as nx
import numpy as np
import pickle as pk

from network_generator import small_world

plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True


# network parameters
N=10# number of players
mean_degree=3

# game paramters
alpha=0.5
threshold=-1-alpha # less than this and player chooses random
beta=0.01
T=10000 # number of time steps
   

#fig=plt.figure(figsize=(10,10))
fig = plt.figure(figsize=(12,5.5))
fs=15

ax = fig.add_subplot(1,2,1)

m=0
for k in [1,2,3]:
    for p in [0,0.2,0.5]:
        print('k=',k,'p=',p)
        m=m+1
        
        nodes,edges=small_world(N,k,p)


        #print(edges)

        neighbours={}
        for i in range(N):
            neighbours[i]=[]
            
        for edge in edges:
            neighbours[edge[0]].append(edge[1])
            neighbours[edge[1]].append(edge[0])
        
        degree=[len(neighbours[i]) for i in range(N)]

#        dist=[0 for i in range(max(degree)+1)]
#        for d in degree:
#            dist[d]=dist[d]+1
#            
#        
     
        # this disctionary keeeps the disease status of every player
        sickness={}
        choice={}
        # set the initial conditions
        for i in range(N):
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
        for i in range(N):
            disease[i]=[sickness[i]]
  
        for t in range(T):
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

                    choice[i]=neighbours[i][int(len(neighbours[i])*rdm.random())]

        ############### Small network ###################################       
        #ax1 = fig.add_subplot(3,3,m)
        
        #
       
        G=nx.Graph()
        for edge in edges:
            G.add_edge(edge[0],edge[1])
        
        pos={}
        for i in nodes:           
            x=2.8*((m-1)%3)+np.sin(i*2*np.pi/N)
            y=2.8*(2-int((m-1)/3))+np.cos(i*2*np.pi/N)
            pos[i]=np.array([x,y])
        
        nodes=nx.draw_networkx_nodes(G,pos,nodelist=[node for node in disease],node_color=[sum(disease[node])/T for node in disease],node_size=100,linewidth=1,cmap=plt.cm.Greys)
        nodes.set_edgecolor('k')
        

        nx.draw_networkx_edges(G,pos,linewidth=2)
        #plt.axis('off')   

fs=15

plt.xticks([0,3,6],[0,0.2,0.5],size=fs)
plt.yticks([0,3,6],[6,4,2],size=fs)


ax.set_title('Randomized edges',size=fs,pad=10)
plt.ylabel('Mean degree',size=fs)
   
cbar=plt.colorbar(nodes,fraction=0.06, pad=0.01,ticks=[0.5,0.7,0.9],cmap='Greys',orientation='horizontal')
cbar.outline.set_visible(False)
cbar.set_label('Sickness',size=fs)

cbar.ax.tick_params(labelsize=fs, length=0)
cbar.set_clim(0.4, 1)
##################################################

sickness_matrix=pk.load(open('../plot_data/sickness_matrix.p','rb'))


Values=np.array(sickness_matrix)

ax = fig.add_subplot(1,2,2)

s=15
ax.set_title('Randomized edges',size=fs,pad=10)
plt.ylabel('Mean degree',size=fs)

plt.xticks([0,10,20,30,40],[0,0.2,0.4,0.6,0.8],size=fs)
plt.yticks([0,9,19,29,39],[2,20,40,60,80],size=fs)

plt.imshow(Values,cmap='Greys') # displays in color
cbar=plt.colorbar(fraction=0.045, pad=0.01,ticks=[0.1, 0.3,0.5],orientation='horizontal')
cbar.outline.set_visible(False)
cbar.set_label('Neighbour effect',size=fs)
cbar.ax.tick_params(labelsize=fs, length=0)

fig.text(0.51,0.9,'B',size=2*fs)
fig.text(0.08,0.9,'A',size=2*fs)

plt.subplots_adjust(wspace=0.2)
plt.savefig('../figures/network_figure.pdf',format='pdf',bbox_inches='tight',dpi=256)





