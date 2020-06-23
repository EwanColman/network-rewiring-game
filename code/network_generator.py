import numpy as np
import numpy.random as rdm


def heterogeneous_dist(ID_list,mean_degree,heterogeneity):
    degree={}
    for ID in ID_list:
        degree[ID]=mean_degree
        
    while np.std([degree[ID] for ID in ID_list])<heterogeneity:
    #for i in range(switches):
        #choose a node at random
        source_list=[ID for ID in ID_list if degree[ID]>1]
        random_node=source_list[int(rdm.random()*len(source_list))]
                
        # select a node with probability proportional to its degree
        
        r=rdm.random()*len(ID_list)*mean_degree

        k=0
        target_node=ID_list[k]
        slide=degree[target_node]
       
        while slide<r:
            k=k+1
            target_node=ID_list[k]
            slide=slide+degree[target_node]
            
        degree[random_node]=degree[random_node]-1   
        degree[target_node]=degree[target_node]+1
        
    return degree


#ID_list=[i for i in range(20)]
#mean_degree=5
#switches=10
#heterogeneity=2
#
#degree=heterogeneous_dist(ID_list,mean_degree,heterogeneity)
#print(degree)


def power_law(alpha,ID_list,mean_degree):
    
    number_of_stubs=mean_degree*len(ID_list)
   
    fitness=[]    
    for ID in ID_list:
        x=((1-rdm.random())**(1/(1-alpha)))
        x=1
        fitness.append([ID,x])
        print(x)
        
    fitness=sorted(fitness, key=lambda item: item[1])
    total_fitness=sum([f[1] for f in fitness])
    degree={}
    for ID in ID_list:
        degree[ID]=0
        
    for i in range(number_of_stubs):
        
        r=rdm.random()*total_fitness
        j=0
        slide=fitness[j][1]
        while r>slide:
            j=j+1
            slide=slide+fitness[j][1]
            
            
        degree[fitness[j][0]]=degree[fitness[j][0]]+1
#
    print(degree)
    print()
    return degree

def small_world(number_of_nodes,k,p):
    
    ID_list=[i for i in range(number_of_nodes)]

    edge_list=[]
    
    # this to help avoid adding illegal or duplicate edges
    neighbours={}
    #
    for ID in ID_list:
        neighbours[ID]=[ID]
    
    for i in range(1,k+1):
        for ID in ID_list:
            # if degree ID = N-1 then ther is no space 
            if len(neighbours[ID])<number_of_nodes:
            
                r=rdm.random()
                # if its chosen to be random or if a random edge already exists in its place
                if r<p:
                    target=int(rdm.random()*number_of_nodes)
                else:
                    target=(ID+i)%number_of_nodes
                #choose randomly until new one is found
                while target in neighbours[ID]:
                    target=int(rdm.random()*number_of_nodes)
                    
                neighbours[ID].append(target)
                neighbours[target].append(ID)
                
                edge=list(sorted([ID,target]))
                edge_list.append(edge)

    return ID_list,edge_list



def ring_config_model(number_of_nodes,radius,alpha,mean_degree):
    # create list of IDs named 1 to N
    ID_list=[i for i in range(number_of_nodes)]
    
    # for the degee distribution
    #xmin=mean_degree*(alpha-2)/(alpha-1)
    
    #print('xmin',xmin)
    # get the degree for each ID
    degree=heterogeneous_dist(ID_list,mean_degree,alpha)
#    degree={}
#    for ID in ID_list:
#        #power law
#        #degree[ID]=int(xmin*((1-rdm.random())**(1/(1-alpha))))
#        #poisson
#        #degree[ID]=rdm.poisson(mean_degree)
#        # negative binomial
#        degree=rdm.negative_binomial(1-(mean/variance),(mean**2)/(variance-mean))

    #print('mean_degree=',sum([degree[ID]/number_of_nodes for ID in ID_list]))

    edge_list=[]

    #first create a list of stubs 
    for ID in ID_list:
        #print('ID=',ID)
        # create list of potential stubs to connect to
        stubs=[]
        for i in range(1,radius):
            #stubs to the right and to the left 
            for r in [-i,i]:
                #print('r=',r)
                #print('ID+r=',ID+r)
                stub=(ID+r) % number_of_nodes
                #print('stub',stub)
                # add the stubs to the list
                stubs=stubs+[stub for j in range(degree[stub])]
        
        stubs=list(rdm.permutation(stubs))
        
        # choose the stubs that ID will connect to
        chosen_stubs=[]

        
        
        if degree[ID]>len(set(stubs)):
            
            #print('Warning: degree[ID]=',degree[ID],'Distinct stubs=',len(set(stubs)))
            degree[ID]=len(set(stubs))
        # if there are not enough stubs for the degree then 
        
        
        while len(chosen_stubs)<degree[ID]:
            #print('len(chosen_stubs)',len(chosen_stubs))
            stub=stubs.pop()
            if stub not in chosen_stubs:
                #add it to the list if its not already there
                chosen_stubs.append(stub)
                #add the edge
                edge_list.append([ID,stub])
                #print([ID,stub],'added')
                # reduce degree of the one connected to
                degree[stub]=degree[stub]-1
        # put the degree of ID to 0 since all its stubs have been used
        degree[ID]=0
        
    ID_list=list(set([edge[0] for edge in edge_list]+[edge[1] for edge in edge_list]))

       
    return ID_list,edge_list

def connect_stubs(stubs,edge_list):
    
    #randomize the list
    stubs=[(s[0],s[1]) for s in rdm.permutation(stubs)]
   
    #new_edges=[]
    if len(stubs)<2:
        no_more_edges=True
    else:
        no_more_edges=False

    while no_more_edges==False:
#        print()
        #print('Stubs:',stubs)
        #take the first stub from the list
        source=stubs.pop(0)
       
        #print('Source:',source)
        found=False
        stubs_to_check=len(stubs)
        
        while found==False and no_more_edges==False:
            #take the first form the list
            target=stubs.pop(0)
            #print('Target:',target)
            #check that it is ok
            # order the edge so that the reverse one dosen't get put in
            new_edge=sorted([source,target])
            if new_edge in edge_list or source==target:
            #if not then throw it back in the list
                stubs.append(target)
            #    print('try another')
                
                
            # otherwise create the edge
            else:
#                print('New edge:',[stub,target])
                edge_list.append(new_edge)
                found=True
                
            # one less edge to check so deduct from n
            stubs_to_check=stubs_to_check-1
            if stubs_to_check<2:
                no_more_edges=True

    return edge_list

def modular_config_model(module_size,number_of_modules,p,alpha,mean_degree):
    ID_list=[]
    # create list of IDs according to the (module,number) naming convention
    for m in range(number_of_modules):
        for n in range(module_size):
            node_ID=(m,n)
            ID_list.append(node_ID)
    
    degree=heterogeneous_dist(ID_list,mean_degree,alpha) 
    
    # for the degee distribution
    edge_list=[]

    #these are the stubs that connect together across modules  
    inter_stubs=[]

    for m in range(number_of_modules):
        #first create a list of stubs within the module 
        intra_stubs=[]
        for n in range(module_size):
            # power law
            #degree=int(xmin*((1-rdm.random())**(1/(1-alpha))))
            # poisson
            #degree=rdm.poisson(mean_degree)
            # negative binomial
#            degree=rdm.negative_binomial((mean**2)/(variance-mean),1-(mean/variance))
#            print(mean,variance,(mean**2)/(variance-mean),1-(mean/variance),degree)
            
            for i in range(degree[(m,n)]):
                r=rdm.random()
                if r<p:
                    intra_stubs.append((m,n))
                else:
                    inter_stubs.append((m,n))

        # for the intra stubs 
        edge_list=connect_stubs(intra_stubs,edge_list)
                    
        #edge_list=edge_list+new_edges

    #now do the same for the inter stubs
    edge_list=connect_stubs(inter_stubs,edge_list)
        
    
#    total_degree=len(stubs)
#    new_edges=[[tuple(stubs[i]),tuple(stubs[total_degree-1-i])] for i in range(int(total_degree/2))]

    #edge_list=edge_list+new_edges
       
    return ID_list,edge_list

def module_connector_network(module_size,number_of_modules,number_of_connectors,alpha,beta):
    ### PART 1 - create the network #####################
    # create a network of disconnected communities and a few between community connectors
  
    ID_list=[]
    # create list of IDs according to the (module,number) naming convention
    for m in range(number_of_modules):
        for n in range(module_size):
            node_ID=(m,n)
            ID_list.append(node_ID)
    
    
    # start with an empty list and add edges to it
    edge_list=[]
    
    ID_list2=ID_list.copy()
    for ID1 in ID_list:
        ID_list2.remove(ID1)
        for ID2 in ID_list2:
            r=rdm.random()
            if ID1[0]==ID2[0]:
                if r<alpha:
                    edge_list.append([ID1,ID2])
    
    
    for i in range(number_of_connectors):
        node_ID=(m+1,i)
        ID_list.append(node_ID)
        #connect to each community with probability 
        for j in range(number_of_modules):
            r=np.random.random()
            if r<beta:
                edge_list.append([node_ID,(j,i)])
    
    return ID_list,edge_list

def basic_modular_network(module_size,number_of_modules,alpha,beta):
### PART 1 - create the network #####################

    ID_list=[]
    # create list of IDs according to the (module,number) naming convention
    for m in range(number_of_modules):
        for n in range(module_size):
            node_ID=(m,n)
            ID_list.append(node_ID)
    
    # start with an empty list and add edges to it
    edge_list=[]
    
    ID_list2=ID_list.copy()
    for ID1 in ID_list:
        ID_list2.remove(ID1)
        for ID2 in ID_list2:
            r=np.random.random()
            if ID1[0]==ID2[0]:
                if r<alpha:
                    edge_list.append([ID1,ID2])
            else:
                 if r<beta:
                    edge_list.append([ID1,ID2])
    return ID_list,edge_list