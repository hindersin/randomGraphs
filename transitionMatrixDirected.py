import numpy as np 
import networkx as nx
from numpy import binary_repr
from scipy.sparse import coo_matrix # Efficient for incremental construction of sparse matrix, but convert to csr_matrix for operations on matrix.

def createTransitionMatrixDirected(fitnessMutants, popSize, graph, update='BD'):
    '''
    Creates the transition matrix for a given graph and fitness.

    Parameters
    ----------
    fitnessMutants: float
    popSize: int
    graph: networkx graph object
    updating: BD for birth-death (default), DB for death-birth

    Returns
    -------
    transitionMatrix: sparse coo_matrix 

    '''
    

    numStates = 2**popSize
    
    states = np.zeros((numStates, popSize), dtype=int) # Represent states with binary vector.
    
    for i in xrange(0, numStates):
        rep = binary_repr(i, width=popSize)
        states[i,:] = [int(item) for items in rep for item in items.split(',')] # Convert rep from str to int and split them into a list.
    
    
    orderedStates = sorted(states, key=sum) # Ordered by increasing number of mutants, but splits them into single rows.
    
    states = np.vstack(orderedStates)   # Stacking the rows on top of each other again.
    
    # Put the all-wild-type state (0,0,0..) at the end:
    temp = states[0,:]
    states = states[1:numStates,:]
    states = np.append(states, temp)
    states = np.reshape(states, (numStates, popSize), order='C')
    
    data = np.array([], dtype=float)
    row = np.array([], dtype=int)
    col = np.array([], dtype=int)
    
    for i in xrange(0, numStates): 
        for j in xrange(0, numStates):
            state_i = states[i,:]
            state_j = states[j,:]
            
            if update == 'BD':
                
                # Check that state j can be reached within one time step.
                if (sum(abs(state_i - state_j)) == 1):
                    # Increasing the number of mutants. If sum is -1, then j has one more mutant than i.
                    if (sum(state_i - state_j) == -1):
                        changingPosition = np.where((state_i - state_j) == -1)[0][0]
                        mutants = np.where(state_i == 1)[0]
                        numMutants = len(mutants)
                        prob = 0.0
                        for k in xrange(0, numMutants):
                            if len(graph.neighbors(mutants[k])) == 0:
                                prob = prob
                            else:
                                # has_edge returns True if there is a directed links between those nodes. The function graph.neighbors(x) gives successors of node x.
                                prob = prob + float(graph.has_edge(mutants[k], changingPosition)) / float(len(graph.neighbors(mutants[k])))
                                                    
                        
                        prob = prob * fitnessMutants / (fitnessMutants * numMutants + popSize - numMutants)
                    
                    # Decreasing the number of mutants. If the sum is 1, then state j has one less mutant than state i.
                    elif (sum(state_i - state_j) == 1):
                        changingPosition = np.where((state_i - state_j) == 1)[0][0] 
                        wilds = np.where(state_i == 0)[0]
                        numMutants = popSize - len(wilds)
                        prob = 0.0
                        for k in xrange(0, len(wilds)):
                            if len(graph.neighbors(wilds[k])) == 0:
                                prob = prob
                            else:
                                prob = prob + float(graph.has_edge(wilds[k], changingPosition)) / float(len(graph.neighbors(wilds[k])))
                                            
                        prob = prob / (fitnessMutants * numMutants + len(wilds))
                        
                    if (prob != 0.0):
                        data = np.append(data, prob)
                        row = np.append(row, i)
                        col = np.append(col, j)
                    
                
                
            elif update == 'DB':
                # Check that state j can be reached within one time step.
                if (sum(abs(state_i - state_j)) == 1):
                    # If sum is -1, then state j has one more mutant than state i.
                    if (sum(state_i - state_j) == -1):
                        mutants = np.where(state_i == 1)[0]
                        changingPosition = np.where((state_i - state_j) == -1)[0][0]
                        neighborsOfChanging = graph.predecessors(changingPosition)
                        mutantNeighbors = np.intersect1d(mutants, neighborsOfChanging)
                        numMutantNeighbors = len(mutantNeighbors)
                        
                        prob = 1/ float(popSize) * (fitnessMutants * numMutantNeighbors) / (fitnessMutants * numMutantNeighbors + len(neighborsOfChanging) - numMutantNeighbors)
                
                
                    # If sum is 1, the state j has one less mutant than i.
                    elif (sum(state_i - state_j) == 1):
                        wilds = np.where(state_i == 0)[0]
                        changingPosition = np.where((state_i - state_j) == 1)[0][0]
                        neighborsOfChanging = graph.predecessors(changingPosition)
                        wildtypeNeighbors = np.intersect1d(wilds, neighborsOfChanging)
                        numWildtypeNeighbors = len(wildtypeNeighbors)
                        
                        prob = 1/ float(popSize) * (numWildtypeNeighbors) / (fitnessMutants * (len(neighborsOfChanging) - numWildtypeNeighbors) + numWildtypeNeighbors)
                
                
                
                    if (prob != 0.0):
                        data = np.append(data, prob)
                        row = np.append(row, i)
                        col = np.append(col, j)
    
    
        
    for i in xrange(0, numStates):
        sumRow = np.sum(data[row==i]) 
        data = np.append(data, 1-sumRow)
        row = np.append(row, i)
        col = np.append(col, i)   
     
        
    transitionMatrix = coo_matrix((data, (row,col)), shape = (numStates, numStates), dtype=float)
    
    
    return transitionMatrix
    
    
    
    
    
    
    
    
    
    
    