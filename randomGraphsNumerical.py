import sys
import numpy as np
import scipy as sp
import networkx as nx
from transitionMatrix import createTransitionMatrix
from transitionMatrixDirected import createTransitionMatrixDirected
from numpy import binary_repr
import scipy.sparse as sps
from scipy.sparse import csr_matrix # compressed sparse row matrix
from scipy.sparse import linalg # for using eigs(), spsolve

    
    
def fixationProb(popSize, transitionMatrix):
    '''
    Calculates the fixation probability for a given transition matrix.


    Parameters
    ----------
    popSize: int
    transitionMatrix: sparse coo_matrix

    Returns
    -------
    probOne: float (this is the fixation probability for a single mutant starting at a random node)

    '''
    
    numStates = 2**popSize
    
    transitionMatrix = csr_matrix(transitionMatrix) #convert to compressed sparse row format
    
    fixatingColumn = transitionMatrix[0:numStates-2,numStates-2] * (-1) # this is the second to last column, it refers to the all-mutant-state
    
    # cut off the last two rows and columns (absorbing states)
    transitionMatrix = transitionMatrix[0:numStates-2,0:numStates-2]
    Id = sps.eye(numStates-2, format='csr')
    minusId = transitionMatrix - Id
    
    
    probs = linalg.spsolve(minusId, fixatingColumn) # solves the matrix equation minusId*probs=fixatingColumn
    
    
    # average over all states that start from one mutant
    probOne = sum(probs[0:popSize]) / float(popSize)
    
    
    return probOne


def isIsothermal(graph, popSize, update, directed=False):
    """
    Check whether a graph is isothermal and return the degree of the first node.
    
    
    Parameters
    ----------
    graph: networkx graph object
    popSize: int
    update: str (either 'BD' or 'DB')

    Returns
    -------
    boolean: bool (True, if graph is isothermal)
    degree: int
    
    """
    
    if directed == False:
        degDict = graph.degree()
        degList = degDict.values()
        boolean = degList.count(degList[0]) == len(degList) # is the number of elements that are equal to the first one, identical to the length?
        degree = degList[0]
        
    else:
        temperature = np.zeros(popSize)
    
        if update == 'BD':
        
            for node in graph.nodes():
                predecessors = np.array(graph.predecessors(node), dtype=float) #incoming neighbors
                predOutDegree = np.array(graph.out_degree(predecessors).values(), dtype = float) #out-degree of incoming nodes
                ones = np.ones_like(predOutDegree, dtype=float)
                inverseOutDegree = ones / predOutDegree #inverse of the out-degree of predecessors, i.e. incoming nodes.
                temperature[node] = np.sum(inverseOutDegree)
            
        
        elif update == 'DB':
        
            for node in graph.nodes():
                successors = np.array(graph.successors(node), dtype=float) #outgoing neighbors
                sucInDegree = np.array(graph.in_degree(successors).values(), dtype = float) #in-degree of outgoing neighbors
                ones = np.ones_like(sucInDegree, dtype=float)
                inverseInDegree = ones / sucInDegree #inverse of the in-degree of successors, i.e. outgoing nodes.
                temperature[node] = np.sum(inverseInDegree)        
       
        boolean = (min(temperature) == max(temperature)) #True if all elements are equal
    
        degree = 0
        if boolean == True:
            degree = graph.in_degree(1)
    
  
    return boolean, degree


def main():
    popSize = int(sys.argv[1])
    probLinkConnection = float(sys.argv[2])
    #seed = int(sys.argv[3])
    graphNumber = int(sys.argv[3])
    update = str(sys.argv[4])
    direction = str(sys.argv[5]) 
    
    if direction == 'directed':
        directed = True
    elif direction == 'undirected':
        directed = False
    else:
        print 'direction must be either directed or undirected'
    
    
    # create a random graph without self-loops:
    graph = nx.gnp_random_graph(popSize, probLinkConnection, seed=None, directed=directed)
    
    
    file = open('output/Size_'+str(popSize)+'_'+str(direction)+'_'+'prob_'+str(probLinkConnection)+'_update_'+str(update)+'_graph_'+str(graphNumber)+'.txt', 'wb')

    
    
    if directed == False:
        connected = (nx.is_connected(graph))
        
    
    else:
        connected = (nx.is_weakly_connected(graph)) #works only for directed graphs
    
        inDegree = list(graph.in_degree().values())
        count = 0
        oneRoot = False
        multiRoots = False
        for i in inDegree:
            if i == 0:
                count = count + 1
        if count == 1:
            oneRoot = True
        elif count > 1:
            multiRoots = True
        
    # Check for isothermality:    
    isIsothermalGraph, degree = isIsothermal(graph, popSize, update, directed)
    
    
  
    
    fitnessArray = np.array([0.75, 1, 1.25, 1.5, 1.75], dtype=float)
    
    

    #discard graphs that are disconnected or have multiple roots
    if (connected == False):
        prob = 0
        for fit in fitnessArray:
            file.write(str(prob)+' ')
        file.write(str(graph.edges()))
    elif  (directed == True) and (oneRoot == True):
        prob = -3
        for fit in fitnessArray:
            file.write(str(prob)+' ')
        file.write(str(graph.edges()))
    elif  (directed == True) and (multiRoots == True):
        prob = 3
        for fit in fitnessArray:
            file.write(str(prob)+' ')
        file.write(str(graph.edges()))
    
    else:
        if isIsothermalGraph == True:
            if update == 'BD' or (update == 'DB' and degree == popSize-1):
                prob = 2 
                for fit in fitnessArray:
                    file.write(str(prob)+' ')
                file.write(str(graph.edges()))
            
            if update == 'DB' and degree < popSize-1: #because for DB updating, the isothermal theorem doesn't hold
                for fitnessMutants in fitnessArray:
                    
                    if directed == False:
                        transition = createTransitionMatrix(fitnessMutants, popSize, graph, update)
                        
                    else:
                        transition = createTransitionMatrixDirected(fitnessMutants, popSize, graph, update)
    
                    prob = fixationProb(popSize, transition)
        
                    file.write(str(prob)+' ')
      
                file.write(str(graph.edges())+' regular of degree '+str(degree))
                                
        
        else:
            
            for fitnessMutants in fitnessArray:
        
                if directed == False:
                    transition = createTransitionMatrix(fitnessMutants, popSize, graph, update)
                    
                else:
                    transition = createTransitionMatrixDirected(fitnessMutants, popSize, graph, update)
    
                prob = fixationProb(popSize, transition)
        
                file.write(str(prob)+' ')
      
            file.write(str(graph.edges()))  
         
    file.close()
            
    
    
    
    
    
    pass
    
    
    
if __name__ == '__main__': 
    main()
    pass
    
    
            
            
            
            
            
            
            
            
            
            
            
            
            