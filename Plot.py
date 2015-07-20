import numpy as np
import matplotlib.pyplot as plt
#import sys
from sympy import binomial, Symbol
from sympy.mpmath import *

def ProbConnected(n,p):
    '''
    Gives the probability that an undirected graph with n nodes is connected, if p is the probability for each link to be independently present.

    Parameters
    ----------
    n: int, number of nodes
    p: float, probability for each link to be present in the graph

    Returns
    -------
    probability for the graph to be connected 

    '''
    mp.pretty = True #convert the numbers of nsum back to real from "mpf"
    if n == 1:
        return 1
    else:
        return 1 - nsum(lambda k: binomial(n-1,k-1) * ProbConnected(k,p) * (1-p)**(k*(n-k)), [1,n-1])

def Plot(popSize, rows, stepSize, update, direction):
    
    probVec = np.arange(0.0, 1.0 + stepSize, stepSize) #array from 0 to 1 in stepSize steps
    
    Iso = np.array([], dtype=int)
    Ampl = np.array([], dtype=int)
    Neither = np.array([], dtype=int)
    Supp = np.array([], dtype=int)
    Disconnected = np.array([], dtype=int)
    MultiRoots = np.array([], dtype=int)
    OneRoot = np.array([], dtype=int)

    for probLinkConnection in probVec:
        reader = open('output/AmplSize_'+str(popSize)+'_'+str(direction)+'_'+'prob_'+str(probLinkConnection)+'_update_'+str(update)+'.txt', 'r')
        line = reader.readline()
        
        totalIso = 0
        totalAmpl = 0
        totalNeither = 0
        totalSupp = 0
        totalDisconnected = 0
        totalMultiRoots = 0
        totalOneRoot = 0
        while line != '':
            
            if line == '2 \n':
                totalIso += 1
            elif line == '1 \n':
                totalAmpl += 1
            elif line == '0 \n':
                totalNeither += 1
            elif line == '-1 \n':
                totalSupp += 1
            elif line == '-2 \n':
                totalDisconnected += 1
            elif line == '3 \n':
                totalMultiRoots += 1
            elif line == '-3 \n':
                totalOneRoot += 1
            line = reader.readline()
        reader.close()
    
        Iso = np.append(Iso, totalIso)
        Ampl = np.append(Ampl, totalAmpl)
        Neither = np.append(Neither, totalNeither)
        Supp = np.append(Supp, totalSupp)
        Disconnected = np.append(Disconnected, totalDisconnected)
        MultiRoots = np.append(MultiRoots, totalMultiRoots)
        OneRoot = np.append(OneRoot, totalOneRoot)
        
        
    width = 0.05
    
    opac = 0.5
    lines = 0.1
    fig, ax = plt.subplots(1)
    
    if direction == 'undirected':
        p0 = plt.bar(probVec, Disconnected, width, color='w', linewidth = lines)    
        p1 = plt.bar(probVec, Supp, width, color='r', bottom=Disconnected, alpha=opac, linewidth = lines)
        p2 = plt.bar(probVec, Neither, width, color='#FFEE00', bottom=Supp+Disconnected, alpha=opac, linewidth = lines)
        p3 = plt.bar(probVec, Ampl, width, color='#33AA22', bottom=Neither+Supp+Disconnected, alpha=opac, linewidth = lines)
        p4 = plt.bar(probVec, Iso, width, color='b', bottom=Ampl+Neither+Supp+Disconnected, alpha=opac, linewidth = lines)
    
    elif direction == 'directed':
        p0 = plt.bar(probVec, Disconnected, width, color='w', linewidth = lines)   
        p01 = plt.bar(probVec, MultiRoots, width, color='gray', bottom = Disconnected, alpha = 0.2, linewidth = lines) 
        p02 = plt.bar(probVec, OneRoot, width, color='orange', bottom = Disconnected+MultiRoots, alpha = opac, linewidth = lines)
        p1 = plt.bar(probVec, Supp, width, color='r', bottom = Disconnected+MultiRoots+OneRoot, alpha = opac, linewidth = lines)
        p2 = plt.bar(probVec, Neither, width, color='#FFEE00', bottom = Disconnected+MultiRoots+OneRoot+Supp, alpha = opac, linewidth = lines)
        p3 = plt.bar(probVec, Ampl, width, color='#33AA22', bottom = Disconnected+MultiRoots+OneRoot+Supp+Neither, alpha=opac, linewidth = lines)
        p4 = plt.bar(probVec, Iso, width, color='b', bottom = Disconnected+MultiRoots+OneRoot+Supp+Neither+Ampl, alpha=opac, linewidth = lines)
    
    p = np.arange(0.0, 1.001, 0.001)
    
    shift = 0.025
    linewidth = 1.3
   
    if direction == 'undirected':
        ProbConnectedGraph = ProbConnected(popSize,p)
        plt.plot(p+shift, rows*(1 - ProbConnectedGraph), 'k', linewidth = linewidth)
        if popSize == 4:
            plt.plot(p+shift, rows*(1-(p**6 + p**4 * (1-p)**2 * 3.0)), 'b', linewidth = linewidth)

        else:
            plt.plot(p+shift, rows*(1-p**(popSize*(popSize-1)/2.0)), 'b', linewidth = linewidth)
        
        if update == 'BD':    
            plt.legend((p4[0], p3[0], p2[0], p1[0], p0[0]), ('Isothermal','Amplifier', 'Unclassified', 'Suppressor','Disconnected'),loc=3)
        elif update == 'DB':    
            plt.legend((p4[0], p3[0], p2[0], p1[0], p0[0]), ('Complete','Amplifier', 'Unclassified', 'Suppressor','Disconnected'),loc=3)
    
    elif direction == 'directed':
        probRooted = 1 - (1 - (1 - p)**(popSize-1))**popSize
        plt.plot(p+shift, rows*probRooted, 'r', linewidth = linewidth)
    
        if update == 'BD':
            plt.legend((p4[0], p3[0], p2[0], p1[0], p02[0], p01[0], p0[0]), ('Isothermal','Amplifier', 'Unclassified', 'Suppressor', 'One-rooted', 'Multi-rooted','Disconnected'),loc=3)
        elif update == 'DB':
            plt.legend((p4[0], p3[0], p2[0], p1[0], p02[0], p01[0], p0[0]), ('Complete','Amplifier', 'Unclassified', 'Suppressor', 'One-rooted', 'Multi-rooted','Disconnected'),loc=3)
        

    if update == 'BD':
        plt.title('Birth-death on '+str(rows)+' random graphs of size '+str(popSize))  
    elif update == 'DB':
        plt.title('death-Birth on '+str(rows)+' random graphs of size '+str(popSize))  
        
    
    plt.ylabel('Number of graphs')
    plt.xlabel('Probability of link connection')
    
    ticks = probVec+width/2.
    ticksLess = ticks[::2]
    plt.xticks(ticksLess, np.arange(0.0,1.05,0.1))
    
    
    plt.margins(0.01)    
    
    plt.savefig('PlotSize_'+str(popSize)+'_'+str(direction)+'_'+str(update)+'.pdf',dpi=300)
    plt.show()
    
        
    
    return 0