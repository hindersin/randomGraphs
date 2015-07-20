import sys
import numpy as np



def main():
    
    popSize = int(sys.argv[1])
    rows = int(sys.argv[2])
    stepSize = float(sys.argv[3])
    update = str(sys.argv[4])
    direction = str(sys.argv[5])
    
    probLinkConnectionArray = np.arange(0.0, 1.0 + stepSize, stepSize) #array from 0 to 1 in stepSize steps
    
    fitnessArray = np.array([0.75, 1, 1.25, 1.5, 1.75], dtype=float)
    cols = len(fitnessArray)
    
    loadData = np.zeros((rows, cols), dtype=float)
    
    
    for probLinkConnection in probLinkConnectionArray:
        
        for graphNumber in xrange(rows):
            loadData[graphNumber,:] = np.loadtxt('output/Size_'+str(popSize)+'_'+str(direction)+'_'+'prob_'+str(probLinkConnection)+'_update_'+str(update)+'_graph_'+str(graphNumber)+'.txt', usecols = range(cols))

        
        writerAmpl = open('output/AmplSize_'+str(popSize)+'_'+str(direction)+'_'+'prob_'+str(probLinkConnection)+'_update_'+str(update)+'.txt', 'wb')
        

        for i in xrange(rows):
        
            if loadData[i,0] == 0: #disconnected graph
                writerAmpl.write('-2 \n')
            
            elif loadData[i,0] == 2: #isothermal graph
                writerAmpl.write('2 \n')
                
            elif loadData[i,0] == 3: #multiple rooted graph
                writerAmpl.write('3 \n')
                
            elif loadData[i,0] == -3: #one-rooted graph, has prob=1/N
                writerAmpl.write('-3 \n') 
            
            else:
            
                amplifierQ = np.array([], dtype=bool)
    
        

                for r in fitnessArray:
                    fitnessMutant = r
            
                    col = np.where(fitnessArray==r)[0][0]
        
                    fixProb = loadData[i, col]

                    if fitnessMutant != 1:
                        if update == 'BD':
                            rhoMix = (1 - 1/fitnessMutant) / (1 - 1/(fitnessMutant**float(popSize)))
                        elif update == 'DB':
                            rhoMix = ((float(popSize)-1) * (fitnessMutant-1)) / (float(popSize) * fitnessMutant * (1 - fitnessMutant**(1-float(popSize))) )

                        if fitnessMutant > 1:
                            if fixProb > rhoMix:                           
                                amplifierQ = np.append(amplifierQ, True)

                            else:
                                amplifierQ = np.append(amplifierQ, False)

                        if fitnessMutant < 1:
                            if fixProb < rhoMix:
                                amplifierQ = np.append(amplifierQ, True)
                            else:
                                amplifierQ = np.append(amplifierQ, False)


                    else:
                        rhoMix = 1 / float(popSize)
                        if (abs(fixProb - rhoMix) > 0.01):
                            print 'something wrong with the neutral case'


        

                if (sum(amplifierQ) == len(amplifierQ)):
                    writerAmpl.write('1 \n')

                elif (sum(amplifierQ) == 0):
                    writerAmpl.write('-1 \n')

                else:
                    writerAmpl.write('0 \n')
                    
                    
                      
    writerAmpl.close()
    
    pass
    
    
    
if __name__ == '__main__': 
    main()
    pass
    
    