import numpy as np





def HPChangeByTime(globalObject, SOE, tempVec):
    dt = globalObject.timeStep
    PGlobal = SOE.PGlobal
    CGlobal = SOE.CGlobal
    
    if SOE.HCalc is None:
        SOE.HCalc = SOE.HGlobal + (CGlobal / dt)
    HGlobal = SOE.HCalc
    
    PGlobal = -( -np.matmul(tempVec, (CGlobal / dt)) + PGlobal)
    
    return (HGlobal, PGlobal)





def simulate(globalObject, SOE, Nodes):
    tempVec = np.zeros((1, len(Nodes.nodes)))
    
    for i, elem in enumerate(Nodes.nodes[:, 2]):
        tempVec[0, i] = elem
        
    iterNum = int(globalObject.simulationTime / globalObject.timeStep)
    
    
    for i in range(1, iterNum+1):
        H, P = HPChangeByTime(globalObject, SOE, tempVec)
        
        tempVec = np.linalg.solve(H, P.T).T
        
        print('{:.0f}\t{:.3f}  {:.3f}'.format(globalObject.timeStep * i, np.min(tempVec), np.max(tempVec)))
