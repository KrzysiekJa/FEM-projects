import numpy as np




class Mocks(object):
    width = 2  # (x, y)
    
    def __init__(self, globalObj):
        self.height = globalObj.nNodes
        self.nodes  = np.zeros((self.height, self.width), dtype=int)
        
        x, y, k = 0, 0, 0
        delta_x = globalObj.Height/globalObj.numH
        delta_y = globalObj.Width/globalObj.numW
        
        
        for i in range(globalObj.numW):
            x = i * delta_x
            
            for j in range(globalObj.numH):
                y = j * delta_y
                
                self.nodes[k, 0] = x
                self.nodes[k, 1] = y
                k += 1



class Elements(object):
    width = 4 # number of nodes for grid's cells
    
    def __init__(self, height, nElems, numH):
        self.height = height
        self.ID     = np.zeros((height, self.width), dtype=int)
        
        j = 0; r = numH - 1
        
        for i in range(nElems):
            self.ID[i, 0] = j
            self.ID[i, 1] = self.ID[i, 0] + numH
            self.ID[i, 2] = self.ID[i, 1] + 1
            self.ID[i, 3] = self.ID[i, 0] + 1
            j += 1
            
            if not (j % r):
                 j += 1
                 r += j



class GlobalData(object):
    
    def __init__(self, filepath):
        self.Height = 0
        self.Width  = 0
        self.numH = 0
        self.numW = 0
        self.nElems = 0
        self.nNodes = 0
        
        try:
            with open(filepath, 'r') as file:
                lines = file.readlines()
            
            self.Height = float(lines[0])
            self.Width  = float(lines[1])
            self.numH   = int(lines[2])
            self.numW   = int(lines[3])
            
        except (TypeError, FileNotFoundError) as err:
            print("Occured:", err)
        
        finally:
            file.close()
        
        
        self.nElems = int((self.numH - 1)*(self.numW - 1))
        self.nNodes = int(self.numH * self.numW)




def quadratureCalc(fun, x1, x2, rank):
    period = 1./(1 -(-1))    # normalization
    det_J  = np.abs(x1 - x2) * period
    
    
    if rank == 2:
        Node1 = (1 + (1/np.sqrt(3))) * period
        Node2 = (1 - (1/np.sqrt(3))) * period
        
        pc1 = Node1*x1 + Node2*x2
        pc2 = Node1*x2 + Node2*x1
        
        return (fun(pc1) *1 + fun(pc2) *1)* det_J
    
    
    if rank == 3:
        Node1 = (1 + np.sqrt(3/5)) * period
        Node2 = period
        Node3 = (1 - np.sqrt(3/5)) * period
        
        pc1 = Node1*x1 + Node3*x2
        pc2 = Node2 * (x1 + x2)
        pc3 = Node1*x2 + Node3*x1
        
        return (fun(pc1) *(5/9) + fun(pc2) *(8/9) + fun(pc3) *(5/9))* det_J
    
    
    else:
        
        return "Cannot calculate this rank."
    