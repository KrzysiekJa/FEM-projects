import numpy as np




class Mocks(object):
    width = 2  # (x, y)
    
    def __init__(self, globalObj):
        self.height = globalObj.nNodes
        self.nodes  = np.zeros((self.height, self.width), dtype=float)
        
        x, y, k = 0, 0, 0
        delta_x = globalObj.Height/globalObj.nNodesH
        delta_y = globalObj.Width/globalObj.nNodesW
        
        
        for i in range(globalObj.nNodesW):
            x = i * delta_x
            
            for j in range(globalObj.nNodesH):
                y = j * delta_y
                
                self.nodes[k, 0] = x
                self.nodes[k, 1] = y
                k += 1




class Elements(object):
    width = 4 # number of nodes for grid's cells
    
    def __init__(self, height, nElems, nNodesH):
        self.height = height
        self.ID     = np.zeros((height, self.width), dtype=int)
        
        j = 0; r = nNodesH - 1
        
        for i in range(nElems):
            self.ID[i, 0] = j
            self.ID[i, 1] = self.ID[i, 0] + nNodesH
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
        self.nNodesH = 0
        self.nNodesW = 0
        self.nElems = 0
        self.nNodes = 0
        
        try:
            with open(filepath, 'r') as file:
                lines = file.readlines()
            
            self.Height = float(lines[0])
            self.Width  = float(lines[1])
            self.nNodesH = int(lines[2])
            self.nNodesW = int(lines[3])
            
        except (TypeError, FileNotFoundError) as err:
            print("Occured:", err)
        
        finally:
            file.close()
        
        
        self.nElems = int((self.nNodesH - 1)*(self.nNodesW - 1))
        self.nNodes = int(self.nNodesH * self.nNodesW)




class Elements4(object):
    
    def __init__(self):
        dN_deta = np.array([-1, -1, 1, 1])
        dN_dksi = np.array([-1, 1, 1, -1])
        diff    = np.array([1, -1, 1, -1])/ np.sqrt(3)
        
        self.integrPointsMatrixX = 0.25 * np.array([dN_dksi - diff, dN_dksi - diff, 
                                               dN_dksi + diff, dN_dksi + diff])
        self.integrPointsMatrixY = 0.25 * np.array([dN_deta - diff, dN_deta + diff,
                                               dN_deta + diff, dN_deta - diff])
        
        
    def calcHMatrix(self, element, index, k):
        dx_dksi = np.dot(element[:,0], self.integrPointsMatrixX[index, :])
        dy_dksi = np.dot(element[:,1], self.integrPointsMatrixX[index, :])
        dx_deta = np.dot(element[:,0], self.integrPointsMatrixY[index, :])
        dy_deta = np.dot(element[:,1], self.integrPointsMatrixY[index, :])
        
        self.jacobian = np.array([[dx_dksi, dy_dksi], [dx_deta, dy_deta]])
        self.detJacobian  = np.linalg.det(self.jacobian)
        self.resultMatrix = (1/self.detJacobian) * np.linalg.inv(self.jacobian)
        
        
        result = np.zeros((4,4))
        
        for i in range(len(self.integrPointsMatrixX[:, 0])):
            rowForX = self.resultMatrix[0,0] * self.integrPointsMatrixX[i, :]
            rowForY = self.resultMatrix[1,1] * self.integrPointsMatrixY[i, :]
            result  += np.power(self.detJacobian, 3) * (np.outer(rowForX, rowForX.T) + np.outer(rowForY, rowForY.T))
        
        return k * result




def quadratureCalc(fun, x1, x2, rank):
    period = 1./(1 -(-1))    # normalization
    det_J  = np.abs(x1 - x2) * period
    
    
    if rank == 2:
        Node1 = (1 + (1/np.sqrt(3))) * period
        Node2 = (1 - (1/np.sqrt(3))) * period
        
        integrationPoint1 = Node1*x1 + Node2*x2
        integrationPoint2 = Node1*x2 + Node2*x1
        
        return (fun(integrationPoint1) *1 + fun(integrationPoint2) *1)* det_J
    
    
    if rank == 3:
        Node1 = (1 + np.sqrt(3/5)) * period
        Node2 = period
        Node3 = (1 - np.sqrt(3/5)) * period
        
        integrationPoint1 = Node1*x1 + Node3*x2
        integrationPoint2 = Node2 * (x1 + x2)
        integrationPoint3 = Node1*x2 + Node3*x1
        
        return (fun(integrationPoint1) *(5/9) + fun(integrationPoint2) *(8/9) + fun(integrationPoint3) *(5/9))* det_J
    
    else:
        
        return "Cannot calculate this rank."
    