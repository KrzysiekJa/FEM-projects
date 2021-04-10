import numpy as np



class GlobalData(object):
    
    def __init__(self, filepath):
        self.Height = 0
        self.Width  = 0
        self.nNodesH = 0
        self.nNodesW = 0
        self.nElems = 0
        self.nNodes = 0
        self.kParameter = 0 # anisotropic temperature factor
        self.initTemp = 0
        self.simulationTime = 0
        self.timeStep = 0
        self.specificHeat = 0
        self.density = 0
        self.alpha = 0
        self.ta = 0
        
        try:
            with open(filepath, 'r') as file:
                lines = file.readlines()
            
            self.Height = float(lines[0])
            self.Width  = float(lines[1])
            self.nNodesH = int(lines[2])
            self.nNodesW = int(lines[3])
            self.kParameter = float(lines[4])
            self.initTemp = float(lines[5])
            self.simulationTime = float(lines[6])
            self.timeStep = float(lines[7])
            self.specificHeat = int(lines[8])
            self.density = int(lines[9])
            self.alpha = int(lines[10])
            self.tAlpha = int(lines[11])
            
        except (TypeError, FileNotFoundError) as err:
            print("Occured:", err)
        
        finally:
            file.close()
        
        
        self.nElems = int((self.nNodesH - 1)*(self.nNodesW - 1))
        self.nNodes = int(self.nNodesH * self.nNodesW)





class Nodes(object):
    width = 4  # point: (x, y)
    
    def __init__(self, globalObj):
        self.height = globalObj.nNodes
        self.nodes  = np.zeros((self.height, self.width), dtype=float)
        
        x, y, k = 0, 0, 0
        delta_x = globalObj.Height/(globalObj.nNodesH - 1)
        delta_y = globalObj.Width/(globalObj.nNodesW - 1)
        
        
        for i in range(globalObj.nNodesW):
            x = i * delta_x
            
            for j in range(globalObj.nNodesH):
                y = j * delta_y
                
                self.nodes[k, 0] = x
                self.nodes[k, 1] = y
                self.nodes[k, 2] = globalObj.initTemp
                
                if i == 0 or j == 0 or i%(globalObj.nNodesW-1) == 0 or j%(globalObj.nNodesH-1) == 0:
                    self.nodes[k, 3] = 1
                else:
                    self.nodes[k, 3] = 0
                k += 1





class Elements(object):
    width = 4  # number of nodes for grid's cells
    
    def __init__(self, height, nNodesH):
        self.height = height
        self.ID     = np.zeros((height, self.width), dtype=int)
        
        j = 0; r = nNodesH - 1
        
        for i in range(height):
            self.ID[i, 0] = j
            self.ID[i, 1] = self.ID[i, 0] + nNodesH
            self.ID[i, 2] = self.ID[i, 1] + 1
            self.ID[i, 3] = self.ID[i, 0] + 1
            j += 1
            
            if not (j % r):
                j += 1
                r += nNodesH
                
    
    
    def fillHCPMatrixTables(self, nodes, globalObj, rank):
        self.HMatrixTable = np.zeros((self.height, self.width, self.width))
        self.CMatrixTable = np.zeros((self.height, self.width, self.width))
        self.PMatrixTable = np.zeros((self.height, self.width))
        elem4 = Elements4(rank)
        nodesTable = np.zeros((self.width, 4))
               
        for i in range(self.height):
            for j in range(self.width):
                nodesTable[j] = nodes.nodes[self.ID[i,j],:]
                
            self.HMatrixTable[i] = elem4.calcHMatrix(nodesTable, globalObj.kParameter) \
                                 + elem4.calcHbcMatrix(nodesTable, globalObj.alpha)
            self.CMatrixTable[i] = elem4.calcCMatrix(nodesTable, globalObj.specificHeat, globalObj.density)
            self.PMatrixTable[i] = elem4.calcPMatrix(nodesTable, globalObj.alpha, globalObj.tAlpha)





class SOE(object):
    
    def __init__(self, elements, nNodes):
        self.HGlobal = np.zeros((nNodes, nNodes))
        self.CGlobal = np.zeros((nNodes, nNodes))
        self.PGlobal = np.zeros((nNodes))
        self.HCalc   = None
        
        
        for k in range(elements.height):
            for i in range(len(elements.ID[0])):
                self.PGlobal[elements.ID[k][i]] += elements.PMatrixTable[k][i].T
                for j in range(len(elements.ID[1])):
                    self.HGlobal[elements.ID[k][i]][elements.ID[k][j]] += elements.HMatrixTable[k][i][j]
                    self.CGlobal[elements.ID[k][i]][elements.ID[k][j]] += elements.CMatrixTable[k][i][j]
                    




class Elements4(object):
    
    def __init__(self, rank):
        self.rank = rank
        
        if self.rank == 2:
            self.sqrt = - 1/np.sqrt(3)
            ksi = np.zeros((rank * rank))
            eta = np.zeros((rank * rank))
            
            for i in range(rank):
                for j in range(rank):
                    ksi[i*rank + j] = self.sqrt - 2 * j * self.sqrt
                    eta[i*rank + j] = self.sqrt - 2 * i * self.sqrt
            
            self.shapeFunctions = np.zeros((rank * rank, rank * rank))
            for i in range(rank * rank):
                self.shapeFunctions[i][0] = self.calculateNode1(ksi[i], eta[i])
                self.shapeFunctions[i][1] = self.calculateNode2(ksi[i], eta[i])
                self.shapeFunctions[i][2] = self.calculateNode3(ksi[i], eta[i])
                self.shapeFunctions[i][3] = self.calculateNode4(ksi[i], eta[i])
                
            dN_dksi = np.array([-1, 1, 1, -1])
            dN_deta = np.array([-1, -1, 1, 1])
            diff    = np.array([1, -1, 1, -1])* self.sqrt
        
            self.derivativesKsi = 0.25 * np.array([dN_dksi - diff, dN_dksi - diff,
                                                   dN_dksi + diff, dN_dksi + diff])
            self.derivativesEta = 0.25 * np.array([dN_deta - diff, dN_deta + diff,
                                                   dN_deta + diff, dN_deta - diff])
        
        if self.rank == 3:
            dN_deta = np.array([-1, -1, 1, 1])
            dN_dksi = np.array([-1, 1, 1, -1])
            diff    = np.array([1, -1, 1, -1])* np.sqrt(3/5)
        
            self.derivativesKsi = 0.25 * np.array([dN_dksi - diff, dN_dksi - diff,
                                                   dN_dksi + diff, dN_dksi + diff])
            self.derivativesMid = 0.25 * np.array([dN_dksi, dN_dksi, dN_dksi, dN_dksi])
            self.derivativesEta = 0.25 * np.array([dN_deta - diff, dN_deta + diff,
                                                   dN_deta + diff, dN_deta - diff])
        
        
        
    def calcHMatrix(self, element, kParameter):
        if not ('rank' in dir(self)):
            raise AttributeError("Rank artibute not defined.")
        
        if self.rank == 2:
            dx_dksi = np.dot(element[:,0], self.derivativesKsi[0, :])
            dy_dksi = np.dot(element[:,1], self.derivativesKsi[0, :])
            dx_deta = np.dot(element[:,0], self.derivativesEta[0, :])
            dy_deta = np.dot(element[:,1], self.derivativesEta[0, :])
            self.jacobian = np.array([[dx_dksi, dy_dksi], [dx_deta, dy_deta]])
        
        
        if self.rank == 3:
            mid = (element[:,0] + element[:,1])*0.5
            dx_dksi = np.dot(element[:,0], self.derivativesKsi[0, :])
            dm_dksi = np.dot(mid, self.derivativesKsi[0, :])
            dy_dksi = np.dot(element[:,1], self.derivativesKsi[0, :])
            dx_dmid = np.dot(element[:,0], self.derivativesMid[0, :])
            dm_dmid = np.dot(mid, self.derivativesMid[0, :])
            dy_dmid = np.dot(element[:,1], self.derivativesMid[0, :])
            dx_deta = np.dot(element[:,0], self.derivativesEta[0, :])
            dm_deta = np.dot(mid, self.derivativesEta[0, :])
            dy_deta = np.dot(element[:,1], self.derivativesEta[0, :])
            self.jacobian = np.array([[dx_dksi, dm_dksi, dy_dksi], [dx_dmid, dm_dmid, dy_dmid], [dx_deta, dm_deta, dy_deta]])
        
        
        self.detJacobian  = np.linalg.det(self.jacobian)
        self.resultMatrix = (1/self.detJacobian) * np.linalg.inv(self.jacobian)
        
        lenght = len(self.derivativesKsi[:, 0])
        
        result = np.zeros((lenght, lenght))
        
        if self.rank == 2:
            for i in range(lenght):
                rowForX = self.resultMatrix[0,0] * self.derivativesKsi[i, :]
                rowForY = self.resultMatrix[1,1] * self.derivativesEta[i, :]
                result  += np.power(self.detJacobian, 3) * (np.outer(rowForX, rowForX.T) + np.outer(rowForY, rowForY.T))
                        
        if self.rank == 3:
            for i in range(lenght):
                rowForX   = self.resultMatrix[0,0] * self.derivativesKsi[i, :]
                rowForMid = self.resultMatrix[1,1] * self.derivativesMid[i, :]
                rowForY   = self.resultMatrix[2,2] * self.derivativesEta[i, :]
                result  += np.power(self.detJacobian, 3) * ((5/9)* np.outer(rowForX, rowForX.T)\
                                    + (8/9)* np.outer(rowForY, rowForY.T) + (5/9)* np.outer(rowForY, rowForY.T))
        
        return kParameter * result
    
    
    
    def calcCMatrix(self, element, heat, density):
        if not ('rank' in dir(self)):
            raise AttributeError("Rank artibute not defined.")
        
        lenght = len(self.derivativesKsi[:, 0])
        result = np.zeros((lenght, lenght))
        
        jacobianC   = (np.array([self.calculateJacobian(i, element, self.derivativesKsi, self.derivativesEta) for i in range(lenght)]))
        jacobianDet = np.array([np.linalg.det(jacobianC[i]) for i in range(lenght)])
        
        if self.rank == 2:
            for i in range(lenght):
                result += jacobianDet[i] * np.outer(self.shapeFunctions[i].T, self.shapeFunctions[i])
        
        return result * heat * density
    
    
    
    def calcHbcMatrix(self, Nodes, alpha):
        if not ('rank' in dir(self)):
            raise AttributeError("Rank artibute not defined.")
        
        resultHbc = np.zeros((4,4))
        
        if self.rank == 2:
            if Nodes[0, 3] == 1 and Nodes[1, 3] == 1:
                resultHbc += self.calculateHbc(Nodes[0], Nodes[1], [[self.sqrt, -1],[-self.sqrt, -1]], alpha)
            if Nodes[1, 3] == 1 and Nodes[2, 3] == 1:
                resultHbc += self.calculateHbc(Nodes[1], Nodes[2], [[1, self.sqrt],[1, -self.sqrt]], alpha)
            if Nodes[2, 3] == 1 and Nodes[3, 3] == 1:
                resultHbc += self.calculateHbc(Nodes[2], Nodes[3], [[self.sqrt, 1],[-self.sqrt, 1]], alpha)
            if Nodes[3, 3] == 1 and Nodes[0, 3] == 1:
                resultHbc += self.calculateHbc(Nodes[3], Nodes[0], [[-1, self.sqrt],[-1, -self.sqrt]], alpha)
            
        return resultHbc
    
    
    
    def calcPMatrix(self, Nodes, alpha, tAlpha):
        if not ('rank' in dir(self)):
            raise AttributeError("Rank artibute not defined.")
        
        result = np.zeros((1, 4))
        
        if self.rank == 2:
            if Nodes[0, 3] == 1 and Nodes[1, 3] == 1:
                result += self.calculateP(Nodes[0], Nodes[1], [[self.sqrt, -1],[-self.sqrt, -1]], alpha, tAlpha)
            if Nodes[1, 3] == 1 and Nodes[2, 3] == 1:
                result += self.calculateP(Nodes[1], Nodes[2], [[1, self.sqrt],[1, -self.sqrt]], alpha, tAlpha)
            if Nodes[2, 3] == 1 and Nodes[3, 3] == 1:
                result += self.calculateP(Nodes[2], Nodes[3], [[self.sqrt, 1],[-self.sqrt, 1]], alpha, tAlpha)
            if Nodes[3, 3] == 1 and Nodes[0, 3] == 1:
                result += self.calculateP(Nodes[3], Nodes[0], [[-1, self.sqrt],[-1, -self.sqrt]], alpha, tAlpha)
            
        return result





    def calculateNode1(self, ksi, eta):
        return 0.25* (1 - ksi) * (1 - eta)
    
    def calculateNode2(self, ksi, eta):
        return 0.25* (1 + ksi) * (1 - eta)
    
    def calculateNode3(self, ksi, eta):
        return 0.25* (1 + ksi) * (1 + eta)
    
    def calculateNode4(self, ksi, eta):
        return 0.25* (1 - ksi) * (1 + eta)
    
    def calculateJacobian(self, iterNum, element, derivativesKsi, derivativesEta):
        # Jacobian for 2D
        
        result = np.zeros((2,2))
        for i in range(4):
            result[0][0] += self.derivativesKsi[iterNum][i] * element[i, 0]
            result[0][1] += self.derivativesKsi[iterNum][i] * element[i, 1]
            result[1][0] += self.derivativesEta[iterNum][i] * element[i, 0]
            result[1][1] += self.derivativesEta[iterNum][i] * element[i, 1]

        return result


    def calculateHbc(self, node1, node2, vec, alpha):
        dist   = 0.5 * np.linalg.norm(node1[:2] - node2[:2])
        result = np.zeros((4, 4))
        values = np.zeros((len(vec), 4))
    
        for i in range(len(vec)):
            values[i][0] = self.calculateNode1(vec[i][0], vec[i][1])
            values[i][1] = self.calculateNode2(vec[i][0], vec[i][1])
            values[i][2] = self.calculateNode3(vec[i][0], vec[i][1])
            values[i][3] = self.calculateNode4(vec[i][0], vec[i][1])
    
        for i in range(len(vec)):
            result += np.matmul(values.T[:, i:i+1], values[i:i+1, :])
    
        return result * alpha * dist


    def calculateP(self, node1, node2, vec, alpha, tAlpha):
        dist   = 0.5 * np.linalg.norm(node1[:2] - node2[:2])
        result = np.zeros((4, 1))
        values = np.zeros((len(vec), 4))
    
        for i in range(len(vec)):
            values[i][0] = self.calculateNode1(vec[i][0], vec[i][1])
            values[i][1] = self.calculateNode2(vec[i][0], vec[i][1])
            values[i][2] = self.calculateNode3(vec[i][0], vec[i][1])
            values[i][3] = self.calculateNode4(vec[i][0], vec[i][1])
    
        for i in range(len(vec)):
            result += values.T[:, i:i+1]
    
        return (-alpha * result * tAlpha * dist).T








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
