from fem.FEM import *




if __name__ == "__main__":
    
    global_data = GlobalData("data.txt")
    
    print(global_data.Height, global_data.Width, global_data.nNodesH, \
            global_data.nNodesW, global_data.nElems, global_data.nNodes)
    
    
    mocks = Mocks(global_data)
    print('Mesh nodes:\n', mocks.nodes)
    
    
    elemenets = Elements(global_data.nElems, global_data.nElems, global_data.nNodesH)
    print('Elements:\n', elemenets.ID, "\n\n")
    
    
    print('Testing quadrature function:', quadratureCalc(lambda x: 2* x**2 + 0.1* x + 3, 2.0, 10.0, 3), "\n")
    
    
    elem4 = Elements4()
    print(elem4.integrPointsMatrixX, '\n', elem4.integrPointsMatrixY, '\n')
    
    
    
    array = np.array([[0,0], [4,0], [4,6], [0,6]])
    H_matrix = elem4.calcHMatrix(array, 0, 30)
    print('Jacobian determinant:\n', elem4.detJacobian, '\nJacobian: \n', elem4.jacobian, '\n')
    print('H matrix:\n', H_matrix, '\n')
    
    
    