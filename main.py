from fem.FEM import *
from fem.simulation import simulate




if __name__ == "__main__":
    
    global_data = GlobalData("data_for_small_simulation.txt")
    rank = 2
    
    print(global_data.Height, global_data.Width, global_data.nNodesH, \
            global_data.nNodesW, global_data.nElems, global_data.nNodes, global_data.kParameter)
    
    nodes = Nodes(global_data)
    with np.printoptions(precision=3, suppress=True):
        print('Mesh nodes:\n', nodes.nodes)
    
    elements = Elements(global_data.nElems, global_data.nNodesH)
    print('\nElements:\n', elements.ID, "\n")
    print("\n--------Testing part --------\n")
    
    
    #print('Testing quadrature function:', quadratureCalc(lambda x: 2* x**2 + 0.1* x + 3, 2.0, 10.0, rank), "\n")
    
    
    
    #elem4 = Elements4(rank)
    #array = np.array([[0,0], [4,0], [4,6], [0,6]])
    #H_matrix = elem4.calcHMatrix(array, 30)
    #print('Jacobian determinant:\n', elem4.detJacobian, '\nJacobian: \n', elem4.jacobian, '\n')
    #print('H matrix:\n', H_matrix, '\n')
    
    
    print("\n-------Simulation part-------\n")
    elements.fillHCPMatrixTables(nodes, global_data, rank)
    
    soe = SOE(elements, global_data.nNodes)
    #with np.printoptions(precision=3, suppress=True):
    #    print(soe.HGlobal, '\n\n')
    #    print(soe.CGlobal, '\n\n')
    #    print(soe.PGlobal, '\n\n')
    
    
    simulate(global_data, soe, nodes)
    