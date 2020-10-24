from fem.FEM import *




if __name__ == "__main__":
    
    global_data = GlobalData("data.txt")
    
    print(global_data.Height, global_data.Width, global_data.nNodesH, \
            global_data.nNodesW, global_data.nElems, global_data.nNodes)
    
    
    mocks = Mocks(global_data)
    
    print(mocks.nodes)
    
    
    elemenets = Elements(global_data.nElems, global_data.nElems, global_data.nNodesH)
    
    print(elemenets.ID)
    
    
    print(quadratureCalc(lambda x: 2* x**2 + 0.1* x + 3, 2.0, 10.0, 3))
    
    
    elem4 = Elements4(mocks.nodes[4:8])