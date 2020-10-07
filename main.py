import numpy as np
from FEM import Mocks, Elements, GlobalData




if __name__ == "__main__":
    
    global_data = GlobalData("data.txt")
    
    print(global_data.Height, global_data.Width, global_data.numH, \
            global_data.numW, global_data.nElems, global_data.nNodes)
    
    
    mocks = Mocks(global_data.nNodes, global_data.numH, global_data.numW, global_data.Height, global_data.Width)
    
    print(mocks.nodes)
    
    
    elemenets = Elements(global_data.nElems, global_data.nElems, global_data.numH)
    
    print(elemenets.ID)