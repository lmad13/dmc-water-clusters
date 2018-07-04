
def printArray(a):
    string=[]
    for row in a:
        string.append(str(row).replace('[',' ').replace(']',' '))
    return string


def sortList(seq,slist):
    return [slist[i] for (v,i) in sorted((v,i) for (i,v) in enumerate(seq))]


def takeDerivative(x,y,derv,noisy=True,smoothOver=5):
    import matplotlib.pyplot as plt
    import numpy as np

    
    if noisy:

        ptSmoothing=smoothOver
        length=x.shape[0]
        newx=np.zeros(length-ptSmoothing)
        newy=np.zeros(length-ptSmoothing)

        for i in range(length-ptSmoothing):
            newx[i]=np.average(x[i:i+ptSmoothing])
            newy[i]=np.average(y[i:i+ptSmoothing])
        x=newx
        y=newy
        avg=np.average(newy)
        

    dx=newx[1]-newx[0]

    if derv==0:
        return newx, newy

    elif derv==1:
        return newx,np.gradient(newy,dx,edge_order=1)

    elif derv==2:
        return newx, np.gradient(np.gradient(newy,dx,edge_order=1),dx,edge_order=1)

