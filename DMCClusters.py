import sys
import os
import numpy as np
import matplotlib.pyplot as plt
au2wn=219474.63
au2ang=0.529177249


class wavefunction:

    def __init__(self,molecule, nWalkers):
        self.molecule=molecule
        self.moleculeDict=self.loadDict("molecule.dict")
        self.pathDict=self.loadDict("paths.dict")
        print 'initialized', nWalkers,' coordinates for ', self.molecule
        self.nAtoms=int(self.moleculeDict["nAtoms"+self.molecule])
        self.potential=self.loadPES()

        #based on the molecule find the coordinates
        self.x=np.zeros((nWalkers,self.nAtoms,3))
        self.x[:]=self.loadCoordinates(self.molecule,nWalkers)
        
    
    def loadPES(self):
        sys.path.insert(0,self.pathDict["potentialPath"+self.molecule])
        import pes 
        return pes
        
    def loadCoordinates(self,molecule,nWalkers):
        #looks up the path of the coordinates for the starting positions of the walkers and makes nWalkers copies of them and then returns that as an self.nWalkers,self.nAtoms, 3 dimensional array
        
        print 'there are ', self.nAtoms , 'atoms in ', self.molecule
        coords=np.array([  
                [   0.000000000000000 ,  0.000000000000000 ,  0.000000000000000 ],
                [  -2.306185590098382 ,  0.000000000000000 ,  0.000000000000000],
                [  -2.749724314110769 ,  1.765018349357672 ,  0.000000000000000],
                [   2.306185590098382 ,  0.000000000000000 ,  0.000000000000000],
                [   2.749724314110769 ,  1.765018349357672 ,  0.000000000000000]
                ])
        
        return coords

    def loadDict(self,fileName):
        fileIn=open(fileName,'r')
        pathDict={}
        for line in fileIn:
            [key,element]=line.split()
            pathDict[key]=element
        return pathDict


    
#wavefunction
#initialize

#propagate

#calculate V

