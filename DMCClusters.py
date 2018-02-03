import numpy as np
import matplotlib.pyplot as plt
au2wn=219474.63
au2ang=0.529177249


class wavefunction:

    def __init__(self,molecule, nWalkers):
        self.moleculeDict=self.loadDict("molecule.dict")
        self.pathDict=self.loadDict("paths.dict")
        print 'initialized', nWalkers,' coordinates for ', molecule


        #based on the molecule find the coordinates
        self.x=self.loadCoordinates(molecule,nWalkers)

        

    def loadCoordinates(self,molecule,nWalkers):
        #looks up the path of the coordinates for the starting positions of the walkers and makes nWalkers copies of them and then returns that as an self.nWalkers,self.nAtoms, 3 dimensional array
        self.nAtoms=self.moleculeDict["nAtoms"+molecule]
        print 'there are ', self.nAtoms , 'atoms in ', molecule

        return 1

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

