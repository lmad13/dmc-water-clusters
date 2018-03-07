import numpy as np
import sys
import os
global massH
global massO
massH=1.00782503223
massO=15.99491561957
massConversionFactor=1.000000000000000000/6.02213670000e23/9.10938970000e-28#1822.88839 
global Water
Water={'H2O', 'water', 'HOH', 'h2o'}
global WaterDimer
WaterDimer = {'water dimer', 'H4O2'}
global ProtonatedWaterDimer
ProtonatedWaterDimer = {'H5O2+','HHOHOHH','H5O2+','h5o2plus','h5o2'}
global DeprotonatedWaterDimer
DeprotonatedWaterDimer = {'HOHOH','H3O2-', 'h3o2-', 'h3o2', 'H3O2'}
ProtonatedWaterTrimer = {'H7O3+','O3H7+', 'H7O3plus','H7O3', 'O3H7'}
ProtonatedWaterTetramer = {'H9O4+','O4H9+', 'H9O4plus','H9O4', 'O4H9'}


class molecule:
    def __init__(self,moleculeName):
        self.name=moleculeName
        if self.name in DeprotonatedWaterDimer:
            self.nAtoms=5
        self.pathDict=self.loadDict("paths.dict")
        self.potential=self.getPES()
        self.surfaceName="GroundState"
        self.dipole=self.getDIPOLE()
    def setNodalSurface(self,surfaceName,side):
        self.surfaceName=surfaceName
        self.side=side
        self.state=1
    def getPES(self):
        sys.path.insert(0,self.pathDict["potentialPath"+self.name])
        import pes
        if self.name in DeprotonatedWaterDimer:
            pes.prepot()
            potential=pes.getpot
            print 'potential retreived for HOHOH. Be sure to feed one walker in at a time!'
        return potential

    def getDIPOLE(self):
        sys.path.insert(0,self.pathDict["dipolePath"+self.name])
        import h3o2dms2 as dms
        if self.name in DeprotonatedWaterDimer:
            dms.predip()
            dip=dms.mycalcdip
        #since all dipole moment calculations need eckart rotation...loading the reference coordinates
        self.refPos=self.loadRefCoord()
        return dip
            
    def loadDict(self,fileName):
        fileIn=open(fileName,'r')
        pathDict={}
        for line in fileIn:
            [key,element]=line.split()
            pathDict[key]=element
        return pathDict
    def loadRefCoord(self):
        if self.name in DeprotonatedWaterDimer:
            coords=np.array([[ 0.2981678882048853 ,   -2.4557992072743176E-002,  -5.5485232545510215E-002],
                             [ -2.354423404994569 ,     0.000000000000000     ,    0.000000000000000],
                             [ -2.858918674095194 ,     1.111268022307282     ,   -1.352651141853729],
                             [  2.354423404994569 ,     0.000000000000000     ,    0.000000000000000],
                             [  2.671741580470489 ,     1.136107563104921     ,    1.382886181959795]])
        return coords

    def V(self,x):
        v=np.array([self.potential(cart_in) for cart_in in x])
        #for n,pot_mol in enumerate(v):
        #    if np.isnan(pot_mol):
        #        print 'gosh, the ', n,'th is a problem:\n', x[n]

        if self.surfaceName=='SharedProton':
            r=self.calcRn(x)
            if self.side=='Right':
                v[(r<0)]=1000.00
            elif self.side=='Left':
                v[(r>0)]=1000.00
            else:
                donothing=1
        return v

    def calcDipole(self,x):
        eckRotx=self.eckartRotate(x)
        dipoleVectors=self.dipole(eckRotx)
        
        return dipoleVectors


    def eckartRotate(self,pos,specialCond=False):
        if len(pos.shape)<3:
            pos=np.array([pos])
        nMolecules=pos.shape[0]
        Fvec=np.zeros((3,3))
        Fvec2=np.zeros((3,3))
        newCoord=np.zeros(pos.shape)
        rHC=np.zeros((nMolecules,3))
        rHCprime=np.zeros((nMolecules,3))
    
        #Center of Mass     
        mass=self.get_mass()
        com=np.dot(mass,pos)/np.sum(mass)
        rs=0.000
        rsprime=0.000

        #First Translate:                                                                                     
        ShiftedMolecules=pos-com[:,np.newaxis,:]

        #Equation 3.1                                                                                         
        for moli,molecule in enumerate(ShiftedMolecules):
            Fvec=np.zeros((3,3))
            for atom,massa,eckatom in zip(molecule,mass,self.refPos):
                Fvec=Fvec+massa*np.outer(eckatom,atom)
            #F from eqn 3.4b                                                                                  
            FF=np.dot(Fvec,Fvec.transpose())
            #Diagonalize FF                                                                                   
            sortEigValsF,sortEigVecF=np.linalg.eigh(FF)
            sortEigVecFT=-sortEigVecF.transpose()
            if specialCond:
                print 'special condition activated!'
                print 'eigenvals \n',sortEigValsF
                print 'vect \n',sortEigVecFT
            if len(np.where(sortEigValsF<=0)[0])!=0:
                #sortEigVecFT=np.abs(sortEigVecFT)                                                            
                sortEigValsF=np.abs(sortEigValsF)
                invRootDiagF=sortEigValsF
                for e,element in enumerate(sortEigValsF):
                    if element>0:
                        invRootDiagF[e]=1.0/np.sqrt(element)
            #Get the inverse sqrt of diagonalized(FF)                                                         
            else:
                invRootDiagF=1.0/np.sqrt(sortEigValsF)
            # F^{-1/2}                                                                                         
            invRootF=np.dot(invRootDiagF[np.newaxis,:]*-sortEigVecF,sortEigVecFT)
            eckVecs=np.dot(Fvec.transpose(),invRootF)
            newCoord[moli]= np.dot(molecule,eckVecs)

            if len(np.where(np.isnan(newCoord[moli]))[0])!=0:
                print 'whaaaaaT?! nan',np.where(np.isnan(newCoord[moli])),'\ncoords:\n',newCoord[moli]
                print '   molecule number:',moli,'\n   sortEigValsF: \n', sortEigValsF,'\n   molecule: \n', molecule, 
                print '\n   eckVecs \n', eckVecs
                octopus

        return newCoord


    def getInitialCoordinates(self):
        #looks up the path of the coordinates for the starting positions of the walkers and makes nWalkers copies of them and then returns that as an self.nWalkers,self.nAtoms, 3 dimensional array                      
        print 'there are ', self.nAtoms , 'atoms in ', self.name
        if self.name in DeprotonatedWaterDimer:
            coords=np.array([
                    [   0.000000000000000 ,  0.000000000000000 ,  0.000000000000000],
                    [  2.306185590098382 ,  0.000000000000000 ,  0.000000000000000],
                    [  2.749724314110769 ,  -1.765018349357672 ,  0.000000000000000],
                    [   -2.306185590098382 ,  0.000000000000000 ,  0.000000000000000],
                    [   -2.749724314110769 ,  -1.765018349357672 ,  0.000000000000000]
                    ])
            coords=np.array([
                    [0.0000000,   0.000000000000000 ,  0.000000000000000  ],
                    [0.0000000,  2.306185590098382 ,  0.000000000000000   ],
                    [0.0000000,  2.749724314110769 ,  -1.765018349357672  ],
                    [0.0000000,   -2.306185590098382 ,  0.000000000000000 ],
                    [0.0000000,   -2.749724314110769 ,  -1.765018349357672]
                    ])
#            coords=np.array([[ 0.26591125,  0.07072797, -0.02256279],
#WW         WW                [ 0.55610034,  2.83109547,  0.14883552],
#  W       W                  [ 1.50122114,  1.22631416, -0.59092507],
#   W  W  W                   [-0.11962985, -1.87021212,  0.22794889],
#    W   W                    [ 1.20503929, -0.77837156,  0.71051114]])
            
            self.names=['H','O','H','O','H']
        else:
            print 'not implemented!!'
        return coords

    def get_mass(self):
        mass=np.zeros((self.nAtoms))
        if self.name in DeprotonatedWaterDimer:
            mass=np.array([massH,massO,massH,massO,massH])
        return mass*massConversionFactor

    def calcReducedmass(self,x):
        # for the shared proton stretch, page 42 in notebook #2
        massOamu=massO*massConversionFactor
        massHamu=massH*massConversionFactor
        if self.name in ProtonatedWaterDimer and self.state==1 and self.surfaceName=='SharedProton':
            ##COMwat1=1.0/massWater*(massO*x[:,0,:]+massH*(x[:,3,:]+x[:,4,:]))
            ##COMwat2=1.0/massWater*(massO*x[:,1,:]+massH*(x[:,5,:]+x[:,6,:]))
            ##U=x[:,2,:]-COMwat1                                              
            ##V=x[:,2,:]-COMwat2                                                                                                  
            U=x[:,2,:]-x[:,1,:]
            V=x[:,2,:]-x[:,0,:]

            magU=np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2)
            magV=np.sqrt(V[:,0]**2+V[:,1]**2+V[:,2]**2)
            
            costheta= np.diag(np.dot(U,V.T))/(magU*magV)
            #mass=1.0/(2.0*((1.0/(massOamu+massHamu+massHamu))+((1-costheta)/(massHamu))))                                        
            #corresponds to calcrncom                                                                                             

            mass=1.0/(2.0*((1.0/(massOamu))+((1-costheta)/(massHamu))))

            #Mass of water or Mass of O??                                                                                         
        elif self.name in ProtonatedWaterDimer and self.state==1 and self.surfaceName=='StretchAntiIn':

            U=x[:,0,:]-x[:,3,:]
            V=x[:,0,:]-x[:,4,:]

            magU=np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2)
            magV=np.sqrt(V[:,0]**2+V[:,1]**2+V[:,2]**2)
            costhetaWat1= np.diag(np.dot(U,V.T))/(magU*magV)

            U=x[:,1,:]-x[:,5,:]
            V=x[:,1,:]-x[:,6,:]
            magU=np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2)
            magV=np.sqrt(V[:,0]**2+V[:,1]**2+V[:,2]**2)
            costhetaWat2= np.diag(np.dot(U,V.T))/(magU*magV)

            g=( (1.0/massOamu)+(1.0/massHamu) )  +  1.0/2.0*((costhetaWat1/massOamu)+(costhetaWat2/massOamu))

            mass= 1.0/g

        elif self.name in DeprotonatedWaterDimer and self.state==1:

            U=x[:,0,:]-x[:,1,:]
            V=x[:,0,:]-x[:,2,:]

            magU=np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2)
            magV=np.sqrt(V[:,0]**2+V[:,1]**2+V[:,2]**2)
            costheta= np.diag(np.dot(U,V.T))/(magU*magV)
            mass=1.0/(2.0*((1.0/(massOamu))+((1-costheta)/(massHamu))))

        elif self.name in ProtonatedWaterDimer and self.state==0:
            m2=2*(massO)*conversionFactor
            m1=massH*conversionFactor
            mass=m1*m2/(m1+m2)
            print 'why are you calculateing the reduced mass on the ground state?'  , end
        elif self.name in DeprotonatedWaterDimer and self.state==0:
            m2=2*(massO)*conversionFactor
            m1=massH*conversionFactor
            mass=m1*m2/(m1+m2)
            print 'why are you calculateing the reduced mass on the ground state?'  , end

        else:
            print 'not implemented for ', self.name , 'and', self.state, 'and', self.surfaceName,end

        return mass
    
    def calcSharedProtonDisplacement(self,x):
        if self.name in ProtonatedWaterDimer:
            r1=self.bondlength(x,atom1=2, atom2=1)
            r2=self.bondlength(x,atom1=2, atom2=0)
            return r2-r1
        elif self.name in DeprotonatedWaterDimer:
            r1=self.bondlength(x,atom1=0, atom2=1)
            r2=self.bondlength(x,atom1=0, atom2=3)#ha                          
            return r2-r1

    def calcStretchAntiIn(self,x):
        if self.name in ProtonatedWaterDimer:
            r1=self.bondlength(x,atom1=0, atom2=3)
            r2=self.bondlength(x,atom1=0, atom2=4)
            r3=self.bondlength(x,atom1=1, atom2=5)
            r4=self.bondlength(x,atom1=1, atom2=6)
            return 0.5*(r1+r2-r3-r4)

    def calcRn(self,x):
        if self.surfaceName=='SharedProton':
            return self.calcSharedProtonDisplacement(x)
        elif self.surfaceName=='StretchAntiIn':
            return self.calcStretchAntiIn(x)


    def bondlength(self,pos,atom1=1,atom2=2):
        length=(pos[:,atom1,0]-pos[:,atom2,0])**2+(pos[:,atom1,1]-pos[:,atom2,1])**2+(pos[:,atom1,2]-pos[:,atom2,2])**2
        length=np.sqrt(length)
        return length

    def mag(self,xList):
        magnitude=np.zeros(xList.shape[0])
        for i,x in enumerate(xList):
            magnitude[i]=np.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
        return magnitude

