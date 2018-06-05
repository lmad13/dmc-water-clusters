import numpy as np
import matplotlib.pyplot as plt
import os
import usefulFunctions as use
import time


au2wn=219474.63
au2ang=0.529177249
massConversionFactor=1.000000000000000000/6.02213670000e23/9.10938970000e-28
ang2bohr=1.88973
bohr2ang=1.000/ang2bohr

verbose=False
class HarmonicApproxSpectrum(object):
    def __init__(self,wfn,coords,dw,path='dataH3O2/stateGround/DWGround'):
        self.wfn=wfn
        self.coords=coords
        self.dw=dw
        self.nVibs=self.wfn.molecule.nVibs
        self.path=path
    def LoadG(self,GfileName):
        print 'does ', GfileName, 'exist?'
        if not os.path.isfile(GfileName):
                print 'no!'
                gnm=self.calculateG(self.coords,self.dw)
                np.savetxt(GfileName,gnm)
        G=np.loadtxt(GfileName)
        return G


    def calculateG(self,eckartRotatedCoords,descendantWeights):
        #Input is x which is a NAtoms x 3(coordinates) sized array
        #input is also dx, the perturbation size, usually .001                                                
        #output is the G matrix, which is a self.nVibs*self.nVibs sized array (there are self.nVibs internals)

        dx=1e-4

        gnm=np.zeros((self.nVibs,self.nVibs))
        start=time.time()
        sumDescendants=0
        mass=self.wfn.molecule.get_mass()

        internal=self.wfn.molecule.SymInternals(eckartRotatedCoords)
        
        print 'some internals that you might care about!', np.average(internal,weights=descendantWeights,axis=0), '\n std',np.std(internal,axis=0)

        threwOut=0
        print 'summing up the descendants', np.sum(descendantWeights)
        sumDescendants=sumDescendants+np.sum(descendantWeights)
        for atom in range(self.wfn.molecule.nAtoms):
            for coordinate in range(3):
                print 'dx number',atom*3+(coordinate+1), 'atom:',atom, 'coordinate',coordinate
                deltax=np.zeros((eckartRotatedCoords.shape))
                deltax[:,atom,coordinate]=deltax[:,atom,coordinate]+dx #perturbs the x,y,z coordinate of the atom of interest                                        
                coordPlus=self.wfn.molecule.SymInternals(self.wfn.molecule.eckartRotate(eckartRotatedCoords+deltax))

                coordMinus=self.wfn.molecule.SymInternals(self.wfn.molecule.eckartRotate(eckartRotatedCoords-deltax))

                partialderv=(coordPlus-coordMinus)/(2.0*dx)

                timegnm=time.time()

                LastPartialDerv2MassWeighted=0
                for i,pd in enumerate(partialderv):
                    partialderv2=np.outer(pd,pd)
                    #print 'zeros ?',partialderv2prime[i*self.nVibs:(i+1)*self.nVibs,i*self.nVibs:(i+1)*self.nVibs]-partialderv2                                     
                    tempPartialDerv2MassWeighted=partialderv2*descendantWeights[i]/mass[atom]
                    
                    if np.any(tempPartialDerv2MassWeighted>1000000.0*dx):#(gnm[9,9]/(np.sum(self.Descendants[:i]))): $$$$$                                           
                        print 'atom',atom, 'coordinate',coordinate, i,'temp',np.transpose(np.where(tempPartialDerv2MassWeighted>10000.0*dx)),'is too big'
                        
                        print 'tempPartialDerv2MassWeighted', tempPartialDerv2MassWeighted, '\n Descendants', descendantWeights[i]
                        print 'coordinates \n', eckartRotatedCoords[i],'\n', eckartRotatedCoords[i]+deltax[i],'\n',eckartRotatedCoords[i]-deltax[i]
                        print 'eckart rotate \n', coordPlus[i], coordMinus[i]
                        print 'pd \n',pd

                        gnm=gnm+LastPartialDerv2MassWeighted
                        threwOut=threwOut+1
                    else:
                        gnm=gnm+tempPartialDerv2MassWeighted
                        LastPartialDerv2MassWEighted=1.0*tempPartialDerv2MassWeighted
            print 'gnmtiminging:',time.time()-timegnm
        print "THREW OUT ", threwOut, " coordinates :-("        
        print 'timing for G matrix', time.time()-start
        print 'dividing by ',sumDescendants
        gnm=gnm/sumDescendants
        return gnm

    def diagonalizeRootG(self,G):
        w,v=np.linalg.eigh(G)
        vinv=np.linalg.inv(v)
        invRootDiagG=np.diag(1.0/np.sqrt(w))  #1.0/rootDiagG  #MAYBE THIS IS OK??                                                                                       
        for i,ValinvRootDiagG in enumerate(invRootDiagG):#range(self.nVibs):                                                                                            
            if verbose: print 'alp[',i,']',ValinvRootDiagG[i]
        invRootG=np.dot(v,np.dot(invRootDiagG,v.transpose()))
        invG=np.dot(invRootG,invRootG.transpose())
        checkG=np.linalg.inv(invG)
        #self.GHalfInv=invRootG
        return  invRootG

    def calculateSecondMoments(self,x,dw):
        #calculate average internals
        internals=self.wfn.molecule.SymInternals(x)
        averageInternals=np.average(internals,weights=dw,axis=0)
        #calculate moments
        moments=internals-averageInternals

        #calculate second moments
        secondMomentsDWeighted=np.zeros((self.nVibs,self.nVibs))
        secondMoments=[]#np.zeros((self.nVibs,self.nVibs))
        for momp,desc in zip(moments,dw):
                secondMomentsDWeighted=secondMomentsDWeighted+np.outer(momp,momp)*desc
                secondMoments.append(np.outer(momp,momp))
        secondMomentsDWeighted=secondMomentsDWeighted/(2*float(np.sum(dw)))
        ###WHY AM I DIVIDING BY 2? because Anne did??####                                                                     
        # divide by 2 issue...it will change the alpha/eigenvalues by 2 but it won't be read in or used in the                 
        # future anyway so it doesn't actually matter.                                                                         
        secondMoments=np.array(secondMoments)

        return moments,secondMoments#,secondMomentsDWeighted


    def calculateSpectrum(self, coords,dw,GfileName):
        print 'I hope I am calculating the spectrum from ECKART ROTATED coordinates!'
        # eckRotcoords=self.molecule.eckartRotate(coords)
        self.G=self.LoadG(GfileName)
        

        #internals=self.molecule.SymInternals(coords)
        moments,secondMoments=self.calculateSecondMoments(coords,dw)

        q,q2,q4=self.calculateQCoordinates(moments,secondMoments,dw)
        q4ave=np.average(q4,axis=0,weights=dw)
        q2ave=np.average(q2,axis=0,weights=dw)
        qave =np.average(q,axis=0,weights=dw)
        if verbose: print '/\/\/\/\/\/\/\/\/\/\ '
        if verbose: print 'some averages', 
        if verbose: print 'q\n',qave
        if verbose: print 'q^2 \n',q2ave
        if verbose: print 'q^4 \n',q4ave
        if verbose: print '/\/\/\/\/\/\/\/\/\/\ '

        potentialEnergy=self.calculatePotentialEnergy(coords,dw)


        V_0=np.average(potentialEnergy[:,None],axis=0,weights=dw)
        Vq=np.average(potentialEnergy[:,None]*q2,axis=0,weights=dw)
        Vq2d=np.zeros((self.nVibs,self.nVibs))
        q2ave2d=np.zeros((self.nVibs,self.nVibs))
        for i in range(coords.shape[0]):
            Vq2d=Vq2d+np.outer(q2[i],potentialEnergy[i]*q2[i])*dw[i]
            q2ave2d=q2ave2d+np.outer(q2[i],q2[i])*dw[i]
        Vq2d=Vq2d/np.sum(dw)
        
    
        q2ave2d=q2ave2d/np.sum(dw)
        
        if verbose: print 'average v_0',V_0*au2wn
        if verbose: print 'Vq', Vq*au2wn

        #What's the difference with alpha??
        alpha=q2ave/(q4ave-q2ave**2)


        #kineticEnergy= hbar**2 nquanta alpha/(2 mass)
        Tq=1.0**2*1.0*alpha/(2.0*1.0)

        Tq2d=np.zeros((self.nVibs,self.nVibs))
        
        Eq2d=np.zeros((self.nVibs,self.nVibs))

        for ivibmode in range(self.nVibs):
            for jvibmode in range(ivibmode):
                Vq2d[ivibmode,jvibmode]=Vq2d[ivibmode,jvibmode]/q2ave2d[ivibmode,jvibmode]
                Vq2d[jvibmode,ivibmode]=Vq2d[jvibmode,ivibmode]/q2ave2d[jvibmode,ivibmode]
                Tq2d[ivibmode,jvibmode]=Tq[ivibmode]+Tq[jvibmode]
                Tq2d[jvibmode,ivibmode]=Tq2d[ivibmode,jvibmode]
                #Eq2d[ivibmode,jvibmode]=Vq2d[ivibmode,jvibmode]-V_0+Tq2d[ivibmode,jvibmode]
                #Eq2d[jvibmode,ivibmode]=Vq2d[jvibmode,ivibmode]-V_0+Tq2d[jvibmode,ivibmode]
            Vq2d[ivibmode,ivibmode]=Vq2d[ivibmode,ivibmode]*4.0*alpha[ivibmode]**2-Vq[ivibmode]*4.0*alpha[ivibmode]+V_0
            Vq2d[ivibmode,ivibmode]=Vq2d[ivibmode,ivibmode]/(4.0*alpha[ivibmode]**2*q2ave2d[ivibmode,ivibmode]-4.0*alpha[ivibmode]*q2ave[ivibmode]+1.0)
            Tq2d[ivibmode,ivibmode]=Tq[ivibmode]*2.0
            
        Eq2d=(Vq2d+Tq2d)-V_0
        Eq=(Vq/q2ave+Tq)-V_0
        print 'm alpha           Vq            Vq/q2ave      Tq            Eq'

        for i in range(9):
            print i, alpha[i], au2wn*Vq[i], au2wn*Vq[i], Vq[i]/q2ave[i]*au2wn, Eq[i]*au2wn
            
        if verbose: print 'energies!', zip(Vq2d[0]*au2wn,Tq2d[0]*au2wn,Eq2d[0]*au2wn)
        
        print 'Spectrum Info'
        fundamentalFileName=self.path+'fundamentals.data'
        comboFileName=self.path+'combinationOverrtoneBands.data'
        fundamentalFile=open(fundamentalFileName,'w')
        comboFile=open(comboFileName,'w')

        for i in range(9):
            print Eq[i]*au2wn, 'mode ',i 
            fundamentalFile.write(str( Eq[i]*au2wn)+"       "+str(i)+"\n")
        for i in range(9):
            for j in range(i):
                print Eq2d[i,j]*au2wn , 'combination bands' , i,j
                comboFile.write(str( Eq2d[i,j]*au2wn)+"       "+str(i)+" "+str(j)+"\n")
        
        fundamentalFile.close()
        comboFile.close
        

    def calculatePotentialEnergy(self,coords,dw):
        equilibriumEnergy=self.wfn.molecule.getEquilibriumEnergy()
        electronicEnergy=self.wfn.molecule.V(coords)
        relativePotentialEnergy=electronicEnergy-equilibriumEnergy
        #print 'average relative potential energy', np.average(relativePotentialEnergy)
        #print 'descendant weighted average relative potential energy', np.average(relativePotentialEnergy,weights=dw)
        return relativePotentialEnergy

    def calculateQCoordinates(self,moments,secondMoments,dw):
        mu2Ave=np.average(secondMoments,axis=0,weights=dw)
        mu2Ave=mu2Ave/2.000000000000
        GHalfInv=self.diagonalizeRootG(self.G)
        mu2AvePrime=np.dot(GHalfInv,np.dot(mu2Ave,GHalfInv))
        eigval,vects=np.linalg.eigh(mu2AvePrime)

        print 'diagonalized <mu^2>'
        for i in range(self.wfn.molecule.nVibs):
            if verbose: print 'v[',i,',]:',eigval[i]
            if verbose: print vects[:,i]

        #GAAH                                                              
        TransformationMatrix=np.dot(vects.transpose(),GHalfInv)
        #print 'RESULTS'
        Tnorm=1.0*TransformationMatrix # NEED TO DEEP COPY :-0                                                                    
        alpha=1.0*eigval #AS in alpha_j in equation 5 of the JPC A 2011 h5o2 dimer paper                                          
        #save the transformation matrix for future reference                                                                      
        TMatFileName=self.path+'TransformationMatrix.data'
        np.savetxt(TMatFileName,TransformationMatrix)


        for i,(vec,q2) in enumerate(zip(TransformationMatrix,eigval)):
            if verbose: print '\n',i,':<q^2 >=',q2
            #        alpha[i]=alpha[i]/(firstterm-(alpha[i]*alpha[i]))                                                            
            if verbose: print 'vec    ',vec
            Tnorm[i]=np.abs(vec)/np.max(np.abs(vec))
            sortvecidx=np.argsort(abs(vec))[::-1]
            sortedvec=vec[sortvecidx]
            sortedNames=use.sortList(abs(vec),self.wfn.molecule.internalName)[::-1]
            sortedvecLocalVersion=use.sortList(abs(vec),vec)[::-1]
            if verbose: print 'sorted\n',zip(sortvecidx,sortedvec,sortedNames,sortedvecLocalVersion)
            if verbose: print 'assignments\n',zip(sortedNames,sortvecidx,sortedvec)[0:4]



            #calculate q from the moments and the transformation matrix
        q=[]
        for s in moments:
            q.append(np.dot(TransformationMatrix,s))
            
        q=np.array(q)
        q2=q**2
        q4=q**4
        return q,q2,q4
