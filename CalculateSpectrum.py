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
        #Equations are referenced from McCoy, Diken, and Johnson. JPC A 2009,113,7346-7352
        print 'I better be calculating the spectrum from ECKART ROTATED coordinates!'

        #First, What is the G Matrix for this set of walkers based on the SymInternals coordinates
        self.G=self.LoadG(GfileName)
        

        # Now determine the intrinsic linear combinations of SymInternal coordinates of the Wfn based on the second moments
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

        #Now calculate the Potential energy
        potentialEnergy=self.calculatePotentialEnergy(coords,dw)
        # V_0=<0|V|0>
        V_0=np.average(potentialEnergy[:,None],axis=0,weights=dw)

        # Vq= <1_q|V|0> bra:state with one quanta in mode q, ket: ground state
        Vq=np.average(potentialEnergy[:,None]*q2,axis=0,weights=dw)

        #Vq2d=<q^2Vq^2> (<> is an average of the descendant weights)
        Vq2d=np.zeros((self.nVibs,self.nVibs))

        #q2ave2d=<q^2 q^2> (<> is an average of the descendant weights)  
        q2ave2d=np.zeros((self.nVibs,self.nVibs))

        for i in range(coords.shape[0]):
            Vq2d=Vq2d+np.outer(q2[i],potentialEnergy[i]*q2[i])*dw[i]
            q2ave2d=q2ave2d+np.outer(q2[i],q2[i])*dw[i]
        Vq2d=Vq2d/np.sum(dw)
        q2ave2d=q2ave2d/np.sum(dw)        
        
        #Now calculate the kinetic energy
        print 'ZPE: average v_0',V_0*au2wn
        print 'Vq', Vq*au2wn

        alpha=q2ave/(q4ave-q2ave**2) # Equation #11 
        alphaPrime=0.5/q2ave   #Equation in text after #8
        #        print 'how similar are these?', zip(alpha,alphaPrime) Still a mystery to my why there were 2 dfns of alpha


        #kineticEnergy= hbar**2 nquanta alpha/(2 mass)
        Tq=1.0**2*1.0*alpha/(2.0*1.0) #Equation #10
        
        Tq2d=np.zeros((self.nVibs,self.nVibs)) #Tij=Ti+Tj
        #Finish calculate the potential and kinetic energy for the combination and overtone bands 
        #put Equation 8 into equation 4, algebra...
        for ivibmode in range(self.nVibs):
            for jvibmode in range(ivibmode):
                Vq2d[ivibmode,jvibmode]=Vq2d[ivibmode,jvibmode]/q2ave2d[ivibmode,jvibmode]
                Vq2d[jvibmode,ivibmode]=Vq2d[jvibmode,ivibmode]/q2ave2d[jvibmode,ivibmode]
                Tq2d[ivibmode,jvibmode]=Tq[ivibmode]+Tq[jvibmode]
                Tq2d[jvibmode,ivibmode]=Tq2d[ivibmode,jvibmode]
            Vq2d[ivibmode,ivibmode]=Vq2d[ivibmode,ivibmode]*4.0*alpha[ivibmode]**2-Vq[ivibmode]*4.0*alpha[ivibmode]+V_0
            Vq2d[ivibmode,ivibmode]=Vq2d[ivibmode,ivibmode]/(4.0*alpha[ivibmode]**2*q2ave2d[ivibmode,ivibmode]-4.0*alpha[ivibmode]*q2ave[ivibmode]+1.0)
            Tq2d[ivibmode,ivibmode]=Tq[ivibmode]*2.0
            
        #Energy is V+T, don't forget to subtract off the ZPE
            #Eq2d=np.zeros((self.nVibs,self.nVibs))
        Eq2d=(Vq2d+Tq2d)-V_0
        Eq=(Vq/q2ave+Tq)-V_0
        
        #Lastly, the dipole moments to calculate the intensities!
        dipoleMoments=self.wfn.molecule.dipole(coords)
        averageDipoleMoment_0=np.average(dipoleMoments,axis=0,weights=dw)

        #aMu_i=<1_i|dipoleMoment|0>=Sum_n^walkers(Dipole_n*q_i*dw_n)/Sum(dw) #Equation 4 for \psi_n=\psi_1
        aMu=np.zeros((3,self.nVibs))  
        #aMu2d_{i,j}=<1_i,1_j|dipoleMoment|0>=Sum_n^walkers(Dipole_n*q_i*q_j*dw_n)/Sum(dw) #Equation 4 for \psi_n=\psi_{1,1}
        aMu2d=np.zeros((3,self.nVibs,self.nVibs))

        for P,R,d in zip(dipoleMoments,q,dw):
            aMu=aMu+np.outer(P,-R)*d  #off by a sign change for some reason so I changed it to be -R instead of R...I guess P could be backwards? 
            temp=np.outer(R,R)
            temp=temp.reshape((self.nVibs*self.nVibs,1))
            aMu2d=aMu2d+(np.outer(P,temp).reshape((3,self.nVibs,self.nVibs)))*d

        aMu=aMu/np.sum(dw)
        aMu2d=aMu2d/np.sum(dw)

        np.savetxt("WtvsMu-qs-a.data",zip(dw,np.sum(dipoleMoments**2,axis=1),q[:,5],q[:,6]))

        aMu=aMu.transpose()
        aMu2d=aMu2d.transpose()  #switches around the axes

        magAvgMu=np.zeros((self.nVibs))
        #Finding |<1_i|dipoleMoment|0>|^2 
        for m,mode in enumerate(aMu):
            magAvgMu[m]=magAvgMu[m]+(mode[0]**2+mode[1]**2+mode[2]**2)/(q2ave[m])
        ivibmode=0
        
        magMu2d=np.zeros((self.nVibs,self.nVibs))
        #Finding |<1_i,1_j|dipoleMoment|0>|^2
        for ivibmode in range(self.nVibs):
            aMu2d[ivibmode,ivibmode]=2*alpha[ivibmode]*aMu2d[ivibmode,ivibmode]-averageDipoleMoment_0 
            magMu2d[ivibmode,ivibmode]=np.sum(aMu2d[ivibmode,ivibmode]**2)
            magMu2d[ivibmode,ivibmode]=magMu2d[ivibmode,ivibmode]/(q2ave2d[ivibmode,ivibmode]*4.0*alpha[ivibmode]**2-4.0*alpha[ivibmode]*q2ave[ivibmode]+1.0)
            for jvibmode in range(ivibmode):
                magMu2d[ivibmode,jvibmode]=np.sum(aMu2d[ivibmode,jvibmode]**2)
                magMu2d[ivibmode,jvibmode]=magMu2d[ivibmode,jvibmode]/q2ave2d[ivibmode,jvibmode]
                magMu2d[jvibmode,ivibmode]=magMu2d[ivibmode,jvibmode]

        print "m alpha   alpha'=0.5/<q2>        Vq            Vq/q2ave      Tq            Eq          |<1|dipole|0>|"
        for i in range(9):
            print i, alpha[i], alphaPrime[i],au2wn*Vq[i], au2wn*Vq[i], Vq[i]/q2ave[i]*au2wn, Eq[i]*au2wn, magAvgMu[i]
            
        if verbose: print 'energies!', zip(Vq2d[0]*au2wn,Tq2d[0]*au2wn,Eq2d[0]*au2wn)

        print 'Spectrum Info'
        fundamentalFileName=self.path+'fundamentals.data'
        comboFileName=self.path+'combinationOverrtoneBands.data'
        fundamentalFile=open(fundamentalFileName,'w')
        comboFile=open(comboFileName,'w')

        for i in range(self.nVibs):
            print Eq[i]*au2wn, magAvgMu[i], 'mode ',i 
            fundamentalFile.write(str( Eq[i]*au2wn)+"   "+str(magAvgMu[i])+"       "+str(i)+"\n")
        for i in range(self.nVibs):
            for j in range(i):
                print Eq2d[i,j]*au2wn , magMu2d[i,j],'combination bands' , i,j
                comboFile.write(str( Eq2d[i,j]*au2wn)+"   "+str(magMu2d[i,j])+"       "+str(i)+" "+str(j)+"\n")
        for i in range(self.nVibs):
            print Eq2d[i,i]*au2wn , magMu2d[i,i],'overtone bands' , i,i
            comboFile.write(str( Eq2d[i,i]*au2wn)+"   "+str(magMu2d[i,i])+"       "+str(i)+" "+str(i)+"\n")
        fundamentalFile.close()
        comboFile.close
        return Eq*au2wn,magAvgMu, Eq2d*au2wn,magMu2d


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
        alpha=1.0*eigval #AS in alpha_j in equation 5 of the JPC A 2011 h5o2 dier paper                                          
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
