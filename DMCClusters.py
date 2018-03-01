import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import molecularInfo

au2wn=219474.63
au2ang=0.529177249
massConversionFactor=1.000000000000000000/6.02213670000e23/9.10938970000e-28#1822.88839 


class wavefunction:

    def __init__(self,moleculeName, nWalkers, potentialName='Ground State'):
        self.molecule=molecularInfo.molecule(moleculeName)

        self.pathDict=self.loadDict("paths.dict")
        print 'initialized', nWalkers,' coordinates for ', self.molecule.name
        self.nAtoms=self.molecule.nAtoms
        self.nDim=3
        

        #based on the molecule find the coordinates
        self.x=np.zeros((nWalkers,self.nAtoms,self.nDim))
        self.x[:]=self.molecule.getInitialCoordinates()
        self.D=.5
        self.set_dtau()
        self.mass=self.molecule.get_mass()
        self.sigma_dx=(2.0000*self.D*self.dtau/self.mass)**0.5
        self.alpha=.25/self.dtau

        self.recrossing=False

    def setNodalSurface(self,surfaceName,side='Both'):
        self.surfaceName=surfaceName
        self.side=side # Right: x>node
        self.molecule.setNodalSurface(self.surfaceName,self.side)
        self.recrossing=True
    def set_dtau(self):
        self.dtau=10.0
        print 'dtau is ', self.dtau, 'imaginary atomic time units'

    def loadDict(self,fileName):
        fileIn=open(fileName,'r')
        pathDict={}
        for line in fileIn:
            [key,element]=line.split()
            pathDict[key]=element
        return pathDict

    def diffuse(self):
        dx=np.zeros((self.currentPop,self.nAtoms,self.nDim))
        for atom in range(self.nAtoms):
            centerOfGaussian=0.000000
            dx[:,atom:atom+1,:]=np.random.normal(centerOfGaussian,self.sigma_dx[atom],(self.currentPop,1,self.nDim))

        return dx

    def exportCoords(self,x,fileName,nDesc):
        fileout=open(fileName,'w')
        for particle, note in zip(x,nDesc):
            fileout.write(str(self.nAtoms)+' \n'+str(note)+' \n')
            for atomName, atom in zip(self.molecule.names,particle):
                fileout.write(str(atomName)+'   '+str(au2ang*atom[0])+"   "+str(au2ang*atom[1])+"   "+str(au2ang*atom[2])+"\n")
            fileout.write("\n")
        fileout.close()
        return
    def loadCoords(self,fileName):
        filein=open(fileName,'r')
        fdata=filein.readlines()
        nMolecules=len(filedata)/(self.naAtoms+3)
        repeatUnit=self.nAtoms+3
        coord=np.zeros((nMolecules,nAtoms,3))
        descCoord=np.zeros((nMolecules))
        for ln in range(nMolecules):
            coordstemp=[]
            for an in range(self.nAtoms):
                coordstemp.append(float(filedata[ln*repeatUnit+2+an].split()[1]))
                coordstemp.append(float(filedata[ln*repeatUnit+2+an].split()[2]))
                coordstemp.append(float(filedata[ln*repeatUnit+2+an].split()[3]))
            coord[ln,:,:]=np.reshape(coordstemp,(self.nAtoms,3))
            descCoord=float(filedata[ln*repeatUnit+1].split()[0])
        coord=coord/au2ang
        return coord,descCoord

    def propagate(self,x,nSteps,setV_ref=False,ConstantV_ref=0,printCensus=True,initialPop=0):
        #print 'ready to propagate for',nSteps, 'steps on x (shaped:',x.shape,') '

        if initialPop==0:#if it is the default value of zero...it needs to be set.
            initialPop=x.shape[0]
        self.currentPop=x.shape[0]
        if setV_ref:
            v_ref=ConstantV_ref
        else:
            v_ref=np.average(self.molecule.V(x))+(self.alpha*(1-float(self.currentPop)/float(initialPop)))

        population=[]
        population.append(self.currentPop)
        
        vRefList=[]
        vRefList.append(v_ref)

        descendants=np.zeros((self.currentPop))
        whoYaFrom=np.arange(self.currentPop)
        
        printRate=100

        for step in range(nSteps):
            if step%printRate==0 and printCensus: print 'step #',step, 
            dx=self.diffuse()
            x=x+dx
            v=self.molecule.V(x)

            
            #Elimination of a random selection of walkers 
            #in the classically forbidden region
            N_r=np.random.random(self.currentPop)
            P_d=1-np.exp(-(v-v_ref)*self.dtau)
            Diff=N_r-P_d
            mask_survive=(Diff>0)
            nDeaths=np.sum(np.array(Diff<0).astype(int))
            #survivors=x[mask_survive]
            if step%printRate==0 and printCensus: print 'Census: Current Population:',self.currentPop,'Deaths:',nDeaths,

            #LASSIE MODE
            #DEATH BY PEF HOLES
            mask_not_holes=(v>0.0) #mask
            mask_survive=np.logical_and(mask_survive,mask_not_holes)

            # removal of all walkers that cross and a random selection of walkers that are too close to the node
            if self.recrossing:
                # m2=2*(massO+2*massH)*massConversionFactor
                # m1=massH*massConversionFactor
                newReducedMass=self.molecule.calcReducedmass(x)
                oldDist=self.molecule.calcRn(x-dx)
                newDist=self.molecule.calcRn(x)
                crossed=(oldDist*newDist<0)

                oldReducedMass=self.molecule.calcReducedmass(x-dx)
                P_recrossDeath=np.exp(-2.0*(oldDist)*newDist*np.sqrt(oldReducedMass*newReducedMass)/self.dtau) ##mass is reduced mass!
                #caution, this is clever, you are combining the probability of crossing and the probability of recrossing
                #if self.plotting:
                #    plt.scatter(x,P_recrossDeath)
                #    plt.show()
                
                Diff=N_r-P_recrossDeath
                mask_survive_recross=(Diff>0)
                tempRecrossCensus=np.sum(np.array(Diff<0).astype(int))
                mask_survive=np.logical_and(mask_survive, mask_survive_recross)
                if step%printRate==0 and  printCensus: print 'Deaths by recrossing: ',tempRecrossCensus
            
            
            survivors=x[mask_survive]            

#Creation of a random selection of walkers in the classically allowed region and who did not die by recrossing, so the survivors
        
                
    
            P_exp_b=np.exp(-(v-v_ref)*self.dtau)-1.0
            if self.recrossing:
                P_exp_b[crossed]=0.00000  #if it crossed, the probability that it gives birth should be zero

            weight_P_b=P_exp_b.astype(int)
            P_b=P_exp_b-weight_P_b            #classically allowed region                            
            # P_b[np.logical_not(mask_survive)]=0.0                                       
            
            Diff=N_r-P_b
            mask_b = (Diff<0)
            next_gen=x[mask_b]
            new_pop_whoYaFrom=whoYaFrom[mask_b]
            nBirths=np.sum(np.array(Diff<0).astype(int))
            addBirthtot=0
            new_pop=next_gen
            
           #for the additional births                                                   
            for n,(particle,weight) in enumerate(zip(x,weight_P_b)):
                if weight>0: #i.e. the dead can't reproduce                              
                    #print 'weight for tiling:', weight, n
                    if weight>10:
                        #this really shouldn't happen                                    
                        print 'weight is too big, resetting to 10'
                        print v[n]*au2wn,'<',v_ref*au2wn, weight, '\n',x[n]
                        weight=10
                        end
                    addBirthtot=addBirthtot+weight
                    temp=np.array([particle])
                    temp_whoYaFrom=np.array([whoYaFrom[n]])
                    for i in range(weight-1):
                        temp=np.concatenate((temp,np.array([particle])))
                        temp_whoYaFrom=np.concatenate((temp_whoYaFrom,np.array([whoYaFrom[n]])))
#                    temp=np.tile(particle,weight)
#                    temp=np.reshape(temp,(weight,self.nAtoms,self.nDim))
#                    temp_whoYaFrom=np.tile(whoYaFrom[n],weight)

                    new_pop=np.concatenate((new_pop,temp))
                    new_pop_whoYaFrom=np.concatenate((new_pop_whoYaFrom,temp_whoYaFrom))

            next_gen=new_pop
            next_gen_whoYaFrom=new_pop_whoYaFrom
            #collect survivors and next generation                                       
            new_population=np.concatenate((survivors,next_gen))
            whoYaFrom=np.concatenate((whoYaFrom[mask_survive],next_gen_whoYaFrom))
            #print 'shapes for next cycle:',new_population.shape,whoYaFrom.shape
            x=new_population
            
            #let's make sure none of the walkers crossed to the other side...
            position_on_coord=self.molecule.calcRn(x)
            #print 'wrong side!', position_on_coord[(position_on_coord>0)],
            #print 'from these array locations!',np.where(position_on_coord>0)
#            if np.where(position_on_coord>0)[0].shape[0]>1:
#               end
            
            self.currentPop=x.shape[0]
            if step%printRate==0 and printCensus: print 'Births:',nBirths, "Extra Births: ", addBirthtot,

            #readjust V_ref
            v_average=np.average(self.molecule.V(new_population))
            
            if not setV_ref:
                v_ref=v_average+(self.alpha*(1-float(self.currentPop)/float(initialPop)))
            if step%printRate==0 and printCensus: print 'v_ref',v_ref, '=', v_average,'+', (self.alpha*(1-float(self.currentPop)/float(initialPop)))
            vRefList.append(v_ref)
            population.append(self.currentPop)
            if self.currentPop<5:
                print 'massive walker die off!', end

        for anc in whoYaFrom:
            descendants[anc]=descendants[anc]+1

        return  vRefList, population, x, descendants
    
#wavefunction
#initialize

#propagate

#calculate V
