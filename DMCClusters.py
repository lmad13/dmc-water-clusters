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
        self.alpha=0.25/self.dtau

        self.recrossing=True if 'node' in potentialName else False

    def set_dtau(self):
        self.dtau=5.0
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

    def propagate(self,x,nSteps,setV_ref=False,ConstantV_ref=0,printCensus=False,initialPop=0):
        print 'ready to propagate for',nSteps, 'steps on x (shaped:',x.shape,')'
        if initialPop==0:#if it is the default value of zero...it needs to be set.
            initialPop=x.shape[0]
        if setV_ref:
            v_ref=ConstantV_ref
        else:
            v_ref=np.average(self.molecule.V(x))#+(self.alpha*(1-float(N_size_step)/float(nSize)))

        self.currentPop=initialPop
        population=[]
        population.append(self.currentPop)
        print 'and x.size is', initialPop

        vRefList=[]
        vRefList.append(v_ref)

        descendants=np.zeros((initialPop))
        whoYaFrom=np.arange(initialPop)
        
        for step in range(nSteps):
            print 'step #',step
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
            survivors=x[mask_survive]
            if printCensus: print 'Census: Deaths:',nDeaths,

            # Recrossing goes here
            if self.recrossing:
                #m2=2*(massO+2*massH)*massConversionFactor
                #m1=massH*massConversionFactor
                newReducedMass=self.molecule.calcReducedmass(x)
                oldDist=self.molecule.calcRn(x-dx)
                newDist=self.molecule.calcRn(x)
                oldReducedmass=calcReducedmass(x-dx)
                P_recrossDeath=np.exp(-2.0*(oldDist)*newDist*np.sqrt(oldReducedMass*newReducedMass)/self.dtau) ##mass is reduced mass!
                #caution, this is clever, you are combining the probability of crossing and the probability of recrossing
                if self.plotting:
                    plt.scatter(x,P_recrossDeath)
                    plt.show()
                Diff=N_r-P_recrossDeath
                mask_survive_recross=(Diff>0)
                tempRecrossCensus=np.sum(np.array(Diff<0).astype(int))
                mask_survive=np.logical_and(mask_survive, mask_survive_recross)
                survivors=x[mask_survive]            

                if printCensus: print 'Deaths by recrossing: ',tempRecrossCensus
            
            #Creation of a random selection of walkers in the classically allowed region
            P_exp_b=np.exp(-(v-v_ref)*self.dtau)-1.0
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
                    if weight>10:
                        #this really shouldn't happen                                    
                        print 'weight is too big, resetting to 10'
                        print x[n],v(n),'<',v_ref, -(v(n)-v_ref)
                        weight=10
                    addBirthtot=addBirthtot+weight
                    temp=np.tile(particle,weight)
                    temp_whoYaFrom=np.tile(whoYaFrom[p],weight)
                    new_pop=np.concatenate((new_pop,temp))
                    new_pop_whoYaFrom=np.concatenate((new_pop_whoYaFrom,temp_whoYaFrom))
            next_gen=new_pop
            next_gen_whoYaFrom=new_pop_whoYaFrom
            #collect survivors and next generation                                       
            new_population=np.concatenate((survivors,next_gen))
            self.currentPop=x.shape[0]
            if printCensus: print '. Births:',nBirths, ". Add' births: ", addBirthtot

            #readjust V_ref
            v_average=np.average(self.molecule.V(new_population))
            if not setV_ref:
                v_ref=v_average+(self.alpha*(1-float(self.currentPop)/float(initialPop)))
            if printCensus: print '(',N_size_step,'/',nSize,') v_ref',v_ref, '=', v_average,'+', (self.alpha*(1-float(N_size_step)/float(nSize)))
            vRefList.append(v_ref)
            population.append(self.currentPop)

            whoYaFrom=np.concatenate((whoYaFrom[mask_survive],next_gen_whoYaFrom))
            x=new_population

        for anc in whoYaFrom:
            descendants[anc]=descendants[anc]+1

        return  vRefList, population, x, descendants
    
#wavefunction
#initialize

#propagate

#calculate V

