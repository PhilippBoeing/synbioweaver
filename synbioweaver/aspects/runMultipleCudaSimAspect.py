from synbioweaver.core import *
import os, sys

from numpy import *
from numpy.random import *
import matplotlib.pyplot as plt

import cudasim
import cudasim.EulerMaruyama as EulerMaruyama
import cudasim.Gillespie as Gillespie
import cudasim.Lsoda as Lsoda


class RunMultipleCudaSim(Aspect):
    def __init__(self, names, vpars, inits, times, obs):
        super(RunMultipleCudaSim,self).__init__()
        self.names = names
        self.vpars = vpars
        self.inits = inits
        self.times = times
        self.obs = obs
     

    def mainAspect(self):
        self.addWeaverOutput(self.runMultipleCudaSim)
        
    def runMultipleCudaSim(self, weaverOutput):
        
        # We are expecting either a set of reactions or a context
        if getattr(weaverOutput, "getContext", None) != None:
            self.models = weaverOutput.getContext()
        else:
            print " runMultipleCudaSim : getContext() not available. Quitting"
            exit()

        # We are expecting that writeMultipleCudaFileODE has also been run
        if getattr(weaverOutput, "writeMultipleCudaFileODE", None) != None:
            self.cudaStrs = weaverOutput.writeMultipleCudaFileODE()
        else:
            print " runMultipleCudaSim : writeMultipleCudaFileODE is not available. Quitting"
            exit()

        for i,m in enumerate(self.cudaStrs):
            cudaFile = open(self.names[i]+".cu","w")
            print >>cudaFile, m
            cudaFile.close()

        for i,model in enumerate(self.models):
            self.runModel(model,self.names[i],self.vpars[i],self.inits[i],self.obs[i])

    def runModel(self,model,name,vpars,inits,obs):
        npars = len(model.freeparams)
        nspecies = model.nspecies
        print "running model cudasim:", npars, nspecies, obs
        nsims = 500

        parameters = zeros( [nsims,npars+1] )
        species = zeros([nsims,nspecies])

        for i in range(nsims):
            parameters[i,0] = 1
            for j in range(npars):
                parameters[i,j+1] = uniform(vpars[j][0], vpars[j][1])
            for j in range(nspecies):
                species[i,j] = uniform(inits[j][0], inits[j][1])

        modelInstance = Lsoda.Lsoda(self.times, name+".cu" , dt=-1)
        result = modelInstance.run(parameters, species)
        sig = result[:,0,:,:] # beta set to 1
        #print "done running model:" shape(sig)
        
        nobs = len(obs)
        #print "sum over observables:", nobs
        sigO = zeros( (shape(sig)[0], shape(sig)[1], nobs) )
        for j in range(nobs):
            for k in obs[j]:
                #print "\tsumming:", j, k
                sigO[:,:,j] += sig[:,:,k]
        #print shape(sigO)

        # shape should be q x nt x nobs
        res = percentile( sigO, [5,50,95], axis=0 )
        
        for k in range(nobs):
            lines = []
            # plot median
            lines.append(plt.plot(self.times, res[1,:,k], label="o"+str(k+1) )  )
            plt.setp(lines[0], linewidth=2)
            # plot 95 CR
            plt.fill_between(self.times,res[0,:,k],res[2,:,k], alpha=0.5)
            plt.savefig("plot-res-"+name+"-o"+str(k+1)+".pdf",bbox_inches='tight')
            plt.close()

        out = open("data-res-"+name+".txt",'w')
        print >>out, "time", 
        for k in range(nobs):
                print >>out, 'obs'+str(k)+'-5', 'obs'+str(k)+'-50', 'obs'+str(k)+'-95',
        print >>out,  ""

        for j in range(len(self.times)):
            print >>out, self.times[j], 
            for k in range(nobs):
                print >>out, res[0,j,k], res[1,j,k], res[2,j,k],
            print >>out,  ""
        out.close()

        
        #print "##### RunCudaSim Success"
