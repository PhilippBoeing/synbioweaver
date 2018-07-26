from synbioweaver.core import *
import os, sys

from numpy import *
from numpy.random import *
import matplotlib.pyplot as plt

import cudasim
import cudasim.EulerMaruyama as EulerMaruyama
import cudasim.Gillespie as Gillespie
import cudasim.Lsoda as Lsoda


class RunCudaSim(Aspect):
    def __init__(self, name, vpars, inits, times, obs):
        super(RunCudaSim,self).__init__()
        self.name = name
        self.vpars = vpars
        self.inits = inits
        self.times = times
        self.obs = obs
     

    def mainAspect(self):
        self.addWeaverOutput(self.runCudaSim)
        
    def runCudaSim(self, weaverOutput):
        
        # We are expecting either a set of reactions or a context
        if getattr(weaverOutput, "getContext", None) != None:
            self.model = weaverOutput.getContext()
        else:
            if getattr(weaverOutput, "getReactions", None) != None:
                        self.model = weaverOutput.getReactions()
            else:
                print " RunCudaSim : Neither getContext() or getReactions() is available. Quitting"
                exit()

        # We are expecting that writeCudaFileODE has also been run
        if getattr(weaverOutput, "writeCudaFileODE", None) != None:
            self.cuda_code = weaverOutput.writeCudaFileODE()
        else:
            print " RunCudaSim : writeCudaFileODE is not available. Quitting"
            exit()

        cudaFile = open(self.name+".cu","w")
        print >>cudaFile, self.cuda_code
        cudaFile.close()
        
        npars = len(self.model.parameters)
        nspecies = self.model.nspecies
        nsims = 500

        parameters = zeros( [nsims,npars+1] )
        species = zeros([nsims,nspecies])

        for i in range(nsims):
            parameters[i,0] = 1
            for j in range(npars):
                parameters[i,j+1] = uniform(self.vpars[j][0], self.vpars[j][1])
            for j in range(nspecies):
                species[i,j] = uniform(self.inits[j][0], self.inits[j][1])

        modelInstance = Lsoda.Lsoda(self.times, self.name+".cu" , dt=-1)
        result = modelInstance.run(parameters, species)
        sig = result[:,0,:,:] # beta set to 1
        #print "done running model:" shape(sig)
        
        nobs = len(self.obs)
        #print "sum over observables:", nobs
        sigO = zeros( (shape(sig)[0], shape(sig)[1], nobs) )
        for j in range(nobs):
            for k in self.obs[j]:
                #print "\tsumming:", j, k
                sigO[:,:,j] += sig[:,:,k]
        #print shape(sigO)

        # shape should be q x nt x nobs
        res = percentile( sigO, [5,50,95], axis=0 )
        
        for k in range(nobs):
            lines = []
            lines.append(plt.plot(self.times, res[1,:,k], label="o"+str(k+1) )  )
            plt.setp(lines[0], linewidth=2)
            plt.fill_between(self.times,res[0,:,k],res[2,:,k], alpha=0.5)
            plt.savefig("plot-res-"+self.name+"-o"+str(k+1)+".pdf",bbox_inches='tight')
            plt.close()

        out = open("data-res-"+self.name+".txt",'w')
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
