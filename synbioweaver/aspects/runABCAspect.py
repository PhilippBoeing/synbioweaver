from synbioweaver.core import *
import os

from abcsysbio import model_cu
from abcsysbio import data
from abcsysbio import abcsmc
from abcsysbio import input_output
from abcsysbio import kernels
from abcsysbio import euclidian

# This currently fixes a lot of options
# Assumes: single model
#          ODE integration on CUDA, nbatch = 25000
#          distance file is in distance.py
#          nparticles = 128

class RunABC(Aspect):

    def __init__(self, name, priors, inits, times, fit, feps):
        super(RunABC,self).__init__()
        self.name = name
        self.priors = priors
        self.inits = inits
        self.times = times
        self.fit = fit
        self.feps = feps

    def mainAspect(self):
        self.addWeaverOutput(self.runABC)
        
    def runABC(self, weaverOutput):
        
        # We are expecting either a set of reactions or a context
        if getattr(weaverOutput, "getContext", None) != None:
            self.model = weaverOutput.getContext()
        else:
            if getattr(weaverOutput, "getReactions", None) != None:
                        self.model = weaverOutput.getReactions()
            else:
                print " RunABC : Neither getContext() or getReactions() is available. Quitting"
                exit()

        # We are expecting that writeCudaFileODE has also been run
        if getattr(weaverOutput, "writeCudaFileODE", None) != None:
            self.cuda_code = weaverOutput.writeCudaFileODE()
        else:
            print " RunABC : writeCudaFileODE is not available. Quitting"
            exit()

        nparticles = 128

        cudaFile = open(self.name+".cu","w")
        print >>cudaFile, self.cuda_code
        cudaFile.close()

        # create a new model
        models = []
        new_model = model_cu.cuda_model( name=self.name, nspecies=len(self.inits), nparameters=len(self.priors), 
                                          prior=self.priors, x0prior=self.inits,
                                          source=self.name, integration="ODE", fit=self.fit,
                                          dt=-1.0, beta=1, timepoints=self.times, logp=False)
        models.append( new_model )

        # create a data model
        # since this is running in design mode the data is all in the distance file
        data_new = data.data( self.times, [] )
        
        # IO
        io = input_output.input_output(folder="res_abc", restart=False, diagnostic=True, plotDataSeries=True, havedata=False )
        # creates separate folders for each model (here we are using only one)
        io.create_output_folders([self.name], numOutput=nparticles, pickling=True, simulation=False)
        
        # setting up ABC SMC
        debug = 1
        kernelfn = kernels.getKernel
        kernelpdffn = kernels.getPdfParameterKernel
        perturbfn = kernels.perturbParticle

        customABC = __import__("distance")
        distancefn = customABC.distance

        algorithm = abcsmc.abcsmc(models=models, nparticles=nparticles, modelprior=[1.0], data=data_new, 
                                  beta = 1, nbatch=25000,
                                  modelKernel=0.7, debug=debug, timing=True, distancefn=distancefn, 
                                  kernel_type=1, kernelfn=kernelfn, kernelpdffn=kernelpdffn, perturbfn=perturbfn )

        # go
        algorithm.run_automated_schedule(final_epsilon=self.feps, alpha=0.3, io=io)
        

        print "##### RunABC Success"
