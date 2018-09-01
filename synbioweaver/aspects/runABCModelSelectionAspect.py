from synbioweaver.core import *
import os

from abcsysbio import model_cu
from abcsysbio import data
from abcsysbio import abcsmc
from abcsysbio import input_output
from abcsysbio import kernels
from abcsysbio import euclidian

# This currently fixes a lot of options
# Assumes: ODE integration on CUDA, nbatch = 25000
#          distance file is in distance.py
#          nparticles = 128

class RunABCModelSelection(Aspect):

    def __init__(self, names, priors, inits, times, fit, feps, distance):
        super(RunABCModelSelection,self).__init__()
        self.names = names
        self.priors = priors
        self.inits = inits
        self.times = times
        self.fit = fit
        self.feps = feps
        self.distance = distance

    def mainAspect(self):
        self.addWeaverOutput(self.runABCModelSelection)
        
    def runABCModelSelection(self, weaverOutput):
        
        # We are expecting a context only
        if getattr(weaverOutput, "getContext", None) != None:
            self.models = weaverOutput.getContext()
        else:
            print " runABCModelSelection : getContext() is not available. Quitting"
            exit()

        self.nmodels = len(self.models)

        # We are expecting that writeMultipleCudaFileODE has also been run
        if getattr(weaverOutput, "writeMultipleCudaFileODE", None) != None:
            self.cudaStrs = weaverOutput.writeMultipleCudaFileODE()
        else:
            print " runABCModelSelection : writeMultipleCudaFileODE is not available. Quitting"
            exit()

        for i,m in enumerate(self.cudaStrs):
            cudaFile = open(self.names[i]+".cu","w")
            print >>cudaFile, m
            cudaFile.close()

        nparticles = 128

        # create new models
        abcmodels = []
        for i in range(self.nmodels):
            print "new abc model:", self.names[i], self.priors[i]
            new_model = model_cu.cuda_model( name=self.names[i], nspecies=len(self.inits[i]), nparameters=len(self.priors[i]), 
                                             prior=self.priors[i], x0prior=self.inits[i],
                                             source=self.names[i], integration="ODE", fit=self.fit[i],
                                             dt=-1.0, beta=1, timepoints=self.times, logp=False)
            abcmodels.append( new_model )

        # create a data model
        # since this is running in design mode the data is all in the distance file
        data_new = data.data( self.times, [] )
        
        # IO
        #io = input_output.input_output(folder="res_abc", restart=False, diagnostic=True, plotDataSeries=True, havedata=False )
        io = input_output.input_output(folder="res_abc", restart=False, diagnostic=False, plotDataSeries=False, havedata=False ) 

        # creates separate folders for each model (here we are using only one)
        io.create_output_folders(self.names, numOutput=nparticles, pickling=True, simulation=False)
        
        # setting up ABC SMC
        debug = 1
        kernelfn = kernels.getKernel
        kernelpdffn = kernels.getPdfParameterKernel
        perturbfn = kernels.perturbParticle

        customABC = __import__(self.distance)
        distancefn = customABC.distance

        modelprior = [1/float(self.nmodels) for i in range(self.nmodels) ]

        algorithm = abcsmc.abcsmc(models=abcmodels, nparticles=nparticles, modelprior=modelprior, data=data_new, 
                                  beta = 1, nbatch=25000,
                                  modelKernel=0.7, debug=debug, timing=True, distancefn=distancefn, 
                                  kernel_type=1, kernelfn=kernelfn, kernelpdffn=kernelpdffn, perturbfn=perturbfn )

        # go
        algorithm.run_automated_schedule(final_epsilon=self.feps, alpha=0.3, io=io)
        

        print "##### RunABCModelSelection Success"
