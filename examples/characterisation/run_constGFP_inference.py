from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.lagLogisticGrowthAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *
from synbioweaver.aspects.writeSBMLModelAspect import *
from synbioweaver.aspects.runCudaSim import *
from synbioweaver.aspects.writeCudaFileODE import *
from synbioweaver.aspects.runABCAspect import *

declareNewMolecule('GFP')

# data timepoints
times = array([0, 10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 350, 400, 450]) 

class constGFP(Circuit):
    def mainCircuit(self):
        # constitutive promoter
        self.addPart(ConstitutivePromoter)
        self.addPart(CodingRegion(GFP))

if 1:
    ### simulate some data using cuda-sim: 
    ###    uniform prior ranges over the parameters
    ###    uniform prior ranges over the initial conditions
    ###    observables of the system: for example in model 2 we want to sum over N and Nd

    sim_prior = [[0.05,0.051],[0.008,0.0081],[0.03,0.031],[1.20,1.21],[0.01,0.011],[0.005,0.0051]]
    sim_inits = [[0,0.01],[0,0.01],[0.5,0.51]]
    sim_obs = [[0],[1,2]]

    ### circuit + lag logistic growth
    compiledDesign1 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, LagLogisticGrowth, PrintReactionNetwork,  WriteSBMLModel, WriteCudaFileODE,
                          RunCudaSim(name="growth",vpars=sim_prior, inits=sim_inits, times=times, obs=sim_obs ) ).output()
    compiledDesign1.printReactionNetwork()
    compiledDesign1.runCudaSim()

if 1:
    # 0 means constant / 2 means uniform prior
    inf_prior  = [[0, 1.0, 0], [2, 0.05, 0.1], [2, 0.008, 0.01], [2, 0.03, 0.035], [2, 1.24, 1.35], [2, 0.01, 0.015], [2, 0.005, 0.006]]
    inf_inits = [[0, 0.0, 0], [0, 0.0, 0], [2, 0.5, 0.6]]
    inf_fit = ['samplePoints[:,0]', 'samplePoints[:,1]+samplePoints[:,2]']
    feps = [1.0]

    ### run the same model through ABC inference
    compiledDesign2 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, LagLogisticGrowth, PrintReactionNetwork, WriteSBMLModel,  WriteCudaFileODE,
                             RunABC(name="growth", priors=inf_prior, inits=inf_inits, times=times, fit=inf_fit, feps=feps) ).output()
    compiledDesign2.runABC()





