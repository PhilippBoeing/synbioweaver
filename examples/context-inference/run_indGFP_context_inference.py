from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.MassActionKineticsRNAAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.lagLogisticGrowthAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *
from synbioweaver.aspects.writeSBMLModelAspect import *

from writeMultipleCudaFileODEAspect import *
from printMultipleModelsAspect import *
from runMultipleCudaSimAspect import *
from runABCModelSelectionAspect import *
from contextModelGeneration import *

from numpy import *

declareNewMolecule('GFP')

# data timepoints
times = array([0, 10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 350, 400, 450]) 

class singleSystem(Circuit):
    def mainCircuit(self):
        self.addPart(ConstitutivePromoter)
        self.addPart(RBS)
        self.addPart(CodingRegion(GFP))
        self.addPart(Terminator)


if 1:
    # Use context generation to generate and simulate a model where translation rate decreases in the second context

    # p1, p2, p3, p4: GFP transcription rate, GFP translation, mRNA deg, protein deg
    sim_prior_all  = [ [4.25,4.3], [5,5.1], [0.13,0.14], [0.03,0.035], [4.25,4.3], [1,1.1], [0.13,0.14], [0.03,0.035] ]
    sim_prior_free = [ [4.25,4.3], [5,5.1], [0.13,0.14], [0.03,0.035] ]
    sim_inits = [ [0,1], [0,0.1], [0,1], [0,0.1] ]
    sim_obs = [[1],[3]]

    compiledDesign = Weaver(singleSystem, PromoterMapping, MassActionKineticsRNA, ContextGenerator(context_models=["all","free"]), PrintMultipleModels, WriteMultipleCudaFileODE,
                             RunMultipleCudaSim(names=["model-sim-dep","model-sim-free"],vpars=[sim_prior_all,sim_prior_free], inits=[sim_inits, sim_inits], times=times, obs=[sim_obs, sim_obs] ) 
                      ).output()
    compiledDesign.printMultipleModels()
    compiledDesign.runMultipleCudaSim()

if 1:
    # Now perform Bayesian model selection to formally test contextual dependence
    # M0: no context effect == both submodels have the same parameters
    # M1: context dependence == submodels have different parameters

    # 0 means constant / 2 means uniform prior
    prior_c_free = [ [0, 1.0, 0], [2, 4,4.5], [2, 1,5.5], [2, 0.1,0.15], [2, 0.03,0.04] ]
    prior_c_all = [ [0, 1.0, 0], [2, 4,4.5], [2, 1,5.5], [2, 0.1,0.15], [2, 0.03,0.04], [2, 4,4.5], [2, 1,5.5], [2, 0.1,0.15], [2,0.03,0.04] ]

    init = [[0, 0.0, 0], [0, 0.0, 0], [0, 0.0, 0.0], [0, 0.0, 0.0]]
    fit = ['samplePoints[:,1]', 'samplePoints[:,3]']
    feps = [5.0]

    compiledDesign = Weaver(singleSystem, PromoterMapping, MassActionKineticsRNA, ContextGenerator(context_models=["free","all"]), PrintMultipleModels, WriteMultipleCudaFileODE,
                            RunABCModelSelection(names=["model-context-free","model-context-dep"], priors=[prior_c_free, prior_c_all], 
                                                 inits=[init, init], fit=[fit, fit], feps=feps, times = times, distance="distance_dep")
                        ).output()

    compiledDesign.printMultipleModels()
    compiledDesign.runABCModelSelection()


