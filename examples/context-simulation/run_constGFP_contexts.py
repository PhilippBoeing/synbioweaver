from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.lagLogisticGrowthAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *
from synbioweaver.aspects.writeSBMLModelAspect import *
from synbioweaver.aspects.runCudaSim import *
from synbioweaver.aspects.writeCudaFileODE import *

declareNewMolecule('GFP')

### parameters for cuda-sim: 
###    uniform prior ranges over the parameters
###    uniform prior ranges over the initial conditions
###    observables of the system: for example in model 2 we want to sum over N and Nd

# define time structure
times = arange(0,480,1.0)
# Model 1
m1_prior = [[0.05,0.1],[0.008,0.01]]
m1_inits = [[0,0.01]]
m1_obs = [[0]]

# Model 2
m2_prior = [[0.05,0.1],[0.008,0.01],[0.03,0.012],[1.24,1.35],[0.01,0.011],[0.005,0.0051]]
m2_inits = [[0,0.01],[0,0.01],[0.5,0.6]]
m2_obs = [[0],[1,2]]

class constGFP(Circuit):
    def mainCircuit(self):
        # constitutive promoter
        self.addPart(ConstitutivePromoter)
        self.addPart(CodingRegion(GFP))

print "####### circuit ####### ####### ####### "
compiledDesign1 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork, WriteCudaFileODE, WriteSBMLModel, 
                         RunCudaSim(name="m1",vpars=m1_prior, inits=m1_inits, times=times, obs=m1_obs ) ).output()
compiledDesign1.printReactionNetwork()
compiledDesign1.runCudaSim()

# If you want the cuda code or sbml file writing out
#cudaStr = compiledDesign1.writeCudaFileODE()
#cudaFile = open("constGFP_ODE.cu","w")
#print >>cudaFile, cudaStr
#cudaFile.close()

#modelStr = compiledDesign1.writeSBMLModel()
#sbmlFile = open("constGFP.sbml","w")
#print >>sbmlFile, modelStr
#sbmlFile.close()

print "####### circuit + lag logistic growth ####### ####### #######"
compiledDesign2 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, LagLogisticGrowth, PrintReactionNetwork, WriteCudaFileODE, WriteSBMLModel,
                          RunCudaSim(name="m2",vpars=m2_prior, inits=m2_inits, times=times, obs=m2_obs ) ).output()
compiledDesign2.printReactionNetwork()
compiledDesign2.runCudaSim()

# If you want the cuda code or sbml file writing out
#cudaStr = compiledDesign2.writeCudaFileODE()
#cudaFile = open("constGFP_lag_ODE.cu","w")
#print >>cudaFile, cudaStr
#cudaFile.close()

#modelStr = compiledDesign2.writeSBMLModel()
#sbmlFile = open("constGFP_lag.sbml","w")
#print >>sbmlFile, modelStr
#sbmlFile.close()




