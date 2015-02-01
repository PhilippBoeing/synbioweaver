from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.MassActionKineticsRNAAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.expGrowthAspect import *
from synbioweaver.aspects.logisticGrowthAspect import *
from synbioweaver.aspects.lagLogisticGrowthAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *
from synbioweaver.aspects.writeSBMLModelAspect import *
from synbioweaver.aspects.runCudaSim import *
from synbioweaver.aspects.writeCudaFileODE import *
from synbioweaver.aspects.writeABCInputFileODE import *
from synbioweaver.aspects.writeCudaFileGillespie import *
from synbioweaver.aspects.writeABCInputFileGillespie import *

declareNewMolecule('GFP')

class constGFP(Circuit):
    def mainCircuit(self):
        # constitutive promoter
        self.addPart(ConstitutivePromoter)
        self.addPart(CodingRegion(GFP))

# Take the compiled design and add Type Advice that generates a mass action model
#compiledDesign1 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()
#compiledDesign1.printReactionNetwork()

#compiledDesign1 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, ExpGrowthAspect, PrintReactionNetwork).output()
#compiledDesign1.printReactionNetwork()

#compiledDesign1 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork, WriteSBMLModel).output()
#compiledDesign1.printReactionNetwork()
#modelStr = compiledDesign1.writeSBMLModel()

#print "####### circuit + exp growth"
#compiledDesign1 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork, ExpGrowthAspect, WriteSBMLModel).output()
#compiledDesign1.printReactionNetwork()
#modelStr = compiledDesign1.writeSBMLModel()

# write an SBML model out
#sbmlFile = open("constGFP_exp.sbml","w")
#print >>sbmlFile, modelStr
#sbmlFile.close()

#print "####### circuit + logistic growth"
#compiledDesign1 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork, LogisticGrowthAspect, WriteSBMLModel).output()
#compiledDesign1.printReactionNetwork()
#modelStr = compiledDesign1.writeSBMLModel()

# write an SBML model out
#sbmlFile = open("constGFP_log.sbml","w")
#print >>sbmlFile, modelStr
#sbmlFile.close()

#print "####### circuit + lag logistic growth"
#compiledDesign1 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork, LagLogisticGrowthAspect, WriteSBMLModel).output()
#compiledDesign1.printReactionNetwork()
#modelStr = compiledDesign1.writeSBMLModel()

# write an SBML model out
#sbmlFile = open("constGFP_lag.sbml","w")
#print >>sbmlFile, modelStr
#sbmlFile.close()


print "####### circuit + logistic growth"
compiledDesign5 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork, LogisticGrowthAspect, WriteABCInputFileODE, WriteCudaFileODE, RunCudaSim).output()
compiledDesign5.printReactionNetwork()
compiledDesign5.writeABCInputFileODE()
compiledDesign5.writeCudaFileODE()
compiledDesign5.runCudaSim()
