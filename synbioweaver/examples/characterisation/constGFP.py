from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.MassActionKineticsRNAAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.expGrowthAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *
from synbioweaver.aspects.writeSBMLModelAspect import *

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

compiledDesign1 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork, ExpGrowthAspect, WriteSBMLModel).output()
compiledDesign1.printReactionNetwork()
modelStr = compiledDesign1.writeSBMLModel()

# write an SBML model out
sbmlFile = open("constGFP.sbml","w")
print >>sbmlFile, modelStr
sbmlFile.close()
