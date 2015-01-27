from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.MassActionKineticsRNAAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.expGrowthAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

declareNewMolecule('GFP')

class constGFP(Circuit):
    def mainCircuit(self):
        # constitutive promoter
        self.addPart(ConstitutivePromoter)
        self.addPart(CodingRegion(GFP))

# Take the compiled design and add Type Advice that generates a mass action model
compiledDesign1 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, ExpGrowthAspect, PrintReactionNetwork).output()
print "######################## Design1:"
compiledDesign1.printReactionNetwork()

compiledDesign2 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsRNA,  ExpGrowthAspect, PrintReactionNetwork).output()
print "######################## Design2:"
compiledDesign2.printReactionNetwork()
