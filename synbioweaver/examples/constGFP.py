from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.MassActionKineticsRNAAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

class constGFP(Circuit):
    def mainCircuit(self):

        # constitutive promoter
        declareNewMolecule('GFP')
        declareNewMolecule('RFP')
        declareNewPart('Pc', ConstitutivePromoter)
        self.addPart(Pc)
        self.addPart(CodingRegion(GFP))
        self.addPart(CodingRegion(RFP))

        

# Take the compiled design and add Type Advice that generates a mass action model
compiledDesign1 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()
print "######################## Design1:"
compiledDesign1.printReactionNetwork()

compiledDesign2 = Weaver(constGFP, DesignRules, PromoterMapping, MassActionKineticsRNA, PrintReactionNetwork).output()
print "######################## Design2:"
compiledDesign2.printReactionNetwork()

compiledDesign3 = Weaver(constGFP, DesignRules, PromoterMapping, SheaAckersKineticsRNA, PrintReactionNetwork).output()
print "######################## Design3:"
compiledDesign3.printReactionNetwork()
