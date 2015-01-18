from synbioweaver.core import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.databaseSearchAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

declareNewMolecule('GFP')
declareNewMolecule('aTc')

# constitutive promoter
declareNewPart('PTet', PositivePromoter, [aTc])

class inducibleGFP(Circuit):
    def mainCircuit(self):

        self.createMolecule(aTc)
        
        self.addPart(PTet)
        self.addPart(CodingRegion(GFP))
        
compiledDesign1 = Weaver(inducibleGFP, PromoterMapping, SheaAckersKineticsRNA, PrintReactionNetwork).output()
print "######################## Design1:"
compiledDesign1.printReactionNetwork()

compiledDesign2 = Weaver(inducibleGFP, PromoterMapping, PromoterSearch, SheaAckersKineticsRNA, PrintReactionNetwork).output()
print "######################## Design2:"
compiledDesign2.printReactionNetwork()

