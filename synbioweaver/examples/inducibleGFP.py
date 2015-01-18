from synbioweaver.core import *
from reactionNetworksAspect import *
from databaseSearchAspect import *

declareNewMolecule('GFP')
declareNewMolecule('aTc')

# constitutive promoter
declareNewPart('PTet', PositivePromoter, [aTc])

class inducibleGFP(Circuit):
    def mainCircuit(self):

        self.createMolecule(aTc)
        
        self.addPart(PTet)
        self.addPart(CodingRegion(GFP))
        
compiledDesign1 = Weaver(inducibleGFP, ReactionNetworks, SheaAckersKineticsRNA, PrintReactionNetwork).output()
print "######################## Design2:"
compiledDesign1.printReactionNetwork()

compiledDesign2 = Weaver(inducibleGFP, ReactionNetworks, PromoterSearch, SheaAckersKineticsRNA, PrintReactionNetwork).output()
print "######################## Design3:"
compiledDesign2.printReactionNetwork()

