from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.MassActionKineticsRNAAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

declareNewMolecule('GFP')
declareNewMolecule('LacI')
declareNewMolecule('AraC')
declareNewPart("P1", HybridPromoter, [LacI,AraC], regulatorInfoMap={LacI:False, AraC:True})
declareNewPart("P2", HybridPromoter, [LacI,AraC], regulatorInfoMap={LacI:False, AraC:True})
declareNewPart("P3", HybridPromoter, [LacI,AraC], regulatorInfoMap={LacI:False, AraC:True})

class strickerOsc(Circuit):
    def mainCircuit(self):

        self.createMolecule(GFP)
        self.createMolecule(LacI)
        self.createMolecule(AraC)

        self.addPart(P1)
        self.addPart(CodingRegion(AraC))

        self.addPart(P2)
        self.addPart(CodingRegion(LacI))

        self.addPart(P3)
        self.addPart(CodingRegion(GFP))
        

        
compiledDesign = Weaver(strickerOsc, DesignRules, PromoterMapping, SheaAckersKineticsRNA, PrintReactionNetwork).output()
print "######################## Stricker Design:"
compiledDesign.printReactionNetwork()
