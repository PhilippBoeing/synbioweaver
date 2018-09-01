from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.MassActionKineticsRNAAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

declareNewMolecule('GFP')
declareNewMolecule('RFP')
declareNewPart('Pc', ConstitutivePromoter)

class circuit1(Circuit):
    def mainCircuit(self):
        self.addPart(Pc)
        self.addPart(RBS)
        self.addPart(CodingRegion(GFP))
        self.addPart(Terminator)

        self.addPart(Pc)
        self.addPart(RBS)
        self.addPart(CodingRegion(RFP))
        self.addPart(Terminator)
        

class circuitGFP(Circuit):
    def mainCircuit(self):
        self.addPart(Pc)
        self.addPart(RBS)
        self.addPart(CodingRegion(GFP))
        self.addPart(Terminator)
        
class circuitRFP(Circuit):
    def mainCircuit(self):
        self.addPart(Pc)
        self.addPart(RBS)
        self.addPart(CodingRegion(RFP))
        self.addPart(Terminator)

class circuit2(Circuit):
    def mainCircuit(self):
        self.addCircuit(circuitRFP)
        self.addCircuit(circuitGFP)

# Take the compiled design and add Type Advice that generates a mass action model
compiledDesign1 = Weaver(circuit1, PromoterMapping, MassActionKineticsRNA, PrintReactionNetwork).output()
print "######################## Design1:"
compiledDesign1.printReactionNetwork()

compiledDesign2 = Weaver(circuit2, PromoterMapping, MassActionKineticsRNA, PrintReactionNetwork).output()
print "######################## Design2:"
compiledDesign2.printReactionNetwork()

