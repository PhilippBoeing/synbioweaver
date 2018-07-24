from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsRNAAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

# this example is based on the inducible lux quorum sensing system
declareNewMolecule('AHL')
declareNewMolecule('AHLn4')
declareNewMolecule('R')
declareNewMolecule('R_AHLn4')
declareNewMolecule('GFP')
declareNewMolecule('zero')
declareNewPart('Plux', PositivePromoter, [R_AHLn4])

class LuxOperon(Circuit):
    def mainCircuit(self):
        self.createMolecule(AHL) 
        self.createMolecule(AHLn4) 
        self.createMolecule(R)
        self.createMolecule(R_AHLn4)
        
        # positive promoters
        self.addPart(Plux)
        self.addPart(RBS)
        self.addPart(CodingRegion(R))
        self.addPart(CodingRegion(AHL))
        self.addPart(Terminator)
        
        declareNewPart('Pc', NegativePromoter, [R])
        self.addPart(Pc)
        self.addPart(CodingRegion(GFP))
        self.addPart(Terminator)

        self.reactionFrom(AHL,AHL,AHL,AHL) >> self.reactionTo( AHLn4 )
        self.reactionFrom(R,AHLn4) >> self.reactionTo( R_AHLn4 )
        self.reactionFrom(AHLn4) >> self.reactionTo( zero )


# Take the compiled design and add Type Advice that generates a mass action model
compiledDesign1 = Weaver(LuxOperon, PromoterMapping, MassActionKineticsRNA, PrintReactionNetwork).output()
compiledDesign1.printReactionNetwork()

compiledDesign2 = Weaver(LuxOperon, PromoterMapping, SheaAckersKineticsRNA, PrintReactionNetwork).output()
compiledDesign2.printReactionNetwork()
