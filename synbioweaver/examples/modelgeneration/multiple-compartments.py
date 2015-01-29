from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

declareNewMolecule('TetR')
declareNewMolecule('TetRn2')
declareNewMolecule('AHL')
declareNewMolecule('AHLn2')
declareNewMolecule('AHLn4')
declareNewMolecule('GFP')
declareNewPart('PTet', NegativePromoter, [TetRn2])
declareNewPart('PLux', PositivePromoter, [AHLn4] )

# the sender circuit: luxI is under the control of PTet
class sn(Circuit):
    def mainCircuit(self):
        self.createMolecule(TetR)
        self.createMolecule(TetRn2)
        self.createMolecule(AHL) 
        self.addPart(PTet)
        self.addPart(CodingRegion(AHL))
        self.exportMolecule(AHL)

        self.reactionFrom(TetR, TetR) >> self.reactionTo( TetRn2 )

# the GFP Circuit compartment expresses GFP induced by Molecule A.
# it can import Molecule A from its environment
class rc(Circuit):
    def mainCircuit(self):
        self.createMolecule(GFP)
        self.importMolecule(AHLn4)
        self.addPart(PLux)
        self.addPart(CodingRegion(GFP))

# the System Compartment encapsules the two previous compartments.
# It contains Molecule A, secreted by Express Circuit
class System(Circuit):
    def mainCircuit(self):
        self.createMolecule(AHL) 
        self.createMolecule(AHLn2) 
        self.createMolecule(AHLn4) 
        self.reactionFrom(AHL,AHL) >> self.reactionTo( AHLn2 )
        self.reactionFrom(AHLn2,AHLn2) >> self.reactionTo( AHLn4 )

        self.addCircuit(sn)
        self.addCircuit(rc)

#compiledDesign = Weaver(ExpressCircuit, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()
#compiledDesign = Weaver(GFPCircuit, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()

compiledDesign = Weaver(System, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()
#compiledDesign = Weaver(System, PromoterMapping, SheaAckersKineticsRNA, PrintReactionNetwork).output()

#print compiledDesign
compiledDesign.printReactionNetwork()
