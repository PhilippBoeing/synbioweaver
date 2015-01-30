from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

declareNewMolecule('TetR')
declareNewMolecule('TetRn2')
declareNewMolecule('AHL')
declareNewMolecule('AHLn4')
declareNewMolecule('GFP')
declareNewMolecule('zero')
declareNewPart('PTet', NegativePromoter, [TetRn2])
declareNewPart('PLux', PositivePromoter, [AHLn4] )

# the sender circuit: luxI is under the control of PTet
class sn(Circuit):
    def mainCircuit(self):
        self.createMolecule(TetR)
        self.createMolecule(TetRn2)
        self.createMolecule(AHL) 
        self.createMolecule(AHLn4) 

        self.addPart(PTet)
        self.addPart(CodingRegion(AHL))
       
        self.reactionFrom(AHL, AHL, AHL, AHL) >> self.reactionTo( AHLn4 )
        self.reactionFrom( AHLn4 ) >> self.reactionTo( zero )
        
        self.reactionFrom(TetR, TetR) >> self.reactionTo( TetRn2 )
        self.reactionFrom( TetRn2 ) >> self.reactionTo( zero )
        
        self.exportMolecule(AHLn4)

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
        self.addCircuit(sn)
        self.addCircuit(rc)

#compiledDesign = Weaver(ExpressCircuit, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()
#compiledDesign = Weaver(GFPCircuit, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()

compiledDesign = Weaver(System, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()
#compiledDesign = Weaver(System, PromoterMapping, SheaAckersKineticsRNA, PrintReactionNetwork).output()

#print compiledDesign
compiledDesign.printReactionNetwork()
