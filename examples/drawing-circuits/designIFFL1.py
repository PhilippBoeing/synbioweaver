from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.printStackAspect import *
from synbioweaver.aspects.pigeonOutputAspect import *

declareNewMolecule('LuxR')
declareNewMolecule('TetR')
declareNewMolecule('GFP')

declareNewPart('t1',Terminator)
declareNewPart('t2',Terminator)
declareNewPart('t3',Terminator)
declareNewPart('r1',RBS )
declareNewPart('r2',RBS )
declareNewPart('r3',RBS )
declareNewPart('cLuxR',CodingRegion,moleculesAfter=[LuxR])
declareNewPart('cTetR',CodingRegion,moleculesAfter=[TetR])
declareNewPart('cGFP',CodingRegion,moleculesAfter=[GFP])
declareNewPart('Pc', ConstitutivePromoter)
declareNewPart('P_lux', PositivePromoter, [LuxR] )
declareNewPart('P_luxtet', HybridPromoter, [LuxR,TetR], regulatorInfoMap={LuxR:True,TetR:False} )

declareNewMolecule('aTc')
declareNewMolecule('TetR_aTc')
declareNewMolecule('AHL')
declareNewMolecule('LuxR_AHL')
declareNewMolecule('zero')

class simpleCircuit(Circuit):
    def mainCircuit(self):
        self.createMolecule(aTc)
        self.createMolecule(TetR_aTc)
        self.createMolecule(AHL)
        self.createMolecule(LuxR_AHL)
        self.createMolecule(zero)

        self.addPart(Pc)
        self.addPart(r1)
        self.addPart(cLuxR)
        self.addPart(t1)
        self.addPart(P_lux)
        self.addPart(r2)
        self.addPart(cTetR)
        self.addPart(t2)
        self.addPart(P_luxtet)
        self.addPart(r3)
        self.addPart(cGFP)
        self.addPart(t3)
        
        self.reactionFrom(aTc, TetR) >> self.reactionTo( TetR_aTc )
        self.reactionFrom(TetR_aTc) >> self.reactionTo( zero )
        self.reactionFrom(AHL, LuxR) >> self.reactionTo( LuxR_AHL )

compiledDesign = Weaver(simpleCircuit, PigeonOutput).output()

compiledDesign.printPigeonOutput()
