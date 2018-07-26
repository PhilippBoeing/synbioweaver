from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.printStackAspect import *
from pigeonOutputAspect import *

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
declareNewPart('P_LtetO', NegativePromoter, [TetR] )
declareNewPart('P_luxtet', HybridPromoter, [LuxR,TetR], regulatorInfoMap={LuxR:True,TetR:False} )

class simpleCircuit(Circuit):
    def mainCircuit(self):

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
        
#compiledDesign = Weaver(constGFP, DesignRules, PrintStack, PigeonOutput).output()
compiledDesign = Weaver(simpleCircuit, PigeonOutput).output()

compiledDesign.printPigeonOutput()
