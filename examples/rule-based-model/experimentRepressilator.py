from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import DesignRules
from synbioweaver.aspects.pysbRulesAspect import PYSBmodel

# transcription factor used in the repressilator system
declareNewMolecule('LacI')
declareNewMolecule('TetR')
declareNewMolecule('lCl')

declareNewPart('BBaB0011',Terminator)
declareNewPart('BBaB0034',RBS)
declareNewPart('BBaC0012',CodingRegion,moleculesAfter=[LacI])
declareNewPart('BBaC0040',CodingRegion,moleculesAfter=[TetR])
declareNewPart('BBaC0051',CodingRegion,moleculesAfter=[lCl])
declareNewPart('BBaR0010',NegativePromoter,moleculesBefore=[LacI])
declareNewPart('BBaR0040',NegativePromoter,moleculesBefore=[TetR])
declareNewPart('BBaR0051',NegativePromoter,moleculesBefore=[lCl])




class Repressilator(Circuit):
    def mainCircuit(self):
        self.createMolecule(LacI)
        self.addPart(BBaR0010)
        self.addPart(BBaC0040)
        self.addPart(BBaR0040)
        self.addPart(BBaC0051)
        self.addPart(BBaR0051)
        self.addPart(BBaC0012)


compiledDesign = Weaver(Repressilator, DesignRules, PYSBmodel).output()
print compiledDesign

compiledDesign.runKappaSimulation()

