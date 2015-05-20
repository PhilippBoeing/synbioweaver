from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import DesignRules

declareNewMolecule('MoleculeA')
declareNewMolecule('MoleculeB')
declareNewMolecule('MoleculeAB')
declareNewMolecule('GFP')

class MoleculeCircuit(Circuit):
    def mainCircuit(self):

        # molecule instances in this circuit must be created
        self.createMolecule(MoleculeA)
        self.createMolecule(MoleculeB)
        self.createMolecule(MoleculeAB)

        self.addPart(Promoter)
        self.addPart(CodingRegion(MoleculeA))
        self.addPart(Promoter)
        self.addPart(CodingRegion(MoleculeB))

        # add a reaction A + B -> AB
        self.reactionFrom(MoleculeA, MoleculeB) >> self.reactionTo(MoleculeAB)

        self.addPart(PositivePromoter(MoleculeAB))
        self.addPart(CodingRegion(GFP))




compiledDesign = Weaver(MoleculeCircuit,DesignRules).output()
print compiledDesign
