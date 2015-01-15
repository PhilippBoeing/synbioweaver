from synbioweaver import *
from synbioweaver.designaspects import DesignRules

declareNewMolecule('MoleculeA',Protein)
declareNewMolecule('MoleculeB')
declareNewMolecule('MoleculeAB')
declareNewMolecule('GFP')

class MoleculeCircuit(Circuit):
    def mainCircuit(self):
        self.addPart(Promoter)
        self.addPart(CodingRegion(MoleculeA() + MoleculeB()))
        self.addPart(Promoter)
        self.addPart(CodingRegion(MoleculeB))

        # of course, this reaction could also be added by an aspect
        # for example, an aspect could check if MoleculeAB is needed
        # and insert the reaction
        self.reactionFrom(MoleculeA, MoleculeB) >> self.reactionTo(MoleculeAB)

        # should declare molecule?
        self.addPart(PositivePromoter(MoleculeAB))
        self.addPart(CodingRegion(GFP))




compiledDesign = Weaver(MoleculeCircuit,DesignRules).output()
print compiledDesign
