from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import DesignRules
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

# This example is identical to the tutorial_molecules.py example only a mass action model is generated
# Note the slightly different syntax: promoters and their type and regulators must be defined using declareNewPart


declareNewMolecule('MoleculeA')
declareNewMolecule('MoleculeB')
declareNewMolecule('MoleculeAB')
declareNewMolecule('GFP')
declareNewPart('P', PositivePromoter, [MoleculeAB])
declareNewPart('Pc', ConstitutivePromoter)

class MoleculeCircuit(Circuit):
    def mainCircuit(self):

        # molecule instances in this circuit must be created
        self.createMolecule(MoleculeA)
        self.createMolecule(MoleculeB)
        self.createMolecule(MoleculeAB)

        self.addPart(Pc)
        self.addPart(CodingRegion(MoleculeA))
        self.addPart(Pc)
        self.addPart(CodingRegion(MoleculeB))

        # add a reaction A + B -> AB
        self.reactionFrom(MoleculeA, MoleculeB) >> self.reactionTo(MoleculeAB)

        self.addPart(P)
        self.addPart(CodingRegion(GFP))


compiledDesign = Weaver(MoleculeCircuit,DesignRules, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()
print compiledDesign.printReactionNetwork()
