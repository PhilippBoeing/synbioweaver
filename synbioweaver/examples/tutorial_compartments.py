from synbioweaver import *
from synbioweaver.designaspects import DesignRules

declareNewMolecule('MoleculeA')
declareNewMolecule('MoleculeB')
declareNewMolecule('MoleculeAB')
declareNewMolecule('GFP')

# the express circuit expresses molecule A constitutively and secretes it
class ExpressCircuit(Circuit):
    def mainCircuit(self):
        self.addPart(ConstitutivePromoter)
        self.addPart(CodingRegion(MoleculeA))
        self.exportMolecule(MoleculeA)

# the GFP Circuit compartment expresses GFP induced by Molecule A.
# it can import Molecule A from its environment
class GFPCircuit(Circuit):
    def mainCircuit(self):
        self.importMolecule(MoleculeA)
        self.addPart(NegativePromoter(MoleculeA))
        self.addPart(CodingRegion(GFP))

# the System Compartment encapsules the two previous compartments.
# It contains Molecule A, secreted by Express Circuit
class System(Circuit):
    def mainCircuit(self):
        self.addCircuit(ExpressCircuit)
        self.addCircuit(GFPCircuit)

compiledDesign = Weaver(System,DesignRules).output()
print compiledDesign



