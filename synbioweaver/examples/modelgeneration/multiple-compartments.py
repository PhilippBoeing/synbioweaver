from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

declareNewMolecule('A')
declareNewMolecule('B')
declareNewMolecule('GFP')
declareNewPart('Pn', NegativePromoter, [A])
declareNewPart('Pc', ConstitutivePromoter )

# the express circuit expresses molecule A constitutively and secretes it
class dev1(Circuit):
    def mainCircuit(self):
        self.createMolecule(A)
        self.addPart(Pc)
        self.addPart(CodingRegion(A))
        self.exportMolecule(A)

# the GFP Circuit compartment expresses GFP induced by Molecule A.
# it can import Molecule A from its environment
class dev2(Circuit):
    def mainCircuit(self):
        self.createMolecule(GFP)
        #self.createMolecule(A)
        self.importMolecule(A)
        self.addPart(Pn)
        self.addPart(CodingRegion(GFP))

# the System Compartment encapsules the two previous compartments.
# It contains Molecule A, secreted by Express Circuit
class System(Circuit):
    def mainCircuit(self):
        self.addCircuit(dev1)
        self.addCircuit(dev2)

#compiledDesign = Weaver(ExpressCircuit, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()
#compiledDesign = Weaver(GFPCircuit, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()
compiledDesign = Weaver(System, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()

#print compiledDesign
compiledDesign.printReactionNetwork()
