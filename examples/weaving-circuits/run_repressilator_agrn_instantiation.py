from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from instantiateAbstractMoleculesAspect import *

import itertools

# Abstract repressilator definitions
declareNewMolecule('AbstractMolecule')

declareNewMolecule('A',AbstractMolecule)
declareNewMolecule('B',AbstractMolecule)
declareNewMolecule('C',AbstractMolecule)

class RepAGRN(Circuit):
    # an abstract repressilator circuit
    #       A(C), B(A), C(B)
           
    def declareRepressilatorProteins(self):
        self.repressilationProteins = []
        self.repressilationProteins.append(A)
        self.repressilationProteins.append(B)
        self.repressilationProteins.append(C)
        
    def getSignal(self):
        return A
        
    def mainCircuit(self):
        self.createMolecule(A)
        self.createMolecule(B)
        self.createMolecule(C)
        self.declareRepressilatorProteins()
        lastprotein = self.repressilationProteins[-1]
        
        for protein in self.repressilationProteins:
            self.addPart(NegativePromoter(lastprotein))
            self.addPart(CodingRegion(protein))
            lastprotein = protein
            
repressors = [declareNewMolecule('TetR'),declareNewMolecule('LacL'),declareNewMolecule('cl')]

print Weaver(RepAGRN).output()

print "Instantiate by adding in repressors in coding sequence order"
print Weaver(RepAGRN,InstantiateAbstractMolecules(repressors)).output()

print "Generate all combinations"
for p in list(itertools.permutations(repressors,3)):
    print list(p)
    print Weaver(RepAGRN,InstantiateAbstractMolecules(list(p))).output()
