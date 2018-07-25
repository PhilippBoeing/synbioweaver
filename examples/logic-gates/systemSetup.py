from synbioweaver import *

declareNewMolecule('Cl')
declareNewMolecule('Ara')
declareNewMolecule('aTc')
declareNewMolecule('LasI')
declareNewMolecule('RhlI')
declareNewMolecule('YFP')

def buildAbstractNorGate(circuitSelf,input1,input2,output):
    circuitSelf.addPart(PositivePromoter(input1))
    circuitSelf.addPart(PositivePromoter(input2))
    circuitSelf.addPart(CodingRegion(Cl))
    circuitSelf.addPart(NegativePromoter(Cl))
    circuitSelf.addPart(CodingRegion(output))

class Cell1(Circuit):
    def mainCircuit(self):
        self.importMolecule(Ara)
        self.importMolecule(aTc)
        buildAbstractNorGate(self,Ara,aTc,LasI)
        self.exportMolecule(LasI)

class Cell2(Circuit):
    def mainCircuit(self):
        self.importMolecule(Ara)
        self.importMolecule(LasI)
        buildAbstractNorGate(self,Ara,LasI,RhlI)
        self.exportMolecule(RhlI)

class Cell3(Circuit):
    def mainCircuit(self):
        self.importMolecule(LasI)
        self.importMolecule(aTc)
        buildAbstractNorGate(self,LasI,aTc,RhlI)
        self.exportMolecule(RhlI)

class Cell4(Circuit):
    def mainCircuit(self):
        self.importMolecule(RhlI)
        self.addPart(PositivePromoter(RhlI))
        self.addPart(CodingRegion(YFP))

class SystemCircuit(Circuit):
    def mainCircuit(self):
        self.createMolecule(Ara)
        self.createMolecule(aTc)
        self.addCircuit(Cell1)
        self.addCircuit(Cell2)
        self.addCircuit(Cell3)
        self.addCircuit(Cell4)