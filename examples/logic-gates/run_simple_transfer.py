from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import DesignRules
from transferFunctionAspect import MoleculeExpressionTrace

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

declareNewMolecule('Cl')
declareNewMolecule('Ara')
declareNewMolecule('aTc')
declareNewMolecule('LasI')
declareNewMolecule('RhlI')
declareNewMolecule('YFP')

class SystemCircuit1(Circuit):
    def mainCircuit(self):
        self.createMolecule(aTc)
        self.addPart(PositivePromoter(aTc))
        self.addPart(CodingRegion(YFP))

class SystemCircuit2(Circuit):
    def mainCircuit(self):
        self.createMolecule(aTc)
        self.addPart(PositivePromoter(aTc))
        self.addPart(CodingRegion(Cl))
        self.addPart(NegativePromoter(Cl))
        self.addPart(CodingRegion(YFP))
        
##
print "Running SystemCircuit1"
compiledSystem1 = Weaver(SystemCircuit1,DesignRules,MoleculeExpressionTrace).output()
inputMol1 = compiledSystem1.moleculeList[0]
outputMol = compiledSystem1.moleculeList[1]
print "input/outputs:", inputMol1, outputMol

molInputRange = np.logspace(-4, 3,100, endpoint=True)
yfpOutput= []

for i in range(0,len(molInputRange)):
    inputMol1.signal = molInputRange[i]
    yfpOutput.append( outputMol.calculateTransfer() )

fig = plt.figure()
plt.plot( np.log10(molInputRange), np.log10(yfpOutput) )
plt.xlabel('Input (aTc)')
plt.ylabel('Ouput (YFP)')
plt.savefig("plot-transfer-system1.pdf", bbox_inches='tight')
plt.close()

##
print "Running SystemCircuit2"
compiledSystem2 = Weaver(SystemCircuit2,DesignRules,MoleculeExpressionTrace).output()
inputMol1 = compiledSystem2.moleculeList[0]
outputMol = compiledSystem2.moleculeList[2]
print "input/outputs:", inputMol1, outputMol

molInputRange = np.logspace(-4, 3,100, endpoint=True)
yfpOutput= []

for i in range(0,len(molInputRange)):
    inputMol1.signal = molInputRange[i]
    yfpOutput.append( outputMol.calculateTransfer() )

fig = plt.figure()
plt.plot( np.log10(molInputRange), np.log10(yfpOutput) )
plt.xlabel('Input (aTc)')
plt.ylabel('Ouput (YFP)')
plt.savefig("plot-transfer-system2.pdf", bbox_inches='tight')
plt.close()
