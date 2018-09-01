from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import DesignRules
from transferLibraryRosenfeld2005Aspect import LibraryRosenfeld2005

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def plot_transfer_curve(inputMol,outputMol,filename):
    print "input/output:", inputMol, outputMol

    molInputRange = np.logspace(1, 3,100, endpoint=True)
    Output= []

    for i in range(0,len(molInputRange)):
        inputMol.signal = molInputRange[i]
        Output.append( outputMol.calculateTransfer() )

    fig = plt.figure()
    plt.plot( np.log10(molInputRange), np.log10(Output) )
    plt.xlabel('Input')
    plt.ylabel('Ouput')
    plt.savefig(filename, bbox_inches='tight')
    plt.close()

declareNewMolecule('CFP')
declareNewMolecule('CI')

declareNewPart('PR', NegativePromoter, [CI])
declareNewPart('PR_OR2', NegativePromoter, [CI])

class circuitPR(Circuit):
    def mainCircuit(self):
        self.createMolecule(CI)
        self.addPart(PR)
        self.addPart(CodingRegion(CFP))

class circuitPR_OR2(Circuit):
    def mainCircuit(self):
        self.createMolecule(CI)
        self.addPart(PR_OR2)
        self.addPart(CodingRegion(CFP))

##
print "Running wildtype"
compiledSystem = Weaver(circuitPR,LibraryRosenfeld2005).output()
print compiledSystem
plot_transfer_curve(compiledSystem.moleculeList[0],compiledSystem.moleculeList[1],"plot-rosenfeld-transfer-PR.pdf")

##
print "Running mutant"
compiledSystem = Weaver(circuitPR_OR2,LibraryRosenfeld2005).output()
print compiledSystem
plot_transfer_curve(compiledSystem.moleculeList[0],compiledSystem.moleculeList[1],"plot-rosenfeld-transfer-PR-OR2.pdf")

## simple example of how parts can be mutated
class MutatePR(Aspect):
    def mainAspect(self):
        self.addAdvice(PointCut(PartSignature('*.PR+(CI)'),PointCut.REPLACE),self.insertPR_OR2)

    def insertPR_OR2(self,context):
        self.addPart(PR_OR2)
        return True

print "Running mutation"
compiledSystem = Weaver(circuitPR,MutatePR,LibraryRosenfeld2005).output()
print compiledSystem
plot_transfer_curve(compiledSystem.moleculeList[0],compiledSystem.moleculeList[1],"plot-rosenfeld-transfer-PR-mutated.pdf")
