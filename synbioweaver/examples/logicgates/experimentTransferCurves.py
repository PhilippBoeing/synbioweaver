from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import DesignRules
from synbioweaver.aspects.transferFunctionAspect import MoleculeExpressionTrace
from systemSetup import SystemCircuit

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

compiledSystem = Weaver(SystemCircuit,DesignRules,MoleculeExpressionTrace).output()

# create "inspector" variables for relevant molecules
inputMol1 = compiledSystem.moleculeList[0]
inputMol2 = compiledSystem.moleculeList[1]
outputMol = compiledSystem.wovenCircuitList[3].moleculeList[1]

araInputRange = np.logspace(-4, 2,100, endpoint=True)
atcInputRange = np.logspace(4, -2,100, endpoint=True)

yfpOutput= []

for i in range(0,len(atcInputRange)):
    yfpOutput.append([])
    for j in range(0,len(araInputRange)):
        inputMol1.signal = araInputRange[j]
        inputMol2.signal = atcInputRange[i]
        yfpOutput[i].append(outputMol.calculateTransfer())

fig1 = plt.figure()
p = plt.imshow(yfpOutput,interpolation='none',norm=matplotlib.colors.LogNorm())

# turn off ticks
plt.xticks(np.arange(0,0,10),np.arange(0,0,1))
plt.yticks(np.arange(0,0,10),np.arange(0,0,1))

plt.xlabel('Input 1 (Ara)')
plt.ylabel('Input 2 (aTc)')

plt.jet()
cb = plt.colorbar()
cb.set_label('Output (YFP)')


plt.show()

