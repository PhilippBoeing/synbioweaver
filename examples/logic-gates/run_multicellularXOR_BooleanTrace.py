from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import DesignRules
from booleanTraceAspect import BooleanLogicTrace
from systemSetup import SystemCircuit

compiledSystem = Weaver(SystemCircuit,DesignRules,BooleanLogicTrace).output()

# create "inspector" variables for relevant molecules
inputMol1 = compiledSystem.moleculeList[0]
inputMol2 = compiledSystem.moleculeList[1]
outputMol = compiledSystem.wovenCircuitList[3].moleculeList[1]

# output truth table:
print "Input: False / False"
print "Output: " + str(outputMol.calculateExpressed())

print "Input: True / False"
inputMol1.moleculeExpressed = True
print "Output: " + str(outputMol.calculateExpressed())

print "Input: False / True"
inputMol1.moleculeExpressed = False
inputMol2.moleculeExpressed = True
print "Output: " + str(outputMol.calculateExpressed())


print "Input: True / True"
inputMol1.moleculeExpressed = True
print "Output: " + str(outputMol.calculateExpressed())

