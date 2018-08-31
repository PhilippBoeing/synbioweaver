from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import DesignRules
from transferLibraryTamsir2010 import MoleculeExpressionTrace

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

declareNewMolecule('aTc')
declareNewMolecule('CFP')
declareNewMolecule('YFP')
declareNewMolecule('cI')

declareNewPart('PR', NegativePromoter, [cI])

class WildType(Circuit):
    def mainCircuit(self):
        self.createMolecule(aTc)
        self.addPart(PositivePromoter(aTc))
        self.addPart(CodingRegion(cI))
        self.addPart(CodingRegion(YFP))
        self.addPart(PR)
        self.addPart(CodingRegion(CFP))

##
print "Running SystemCircuit1"
compiledSystem = Weaver(WildType,DesignRules).output()
print compiledSystem
