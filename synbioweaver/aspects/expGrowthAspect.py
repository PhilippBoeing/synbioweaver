from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os, copy

class ExpGrowthAspect(Aspect):
    
    def mainAspect(self):
        self.addWeaverOutput(self.getReactions)

    def getReactions(self,weaverOutput):
        self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix = weaverOutput.getReactions()

        return [self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix]

    
