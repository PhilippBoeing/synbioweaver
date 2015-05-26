from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os
from copy import *

class ExpGrowth(Aspect):
    
    def mainAspect(self):
        ExpGrowth.builtReactions = False
        self.addWeaverOutput(self.getContext)

    def getContext(self,weaverOutput):

        if ExpGrowth.builtReactions == False:
            # first access the existing reactions
            if getattr(weaverOutput, "getReactions", None) != None:
                self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getReactions()
            else:
                print "ExpGrowth : getReactions() is available. Quitting"
                exit()

            self.addExpGrowth()

            ExpGrowth.builtReactions = True

        return [deepcopy(self.nspecies), deepcopy(self.nreactions), deepcopy(self.species), deepcopy(self.reactions), deepcopy(self.stoichiometry_matrix), deepcopy(self.parameters) ]

    
    def addExpGrowth(self):
        # modify the reactions to incorporate exponential growth
        # this essentially means converting any reactions of the form 
        #     0 -> P, rate k
        # to
        #     N -> X + P, rate N*k
        #     N -> 2N
        #     N -> 0

        # in this particular example we should be able to modify the reactions in place
        # this probably isn't true in general

        print self.parameters
        for r in self.reactions:
            if r.process == "proteinExp" or r.process == "rnaTransc":
                r.products.append("N")
                r.reactants.append("N")
                r.rate = r.rate+"*N"

        # add new reactions
        newreac =  Reaction(["N"], ["N","N"], "context")
        newreac.assignMassAction()
        self.reactions.append( newreac )
        self.nreactions += 1
        for k in newreac.param:
            self.parameters.append( k )

        # add new species
        self.species.append("N")
        self.nspecies += 1

        # recalculate stoichiometry matrix
        self.stoichiometry_matrix = stoichiometry(self.nspecies, self.nreactions, self.species, self.reactions)
