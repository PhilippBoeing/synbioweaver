from synbioweaver.core import *
from synbioweaver.aspects.modelDefinitions import *
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
                self.model = weaverOutput.getReactions()
            else:
                print "ExpGrowth : getReactions() is not available. Quitting"
                exit()

            self.addExpGrowth()

            ExpGrowth.builtReactions = True

        return deepcopy(self.model)

    
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

        print self.model.parameters
        for r in self.model.reactions:
            if r.process == "proteinExp" or r.process == "rnaTransc":
                r.products.append( Species("con","N") )
                r.reactants.append( Species("con","N") )
                r.rate = r.rate+"*N"

        # add new reactions
        newreac =  Reaction([Species("con","N")], [Species("con","N"),Species("con","N")], "context")
        newreac.assignMassAction()
        self.model.reactions.append( newreac )
        self.model.nreactions += 1
        for k in newreac.param:
            self.model.parameters.append( k )

        # add new species
        self.model.species.append(Species("con","N"))
        self.model.nspecies += 1

        # recalculate stoichiometry matrix
        self.model.stoichiometry_matrix = self.model.stoichiometry(self.model.nspecies, self.model.nreactions, self.model.species, self.model.reactions)
