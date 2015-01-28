from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os, copy

class LogisticGrowthAspect(Aspect):
    
    def mainAspect(self):
        LogisticGrowthAspect.builtReactions = False
        self.addWeaverOutput(self.getContext)

    def getContext(self,weaverOutput):
        
        # first access the existing reactions
        if getattr(weaverOutput, "getReactions", None) != None:
            self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getReactions()
        else:
            print "LogisticGrowthAspect : getReactions() is available. Quitting"
            exit()

        if LogisticGrowthAspect.builtReactions == False:
            self.addLogisticGrowth()

            LogisticGrowthAspect.builtReactions = True

        return [self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters]

    
    def addLogisticGrowth(self):
        # modify the reactions to incorporate logistic growth
        # this essentially means converting any reactions of the form 
        #     0 -> P, rate k
        # to
        #     N -> N + P, rate kN
        #     N -> 2N,    rate p1*N(1-N/p2)

        # in this particular example we should be able to modify the reactions in place
        # this probably isn't true in general

        # modify exixting reactions : change products, reactants and rate - no new parameters
        for r in self.reactions:
            if r.process == "proteinExp" or r.process == "rnaTransc":
                r.products.append("N")
                r.reactants.append("N")
                r.rate = r.rate+"*N"

        # create new reaction to model the growth
        newreac =  Reaction(["N"], ["N","N"], "context")
        Reaction.param_counter += 1
        par1 = 'p'+str(Reaction.param_counter)
        Reaction.param_counter += 1
        par2 = 'p'+str(Reaction.param_counter)
        newreac.rate = par1 + "*( 1 - N/" + par2 + ")"
        newreac.param.append(par1)
        newreac.param.append(par2)

        # add new reaction to list
        self.reactions.append( newreac )
        self.nreactions += 1
       
        # add new species to list
        self.species.append("N")
        self.nspecies += 1

        # add new parameters to list
        for k in newreac.param:
            self.parameters.append( k )

        # recalculate stoichiometry matrix
        self.stoichiometry_matrix = stoichiometry(self.nspecies, self.nreactions, self.species, self.reactions)
        
