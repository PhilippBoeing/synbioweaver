from synbioweaver.core import *
from synbioweaver.aspects.modelDefinitions import *
import numpy, os
from copy import *

class LogisticGrowth(Aspect):
    
    def mainAspect(self):
        LogisticGrowth.builtReactions = False
        self.addWeaverOutput(self.getContext)

    def getContext(self,weaverOutput):

        if LogisticGrowth.builtReactions == False:
            
            # first access the existing reactions
            if getattr(weaverOutput, "getReactions", None) != None:
                self.model = weaverOutput.getReactions()
            else:
                print "LogisticGrowth : getReactions() is available. Quitting"
                exit()

            self.addLogisticGrowth()

            LogisticGrowth.builtReactions = True

        return deepcopy(self.model)

    
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
        for r in self.model.reactions:
            if r.process == "proteinExp" or r.process == "rnaTransc":
                r.products.append( Species("con","N") )
                r.reactants.append( Species("con","N") )
                r.rate = r.rate+"*N"

        # create new reaction to model the growth
        newreac =  Reaction([Species("con","N")], [Species("con","N"),Species("con","N")], "context")
        Reaction.param_counter += 1
        par1 = 'p'+str(Reaction.param_counter)
        Reaction.param_counter += 1
        par2 = 'p'+str(Reaction.param_counter)
        newreac.rate = par1 + "*N*( 1 - N/" + par2 + ")"
        newreac.param.append(par1)
        newreac.param.append(par2)

        # add new reaction to list
        self.model.reactions.append( newreac )
        self.model.nreactions += 1
       
        # add new species to list
        self.model.species.append( Species("con","N") )
        self.model.nspecies += 1

        # add new parameters to list
        for k in newreac.param:
            self.model.parameters.append( k )

        # recalculate stoichiometry matrix
        self.model.stoichiometry_matrix = self.model.stoichiometry(self.model.nspecies, self.model.nreactions, self.model.species, self.model.reactions)
        
