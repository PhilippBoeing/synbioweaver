from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os, copy

class LagLogisticGrowth(Aspect):
    
    def mainAspect(self):
        LagLogisticGrowth.builtReactions = False
        self.addWeaverOutput(self.getContext)

    def getContext(self,weaverOutput):
        
        # first access the existing reactions
        if getattr(weaverOutput, "getReactions", None) != None:
            self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getReactions()
        else:
            print "LagLogisticGrowth : getReactions() is available. Quitting"
            exit()

        if LagLogisticGrowth.builtReactions == False:
            self.addLagLogisticGrowth()

            LagLogisticGrowth.builtReactions = True

        return [self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters]

    
    def addLagLogisticGrowth(self):
        # modify the reactions to incorporate logistic growth
        # this essentially means converting any reactions of the form 
        #     0 -> P, rate k
        # to
        #     N -> N + P, rate kN
        # plus adding
        #     N -> 2N,    rate p1*N(1-N/p2)
        # here N are duplicating cells
        # we also add Nd cells
        # these can transform into N cells and also die
        #     Nd -> N
        #     Nd -> 0

        # modify exixting reactions : change products, reactants and rate - no new parameters
        for r in self.reactions:
            if r.process == "proteinExp" or r.process == "rnaTransc":
                r.products.append("N")
                r.reactants.append("N")
                r.rate = r.rate+"*N"

        ###
        # create new logistic growth reaction to model the growth
        newreac =  Reaction(["N"], ["N","N"], "context")
        Reaction.param_counter += 1
        par1 = 'p'+str(Reaction.param_counter)
        Reaction.param_counter += 1
        par2 = 'p'+str(Reaction.param_counter)
        newreac.rate = par1 + "*N*( 1 - N/" + par2 + ")"
        newreac.param.append(par1)
        newreac.param.append(par2)

        # add new reaction, species, params to list
        self.reactions.append( newreac )
        self.species.append("N")
        for k in newreac.param:
            self.parameters.append( k )

        ###
        # add Nd species
        self.species.append("Nd")

        # add decay reaction mass action
        newreac2 =  Reaction(["Nd"], [], "context")
        newreac2.assignMassAction()
        for k in newreac2.param:
            self.parameters.append( k )
        self.reactions.append( newreac2 )

        # add transition reaction mass action
        newreac3 =  Reaction(["Nd"], ["N"], "context")
        newreac3.assignMassAction()
        for k in newreac3.param:
            self.parameters.append( k )
        self.reactions.append( newreac3 )

        # update number of species and reactions
        self.nreactions = len(self.reactions)
        self.nspecies = len(self.species)

        # recalculate stoichiometry matrix
        self.stoichiometry_matrix = stoichiometry(self.nspecies, self.nreactions, self.species, self.reactions)
        
