from synbioweaver.core import *
from synbioweaver.aspects.modelDefinitions import *
import numpy, os
from copy import *

class LagLogisticGrowth(Aspect):
    
    def mainAspect(self):
        LagLogisticGrowth.builtReactions = False
        self.addWeaverOutput(self.getContext)

    def getContext(self,weaverOutput):

        if LagLogisticGrowth.builtReactions == False:
        
            # first access the existing reactions
            if getattr(weaverOutput, "getReactions", None) != None:
                self.model = weaverOutput.getReactions()
            else:
                print "LagLogisticGrowth : getReactions() is not available. Quitting"
                exit()

            self.addLagLogisticGrowth()

            LagLogisticGrowth.builtReactions = True

        return deepcopy(self.model)

    
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

        SN = Species("con","N")
        SNd = Species("con","Nd")

        # modify exixting reactions : change products, reactants and rate - no new parameters
        for r in self.model.reactions:
            if r.process == "proteinExp" or r.process == "rnaTransc":
                r.products.append( SN )
                r.reactants.append( SN )
                r.rate = r.rate+"*N"

        ###
        # create new logistic growth reaction to model the growth
        newreac =  Reaction([SN], [SN,SN], "context")
        Reaction.param_counter += 1
        par1 = 'p'+str(Reaction.param_counter)
        Reaction.param_counter += 1
        par2 = 'p'+str(Reaction.param_counter)
        newreac.rate = par1 + "*N*( 1 - N/" + par2 + ")"
        newreac.param.append(par1)
        newreac.param.append(par2)

        # add new reaction, species, params to list
        self.model.reactions.append( newreac )
        self.model.species.append(SN)
        for k in newreac.param:
            self.model.parameters.append( k )

        ###
        # add Nd species
        self.model.species.append(SNd)

        # add decay reaction mass action
        newreac2 =  Reaction([SNd], [], "context")
        newreac2.assignMassAction()
        for k in newreac2.param:
            self.model.parameters.append( k )
        self.model.reactions.append( newreac2 )

        # add transition reaction mass action
        newreac3 =  Reaction([SNd], [SN], "context")
        newreac3.assignMassAction()
        for k in newreac3.param:
            self.model.parameters.append( k )
        self.model.reactions.append( newreac3 )

        # update number of species and reactions
        self.model.nreactions = len(self.model.reactions)
        
        self.model.nspecies = len(self.model.species)
      
        # recalculate stoichiometry matrix
        self.model.stoichiometry_matrix = self.model.stoichiometry(self.model.nspecies, self.model.nreactions, self.model.species, self.model.reactions)
        
