from synbioweaver.core import *
from synbioweaver.aspects.modelDefinitions import *
import numpy, os
from copy import *

def sort_react_by_param( x ):
    tosort = x.param[0]
    return tosort

class ContextGenerator(Aspect):
    def __init__(self, context_models):
        super(ContextGenerator,self).__init__()
        self.context_models = context_models

    def mainAspect(self):
        ContextGenerator.builtReactions = False
        self.addWeaverOutput(self.getContext)

    def getContext(self,weaverOutput):

        if ContextGenerator.builtReactions == False:
            # first access the existing reactions
            if getattr(weaverOutput, "getReactions", None) != None:
                self.originalModel = weaverOutput.getReactions()
            else:
                print "ContextGenerator : getReactions() is not available. Quitting"
                exit()

            self.models = []
            for p in self.context_models:
                self.models.append( self.createModel(p) )

            ContextGenerator.builtReactions = True

        return deepcopy(self.models)

    
    def createModel(self, mgen):
        # create a new model based upon the original

        newModel = deepcopy(self.originalModel)

        # reset the global reaction counter
        Reaction.param_counter = len(set(self.originalModel.parameters))

        # order the reactions based on parameters so that the dual system is assigned the same parameters to reactions
        self.originalModel.reactions = sorted( self.originalModel.reactions, key=sort_react_by_param )

        if mgen == "full-context-dep":
            # create a supermodel with double the number of parameters
            newreacs = []
            for r in self.originalModel.reactions: 
                #newreacs.append(Reaction( ['c'+x for x in r.reactants], ['c'+x for x in r.products], "context"))
                nreact = [ Species("con", 'cc'+x.name) for x in r.reactants ]
                nprods = [ Species("con",'cc'+x.name) for x in r.products ]
                newreacs.append(Reaction( nreact, nprods, "context")) 

                #newreacs.append(Reaction( r.reactants, r.products, "context")) 
                newreacs[-1].assignMassAction()
                #print r.reactionString(), r.param, newreacs[-1].reactionString(), newreacs[-1].param

            for r in newreacs:
                newModel.reactions.append( r ) 
                newModel.nreactions += 1 

                for k in r.param:
                    newModel.parameters.append( k )
                
            newspecies = []
            for s in self.originalModel.species:
                newspecies.append( Species("con",'cc'+s.name) )

            for s in newspecies:
                newModel.species.append(s)
                newModel.nspecies += 1

        elif mgen == "context-ind":
            # create a supermodel with the same number of parameters
            newreacs = []
            for r in self.originalModel.reactions:
                nreact = [ Species("con", 'cc'+x.name) for x in r.reactants ]
                nprods = [ Species("con",'cc'+x.name) for x in r.products ]
                newreacs.append(Reaction( nreact, nprods, "context"))
                  
                #print r.reactionString(), r.paramString(), r.param, r.rate
                newreacs[-1].assignMassActionExisting(r.param[0])

            for r in newreacs:
                newModel.reactions.append( r )
                newModel.nreactions += 1

                for k in r.param:
                    newModel.parameters.append( k )

            newspecies = []
            for s in self.originalModel.species:
                newspecies.append( Species("con",'cc'+s.name) )

            for s in newspecies:
                newModel.species.append(s)
                newModel.nspecies += 1

        elif mgen in ["rnaTransc","rnaDeg","proteinTransl","proteinDeg"]:
            # create a supermodel with different parameters for rnaTransc
            newreacs = []
            for r in self.originalModel.reactions:
                nreact = [ Species("con", 'cc'+x.name) for x in r.reactants ]
                nprods = [ Species("con",'cc'+x.name) for x in r.products ]
                newreacs.append(Reaction( nreact, nprods, "context"))

                if r.process != mgen:
                    newreacs[-1].assignMassActionExisting(r.param[0])
                else:
                    print Reaction.param_counter
                    newreacs[-1].assignMassAction()
                    #freeparams = list(set(self.))
                    
                    
            for r in newreacs:
                newModel.reactions.append( r )
                newModel.nreactions += 1

                for k in r.param:
                    newModel.parameters.append( k )

            newspecies = []
            for s in self.originalModel.species:
                newspecies.append( Species("con",'cc'+s.name) )

            for s in newspecies:
                newModel.species.append(s)
                newModel.nspecies += 1

        else:
            print "ContextGenerator : process unknown :", mgen, "Quitting"
            exit()

        # recalculate stoichiometry matrix
        newModel.stoichiometry_matrix = newModel.stoichiometry(newModel.nspecies, newModel.nreactions, newModel.species, newModel.reactions)

        # calculate unique parameters
        freeparams = list(set(newModel.parameters))
        newModel.freeparams = freeparams

        return newModel
