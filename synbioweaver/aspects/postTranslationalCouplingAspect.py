from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
from synbioweaver.aspects.promoterMappingAspect import *
import numpy, os
from copy import *

class PostTranslationalCoupling(Aspect):
    
    def mainAspect(self):
        PostTranslationalCoupling.builtReactions = False
        self.reactions = []
        self.species = []
        self.parameters = []

        self.addWeaverOutput(self.getContext)

    def getContext(self,weaverOutput):
        
        if PostTranslationalCoupling.builtReactions == False:
            # first access the existing reactions
            if getattr(weaverOutput, "getReactions", None) != None:
                self.model = weaverOutput.getReactions()
            else:
                print "PostTranslationalCoupling : getReactions() is available. Quitting"
                exit()

            #print "running coupling"
            self.reactions = self.model.reactions
            self.species = self.model.species
            self.parameters = self.model.parameters
            
            self.addPostTranslationalCoupling()

            self.newModel = Model( self.species, self.reactions, self.parameters )

            PostTranslationalCoupling.builtReactions = True

        return deepcopy(self.newModel)
    
    def addPostTranslationalCoupling(self):
        # modify the reactions to add enzymatic decay of tagged proteins
       
        tagged_proteins = []
        for sp in self.species:
            tag_marker = "_tag"
            mRNA_marker = "m"

            if tag_marker in sp.name  and mRNA_marker not in sp.name:
                print "PostTranslationalCoupling : Found tag:", sp
                tagged_proteins.append( sp )

        # Create new reactions
        rates, params = assignEnzDeg(tagged_proteins)
           
        for i in range(len(tagged_proteins)):

            newreac =  Reaction([ tagged_proteins[i] ], [], "context")
            print "PostTranslationalCoupling : rates:", rates[i]
            newreac.rate = rates[i]
            newreac.param = params[i]

            for k in newreac.param:
                self.parameters.append( k )
            self.reactions.append( newreac )

        #

# assign enzymatic degradation
def assignEnzDeg(tagged_proteins):
    totalregs = len( tagged_proteins )
        
    #print "tagged_proteins:"
    #print tagged_proteins

    params = [[] for i in range(totalregs)]
    num = "( ";
    denom = "/( 1"

    terms = []
   
    for i in range(0,totalregs):
        Reaction.param_counter += 1
        par = "p"+str(Reaction.param_counter)
        terms.append( par + "*" + tagged_proteins[i].name )
        params[i].append( par )
        
        denom = denom + " + " + par + "*" + tagged_proteins[i].name
    denom = denom + " )"

    rates = []
    for i in range(0,totalregs):
        Reaction.param_counter += 1
        par = "p"+str(Reaction.param_counter)
        params[i].append( par )
        
        rates.append(  par + "* " + terms[i] + denom )

    return [rates, params]
