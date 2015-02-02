from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
from synbioweaver.aspects.promoterMappingAspect import *
import numpy, os, copy

class PostTranslationalCoupling(Aspect):
    
    def mainAspect(self):
        PostTranslationalCoupling.builtReactions = False
        self.addWeaverOutput(self.getContext)

    def getContext(self,weaverOutput):
        
        # first access the existing reactions
        if getattr(weaverOutput, "getReactions", None) != None:
            self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getReactions()
        else:
            print "PostTranslationalCoupling : getReactions() is available. Quitting"
            exit()

        if PostTranslationalCoupling.builtReactions == False:
            self.addPostTranslationalCoupling()

            PostTranslationalCoupling.builtReactions = True

        return [self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters]

    
    def addPostTranslationalCoupling(self):
        # modify the reactions to add enzymatic decay of tagged proteins
       
        tagged_proteins = []
        for sp in self.species:
            tag_marker = "_tag"
            mRNA_marker = "m"

            if tag_marker in sp and mRNA_marker not in sp:
                print "Found tag:", sp
                tagged_proteins.append( sp )

        # Create new reactions
        rates, params = assignEnzDeg(tagged_proteins)
           
        for i in range(len(tagged_proteins)):

            newreac =  Reaction([ tagged_proteins[i] ], [], "context")
            print rates[i]
            newreac.rate = rates[i]
            newreac.param = params[i]

            for k in newreac.param:
                self.parameters.append( k )
            self.reactions.append( newreac )

        # update number of species and reactions
        self.nreactions = len(self.reactions)
        self.nspecies = len(self.species)

        # recalculate stoichiometry matrix
        self.stoichiometry_matrix = stoichiometry(self.nspecies, self.nreactions, self.species, self.reactions)

        #for r in self.reactions:
        #    print r.rate

# assign enzymatic degradation
def assignEnzDeg(tagged_proteins):
    totalregs = len( tagged_proteins )
        
    #Reaction.param_counter += 1
    #par = "p"+str(Reaction.param_counter)

    #self.param.append( par )

    params = [[] for i in range(totalregs)]
    num = "( ";
    denom = "/( 1"

    terms = []
   
    for i in range(0,totalregs):
        Reaction.param_counter += 1
        par = "p"+str(Reaction.param_counter)
        terms.append( par + "*" + tagged_proteins[i] )
        params[i].append( par )
        
        denom = denom + " + " + par + "*" + tagged_proteins[i]
    denom = denom + " )"

    rates = []
    for i in range(0,totalregs):
        Reaction.param_counter += 1
        par = "p"+str(Reaction.param_counter)
        params[i].append( par )
        
        rates.append(  par + "* " + terms[i] + denom )

    #num = num + " )"
    #denom = denom + " )"
    
    #Reaction.param_counter += 1
    #par = "p"+str(Reaction.param_counter)
    #params.append( par )

    #rate = par + "* " +  num+denom + ")"

    return [rates, params]
