from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os, copy

class PrintReactionNetwork(Aspect):
    
    def mainAspect(self):
        self.addWeaverOutput(self.printReactionNetwork)

    def printReactionNetwork(self,weaverOutput):
        # We are expecting either a set of reactions or a context

        if getattr(weaverOutput, "getContext", None) != None:
            self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix = weaverOutput.getContext()
        else:
            if getattr(weaverOutput, "getReactions", None) != None:
                self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix = weaverOutput.getReactions()
            else:
                print "printReactionNetwork : Neither getContext() or getReactions() is available. Quitting"
                exit()

        self.printReactions()
        self.printStoichiometryMatrix()

    def printReactions(self):
        print "\n\nGenerated set of reactions:"
        for i in range(self.nreactions):
            print self.reactions[i].reactionString().rjust(45," "), "\t\t\t", self.reactions[i].rate.rjust(35," "), \
                "\t\t\t", self.reactions[i].param.rjust(10," "), "\t\t\t", self.reactions[i].process.rjust(15," ")
        
    def printStoichiometryMatrix(self):
        print "\n\nGenerated stoichiometry matrix:"
        print "".rjust(45," "),
        for j in range(self.nspecies):
            print repr(self.species[j]).rjust(15," "),
        print ""

        for i in range(self.nreactions):
            print self.reactions[i].reactionString().rjust(45," "),
            for j in range(self.nspecies):
                print repr(self.stoichiometry_matrix[i,j]).rjust(15," "),
            print ""
