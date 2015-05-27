from synbioweaver.core import *
from synbioweaver.aspects.modelDefinitions import *
import numpy, os, copy

class PrintReactionNetwork(Aspect):
    
    def mainAspect(self):
        self.addWeaverOutput(self.printReactionNetwork)

    def printReactionNetwork(self,weaverOutput):
        # We are expecting either a set of reactions or a context

        if getattr(weaverOutput, "getContext", None) != None:
            self.model = weaverOutput.getContext()
        else:
            if getattr(weaverOutput, "getReactions", None) != None:
                self.model = weaverOutput.getReactions()
            else:
                print "printReactionNetwork : Neither getContext() or getReactions() is available. Quitting"
                exit()

        self.printReactions()
        self.printStoichiometryMatrix()

    def printReactions(self):
        print "\n\nGenerated set of reactions:"
        for i in range(self.model.nreactions):
            print self.model.reactions[i].reactionString().rjust(45," "), "\t\t\t", self.model.reactions[i].rate.rjust(35," "), \
                "\t\t\t", self.model.reactions[i].paramString().rjust(10," "), "\t\t\t", self.model.reactions[i].process.rjust(15," ")
        
    def printStoichiometryMatrix(self):
        print "\n\nGenerated stoichiometry matrix:"
        print "".rjust(45," "),
        for j in range(self.model.nspecies):
            print str(self.model.species[j]).rjust(15," "),
        print ""

        for i in range(self.model.nreactions):
            print self.model.reactions[i].reactionString().rjust(45," "),
            for j in range(self.model.nspecies):
                print repr(self.model.stoichiometry_matrix[i,j]).rjust(15," "),
            print ""
