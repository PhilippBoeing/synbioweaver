from synbioweaver.core import *
from synbioweaver.aspects.modelDefinitions import *
import numpy, os, copy

def sort_react_by_param( x ):
    tosort = x.param[0]
    return tosort

class PrintMultipleModels(Aspect):
    
    def mainAspect(self):
        self.addWeaverOutput(self.printMultipleModels)

    def printMultipleModels(self,weaverOutput):
        # We are expecting a context with multiple models

        if getattr(weaverOutput, "getContext", None) != None:
            self.models = weaverOutput.getContext()
        else:
            print "printMultipleModels : getContext() is missing. Quitting"
            exit()

        print "Multiple contextual models generated:", len(self.models) 
        for i,m in enumerate(self.models):
            print "Model", i
            print "species:", [str(x) for x in m.species]
            print "parameters:", len(m.freeparams), m.freeparams
            self.printReactions(m)

    def printReactions(self, model):
        #model.reactions = sorted( model.reactions, key=sort_react_by_param )

        print "\n\nGenerated set of reactions:"
        for i in range(model.nreactions):
            print model.reactions[i].reactionString().rjust(45," "), "\t\t\t", model.reactions[i].rate.rjust(35," "), \
                "\t\t\t", model.reactions[i].paramString().rjust(10," "), "\t\t\t", model.reactions[i].process.rjust(15," ")
        
    
