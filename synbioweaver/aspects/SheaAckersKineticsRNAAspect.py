from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
from synbioweaver.aspects.molecularReactions import *
import numpy, os
from copy import *

class SheaAckersKineticsRNA(Aspect, MolecularReactions):
    
    def mainAspect(self):
        SheaAckersKineticsRNA.builtReactions = False
        self.addWeaverOutput(self.getReactions)

        self.reactions = []
        self.species = []
        self.parameters = []
        self.locatedParts = {}

    def getReactions(self, weaverOutput):

        if getattr(weaverOutput, "buildPromoterMap", None) == None:
            sys.exit("SheaAckersKineticsRNA : buildPromoterMap() is unavailable. Quitting"  )

        if SheaAckersKineticsRNA.builtReactions == False:
            self.promoterMap = weaverOutput.buildPromoterMap()

            if getattr(weaverOutput, "getLocatedParts", None) != None:
                self.locatedParts = weaverOutput.getLocatedParts()
                print "located:", self.locatedParts 
            else:
                self.locatedParts = []

            self.getReactionsSheaAckersRNA()
            self.getMolecules(weaverOutput)
            
            # assign mass action rates and parameters to certain processes
            for r in self.reactions:
                if r.process == "rnaDeg" or r.process == "proteinDeg" or r.process == "proteinTransl" or r.process == "complexDiss"\
                   or r.process == "complexAss" or r.process == "complexDeg":
                    r.assignMassAction()
                for k in r.param:
                    self.parameters.append( k )

            self.newModel = Model( self.species, self.reactions, self.parameters )

            SheaAckersKineticsRNA.builtReactions = True

        return deepcopy(self.newModel)

    def getReactionsSheaAckersRNA(self):
        for key in range( len(self.promoterMap) ):
            mapping = self.promoterMap[key]

            partname = mapping.getId()
            regulators = mapping.getRegulators()
            codings = mapping.getCodings()

            located = False
            if partname in self.locatedParts:
                located = True

            if len(regulators) == 0 or located == True:
                # This is a const promoter
                # need to add:
                # pr -> mX + pr, mX -> X, mX -> 0, X -> 0 

                for p in codings:
                    self.species.append( Species(mapping.getScope(),str(p), "protein") )
                    self.species.append( Species(mapping.getScope(),"m"+str(p), "mRNA") )

                # production of mRNAs
                prods1 = [ Species(mapping.getScope(),"m"+str(x)) for x in codings]
                R1 = Reaction([], prods1, "rnaTransc" )
                #if located == True:
                #    R1.rate = self.locatedParts[partname].transferFunction()
                #else:
                #    R1.assignMassAction()
                self.reactions.append(R1)

                # translation
                for p in codings:
                    # translation
                    self.reactions.append( Reaction( [Species(mapping.getScope(),"m"+str(p))], [Species(mapping.getScope(),str(p))], "proteinTransl" ) )
                    # decay of mRNAs
                    self.reactions.append( Reaction([Species(mapping.getScope(),"m"+str(p))], [], "rnaDeg") )
                    # decay of proteins
                    self.reactions.append( Reaction([Species(mapping.getScope(),str(p))], [], "proteinDeg") )

            else:

                # assign the transcription reaction
                mprods = [Species(mapping.getScope(),"m"+str(x)) for x in codings]
                Rprod = Reaction([], mprods, "rnaTransc")
                Rprod.assignSA(mapping)
                self.reactions.append(Rprod)
                
                for i in range(len(regulators)):
                    regulator = regulators[i]

                    complx = partname + '_' + str(regulator)
                    for p in codings:
                        self.reactions.append( Reaction([Species(mapping.getScope(),str(p))], [], "proteinDeg") )
                        self.reactions.append( Reaction([Species(mapping.getScope(),"m"+str(p))], [], "rnaDeg") )
                        self.reactions.append( Reaction([Species(mapping.getScope(),"m"+str(p))], [Species(mapping.getScope(),str(p))], "proteinTransl") )
                        self.species.append( Species(mapping.getScope(),str(p), "protein" ) )
                        self.species.append( Species(mapping.getScope(),"m"+str(p), "mRNA" ) )

