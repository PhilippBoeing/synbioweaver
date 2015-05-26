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
        self.nreactions = 0
        self.species = []
        self.nspecies = 0
        self.stoichiometry_matrix = 0
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
            
            # Remove duplicate species
            self.species = list( set(self.species) ) 
            
            # Remove zero species
            if "zero" in self.species:
                self.species.remove("zero")

            # Finalise number of species and reactions
            self.nreactions = len(self.reactions)
            self.nspecies = len(self.species)
        
            # assign mass action rates and parameters to certain processes
            for r in self.reactions:
                if r.process == "rnaDeg" or r.process == "proteinDeg" or r.process == "proteinTransl" or r.process == "complexDiss"\
                   or r.process == "complexAss" or r.process == "complexDeg":
                    r.assignMassAction()

                for k in r.param:
                    self.parameters.append( k )

            # calculate stoichiometry
            self.stoichiometry_matrix = stoichiometry(self.nspecies, self.nreactions, self.species, self.reactions)
       
            SheaAckersKineticsRNA.builtReactions = True

        return [deepcopy(self.nspecies), deepcopy(self.nreactions), deepcopy(self.species), deepcopy(self.reactions), deepcopy(self.stoichiometry_matrix), deepcopy(self.parameters)]

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

                prods = deepcopy(codings)

                # add species first
                #self.species.append( partname )
                for p in codings:
                    self.species.append( p )
                    self.species.append( "m"+str(p) )

                # production of mRNAs
                prods1 = ["m"+str(x) for x in prods]
                #prods1.append( partname )
                #self.reactions.append( Reaction([partname], prods1, "rnaTransc" ) )
                R1 = Reaction([], prods1, "rnaTransc" )
                if located == True:
                    R1.rate = self.locatedParts[partname].transferFunction()
                else:
                    R1.assignMassAction()
                
                
                self.reactions.append(R1)

                # translation
                for p in codings:
                    # translation
                    self.reactions.append( Reaction( ["m"+str(p)], [str(p)], "proteinTransl" ) )
                    # decay of mRNAs
                    self.reactions.append( Reaction(["m"+str(p)], [], "rnaDeg") )
                    # decay of proteins
                    self.reactions.append( Reaction([p], [], "proteinDeg") )

            else:

                # assign the transcription reaction
                prods = deepcopy(codings)
                mprods = ["m"+str(x) for x in prods]
                Rprod = Reaction([], mprods, "rnaTransc")
                Rprod.assignSA(mapping)
                self.reactions.append(Rprod)
                
                for i in range(len(regulators)):
                    regulator = regulators[i]

                    complx = partname + '_' + str(regulator)
                    for p in codings:
                        self.reactions.append( Reaction([p], [], "proteinDeg") )
                        self.reactions.append( Reaction(["m"+str(p)], [], "rnaDeg") )
                        self.reactions.append( Reaction(["m"+str(p)], [p], "proteinTransl") )
                        self.species.append( p )
                        self.species.append( "m"+str(p) )

