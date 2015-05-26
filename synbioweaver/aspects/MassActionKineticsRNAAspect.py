from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
from synbioweaver.aspects.molecularReactions import *
import numpy, os
from copy import *

class MassActionKineticsRNA(Aspect, MolecularReactions):
    
    def mainAspect(self):
        MassActionKineticsRNA.builtReactions = False
        self.addWeaverOutput(self.getReactions)

        self.reactions = []
        self.nreactions = 0
        self.species = []
        self.nspecies = 0
        self.stoichiometry_matrix = 0
        self.parameters = []

    def getReactions(self, weaverOutput):
        
        if getattr(weaverOutput, "buildPromoterMap", None) == None:
            sys.exit("MassActionKineticsRNA : buildPromoterMap() is unavailable. Quitting"  )

        if MassActionKineticsRNA.builtReactions == False:
            self.promoterMap = weaverOutput.buildPromoterMap()
          
            self.getReactionsMassActionRNA()
            self.getMolecules(weaverOutput)
          
            # Remove duplicate species
            self.species = list( set(self.species) ) 
            
            # Remove zero species
            if "zero" in self.species:
                self.species.remove("zero")

            # Finalise number of species and reactions
            self.nreactions = len(self.reactions)
            self.nspecies = len(self.species)
        
            # assign mass action rates and parameters
            for r in self.reactions:
                r.assignMassAction()

                for k in r.param:
                    self.parameters.append( k )

            # calculate stoichiometry
            self.stoichiometry_matrix = stoichiometry(self.nspecies, self.nreactions, self.species, self.reactions)
            
            MassActionKineticsRNA.builtReactions = True

        return [deepcopy(self.nspecies), deepcopy(self.nreactions), deepcopy(self.species), deepcopy(self.reactions), deepcopy(self.stoichiometry_matrix), deepcopy(self.parameters) ]

    
    def getReactionsMassActionRNA(self):
        for key in range( len(self.promoterMap) ):
            mapping = self.promoterMap[key]

            partname = mapping.getId()
            regulators = mapping.getRegulators()
            codings = mapping.getCodings()

            if len(regulators) == 0:
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
                self.reactions.append( Reaction([], prods1, "rnaTransc" ) )

                # translation
                for p in codings:
                    # translation
                    self.reactions.append( Reaction( ["m"+str(p)], [str(p)], "proteinTransl" ) )
                    # decay of mRNAs
                    self.reactions.append( Reaction(["m"+str(p)], [], "rnaDeg") )
                    # decay of proteins
                    self.reactions.append( Reaction([p], [], "proteinDeg") )

            else:
                for i in range(len(regulators)):
                    regulator = regulators[i]

                    complx = partname + '_' + str(regulator)
                
                    # the binding/unbinding reactions Pr + R <-> Pr_R
                    self.reactions.append( Reaction([partname, regulator], [complx], "dnaBind" ) )
                    self.reactions.append( Reaction([complx], [partname, regulator], "dnaUnbind") )
                    self.species.append( partname )
                    self.species.append( str(regulator) )
                    self.species.append( complx )

                    # if positive promoter then the bound complex expresses
                    if mapping.getPolarities()[i] == 1:
                        prods = deepcopy(codings)
                    
                        mprods = ["m"+str(x) for x in prods]
                        mprods.append(complx)

                        self.reactions.append( Reaction([complx], mprods, "rnaTransc") )
                        for p in codings:
                            self.reactions.append( Reaction([p], [], "proteinDeg") )
                            self.reactions.append( Reaction(["m"+str(p)], [], "rnaDeg") )
                            self.reactions.append( Reaction(["m"+str(p)], [p], "proteinTransl") )
                            self.species.append( p )
                            self.species.append( "m"+str(p) )

                    # if negative promoter then just the promoter expresses
                    else:
                        prods = deepcopy(codings)
                    
                        mprods = ["m"+str(x) for x in prods]
                        mprods.append(partname)
                        self.reactions.append( Reaction([partname], mprods,"rnaTransc") )

                        for p in codings:
                            self.reactions.append( Reaction([p], [], "proteinDeg") )
                            self.reactions.append( Reaction(["m"+str(p)], [], "rnaDeg") )
                            self.reactions.append( Reaction(["m"+str(p)], [p], "proteinTransl") )
                            self.species.append( p )
                            self.species.append( "m"+str(p) )

