from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
from synbioweaver.aspects.molecularReactions import *
from synbioweaver.aspects.promoterMappingAspect import *
import numpy, os, copy, sys

class MassActionKineticsProtein(Aspect, MolecularReactions):
    
    def mainAspect(self):
        MassActionKineticsProtein.builtReactions = False
        self.addWeaverOutput(self.getReactions)

        self.reactions = []
        self.nreactions = 0
        self.species = []
        self.nspecies = 0
        self.stoichiometry_matrix = 0
        self.parameters = []

    def getReactions(self, weaverOutput):
        
        if getattr(weaverOutput, "buildPromoterMap", None) == None:
            sys.exit("MassActionKineticsProtein : buildPromoterMap() is unavailable. Quitting"  )

        if MassActionKineticsProtein.builtReactions == False:
            
            self.promoterMap = weaverOutput.buildPromoterMap()
            #self.locatedPromoters = weaverOutput.getLocatedParts()

            #print "Multiple comps :", PromoterMapping.multiComp

            self.getReactionsMassActionProtein()
            self.getMolecules(weaverOutput)
          
            # Remove duplicate species
            self.species = list( set(self.species) ) 
            
            # Remove zero species
            if "zero" in self.species:
                self.species.remove("zero")

            # Finalise number of species, reactions, parameters
            self.nreactions = len(self.reactions)
            self.nspecies = len(self.species)
        
            # assign mass action rates and parameters
            for r in self.reactions:
                r.assignMassAction()
                
                for k in r.param:
                    self.parameters.append( k )

            # calculate stoichiometry
            self.stoichiometry_matrix = stoichiometry(self.nspecies, self.nreactions, self.species, self.reactions)
       
            MassActionKineticsProtein.builtReactions = True

        return [self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters]

    def getReactionsMassActionProtein(self):
        for key in range( len(self.promoterMap) ):
            mapping = self.promoterMap[key]

            partname = mapping.getId()
            regulators = mapping.getRegulators()
            codings = mapping.getCodings()
           
            if len(regulators) == 0:
                # This is a const promoter
                # Add the promoter itself, no need for complexes
                prods = copy.deepcopy(codings)
                self.reactions.append( Reaction([], prods, "proteinExp") )

                for p in codings:
                    self.reactions.append( Reaction([p], [], "proteinDeg") )
                    self.species.append( p )
                    
            else:
                for i in range(len(regulators)):
                    regulator = regulators[i]
                    
                    complx = partname + '_' + str(regulator)

                    # the binding/unbinding reactions Pr + R <-> Pr_R
                    self.reactions.append( Reaction([partname, regulator], [complx], "dnaBind") )
                    self.reactions.append( Reaction([complx], [partname, str(regulator) ], "dnaUnbind") )
                    self.species.append( partname )
                    self.species.append( str(regulator) )
                    self.species.append( complx )

                    # if positive promoter then the bound complex expresses
                    if mapping.getPolarities()[i] == 1:
                        prods = copy.deepcopy(codings)
                        prods.append(complx)
                        
                        self.reactions.append( Reaction([complx], prods, "proteinExp") )
                        for p in codings:
                            self.reactions.append( Reaction([p], [], "proteinDeg" ) )
                            self.species.append( p )

                    # if negative promoter then just the promoter expresses
                    else:
                        #print "here", codings
                        prods = copy.deepcopy(codings)
                        prods.append(partname)
                        self.reactions.append( Reaction([partname], prods, "proteinExp") )

                        for p in codings:
                            self.reactions.append( Reaction([p], [], "proteinDeg") )
                            self.species.append( p )


  
