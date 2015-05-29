from synbioweaver.core import *
from synbioweaver.aspects.modelDefinitions import *
from synbioweaver.aspects.molecularReactions import *
from synbioweaver.aspects.promoterMappingAspect import *
import numpy, os, sys
from copy import *

class MassActionKineticsProtein(Aspect, MolecularReactions):
    
    def mainAspect(self):
        MassActionKineticsProtein.builtReactions = False
        self.addWeaverOutput(self.getReactions)

        self.reactions = []
        self.species = []
        self.parameters = []

    def getReactions(self, weaverOutput):
        
        if getattr(weaverOutput, "buildPromoterMap", None) == None:
            sys.exit("MassActionKineticsProtein : buildPromoterMap() is unavailable. Quitting"  )

        if MassActionKineticsProtein.builtReactions == False:
            
            self.promoterMap = weaverOutput.buildPromoterMap()
            #self.locatedPromoters = weaverOutput.getLocatedParts()

            # get molecular reactions first
            mmol = MolecularReactions()
            mreactions, mspecies = mmol.getMolecules(weaverOutput)

             # get the unique namespaces and generate reactions
            nspaces = getNamespaces(self.promoterMap)
            for ns in nspaces:
                reactions = self.getReactionsMassActionProtein(ns)
                
                # assign mass action rates and parameters to all processes
                for r in reactions:
                    r.assignMassAction()
                    for k in r.param:
                        self.parameters.append( k )
                
                self.reactions = self.reactions + reactions

                # assign mass action rates to molecular reactions
                for r in mreactions:
                    # get either reactant or product
                    if len(r.reactants) != 0:
                        scope = r.reactants[0].scope
                    else:
                        scope = r.products[0].scope
        
                    if scope == ns:
                        r.assignMassAction()
                        for k in r.param:
                            self.parameters.append( k )
                        
            
            self.reactions = self.reactions + mreactions 
            self.species = self.species + mspecies

            self.newModel = Model( self.species, self.reactions, self.parameters )
           
            MassActionKineticsProtein.builtReactions = True

        return deepcopy( self.newModel )

    def getReactionsMassActionProtein(self, nspace):
        reactions = []

        for key in range( len(self.promoterMap) ):
            mapping = self.promoterMap[key]

            if mapping.getScope() != nspace:
                continue

            partname = mapping.getId()
            regulators = mapping.getRegulators()
            codings = mapping.getCodings()
            #print partname, regulators, codings
           
            if len(regulators) == 0:
                # This is a const promoter
                # Add the promoter itself, no need for complexes

                # create a list of species for the products
                #prods = deepcopy(codings)
                prods = [ Species(mapping.getScope(), x) for x in codings]

                reactions.append( Reaction([], prods, "proteinExp") )

                for p in codings:
                    reactions.append( Reaction([ Species(mapping.getScope(),p) ], [], "proteinDeg") )
                    self.species.append( Species(mapping.getScope(),p, "protein") )
                    
            else:
                for i in range(len(regulators)):
                    regulator = regulators[i]
                    
                    complx = partname + '_' + str(regulator)

                    # the binding/unbinding reactions Pr + R <-> Pr_R
                    newSpecies1 = Species(mapping.getScope(),partname, "promoter")
                    newSpecies2 = Species(mapping.getScope(),regulator, "protein" )
                    newSpecies3 = Species(mapping.getScope(),complx, "bindingDNA")
                    
                    reactions.append( Reaction([newSpecies1, newSpecies2], [newSpecies3], "dnaBind") )
                    reactions.append( Reaction([newSpecies3], [newSpecies1, newSpecies2], "dnaUnbind") )
                    self.species.append( newSpecies1 )
                    self.species.append( newSpecies2 )
                    self.species.append( newSpecies3 )

                    # if positive promoter then the bound complex expresses
                    if mapping.getPolarities()[i] == 1:
                        #prods = deepcopy(codings)
                        prods = [ Species(mapping.getScope(), x) for x in codings]
                        prods.append( newSpecies3 )
                        
                        reactions.append( Reaction([newSpecies3], prods, "proteinExp") )
                        for p in codings:
                            reactions.append( Reaction([Species(mapping.getScope(),p)], [], "proteinDeg" ) )
                            self.species.append( Species(mapping.getScope(),p, "protein") )

                    # if negative promoter then just the promoter expresses
                    else:
                        #prods = deepcopy(codings)
                        prods = [ Species(mapping.getScope(), x) for x in codings]
                        prods.append(newSpecies1)

                        reactions.append( Reaction([newSpecies1], prods, "proteinExp") )
                        for p in codings:
                            reactions.append( Reaction([Species(mapping.getScope(),p)], [], "proteinDeg") )
                            self.species.append( Species(mapping.getScope(),p, "protein") )
  
        return reactions
