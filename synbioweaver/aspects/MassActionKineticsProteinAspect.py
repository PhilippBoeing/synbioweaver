from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os, copy

class MassActionKineticsProtein(Aspect):
    
    def mainAspect(self):
        self.addWeaverOutput(self.getReactions)

        self.reactions = []
        self.nreactions = 0
        self.species = []
        self.nspecies = 0
        self.stoichiometry_matrix = 0

    def getReactions(self, weaverOutput):
        self.promoterMap = weaverOutput.buildPromoterMap()
        self.locatedPromoters = weaverOutput.getLocatedParts()

        self.getReactionsMassActionProtein()

        for mol in weaverOutput.moleculeList:
            # This should be done using a molecule Type Advice
            #mol.generateReactions()

            # Add all molecules to list then do a unique
            spl = str(mol).split("(")[0]
            self.species.append(spl)
            
            #print mol, mol.before, mol.after

            # up until this point we have captured everything other than additional molecule-molecule reactions
            for j in range(len(mol.before)):
                # this indicates that downstream is not a part but molecules
                if type(mol.before[j]) == type(list()):
                    #print mol.before[j][0], mol.before[j][1], mol
                    self.reactions.append( Reaction([mol.before[j][0], mol.before[j][1]], [mol], "complexAss") )
                    self.reactions.append( Reaction([mol], [mol.before[j][0], mol.before[j][1]], "complexDiss") )
                    self.reactions.append( Reaction([mol], [], "complexDeg" ) )


        # Remove duplicate species
        self.species = list( set(self.species) ) 
            
        # Finalise number of species and reactions
        self.nreactions = len(self.reactions)
        self.nspecies = len(self.species)
        #print "species:", self.species, set(self.species)
        #print "nreactions/nspecies:", self.nreactions, self.nspecies

        # assign mass action rates and parameters
        for r in self.reactions:
            r.assignMassAction()

        # calculate stoichiometry
        self.stoichiometry_matrix = stoichiometry(self.nspecies, self.nreactions, self.species, self.reactions)
       
        return [self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix]

    
    def getReactionsMassActionProtein(self):
        for key in range( len(self.promoterMap) ):
            mapping = self.promoterMap[key]
            partname = str( mapping.prmtr_name )
            regulator = mapping.regulators[0]
            codings = [str(x) for x in mapping.coding ]
            #print partname, regulator, codings

            if regulator == None:
                # This is a const promoter
                # Add the promoter itself, no need for complexes
                prods = copy.deepcopy(codings)
                
                #prods.append( partname )
                #self.reactions.append( Reaction([partname], prods, "proteinExp") )
                #self.species.append( partname )
                self.reactions.append( Reaction([], prods, "proteinExp") )

                for p in codings:
                    self.reactions.append( Reaction([p], [], "proteinDeg") )
                    self.species.append( p )
            else:
                # This is a regulated promoter
                complx = partname + '_' + str(regulator)
                
                # the binding/unbinding reactions Pr + R <-> Pr_R
                self.reactions.append( Reaction([partname, regulator], [complx], "dnaBind") )
                self.reactions.append( Reaction([complx], [partname, regulator], "dnaUnbind") )
                self.species.append( partname )
                self.species.append( str(regulator) )
                self.species.append( complx )

                # if positive promoter then the bound complex expresses
                if mapping.regulators[0] == 1:
                    prods = copy.deepcopy(codings)
                    prods.append(complx)
                    
                    self.reactions.append( Reaction([complx], prods, "proteinExp") )
                    for p in codings:
                        self.reactions.append( Reaction([p], [], "proteinDeg" ) )
                        self.species.append( p )

                # if negative promoter then just the promoter expresses
                else:
                    prods = copy.deepcopy(codings)
                    prods.append(partname)
                    self.reactions.append( Reaction([partname], prods, "proteinExp") )

                    for p in codings:
                        self.reactions.append( Reaction([p], [], "proteinDeg") )
                        self.species.append( p )

