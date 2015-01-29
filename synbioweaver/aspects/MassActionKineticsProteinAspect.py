from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os, copy, sys

class MassActionKineticsProtein(Aspect):
    
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
            
            self.promoterMap, self.circuitMap = weaverOutput.buildPromoterMap()
            #self.locatedPromoters = weaverOutput.getLocatedParts()

            self.getReactionsMassActionProtein()
            self.getMolecules(weaverOutput)
        
            # Remove duplicate species
            self.species = list( set(self.species) ) 
            
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


    def fullRegulatorName(self, regulator, prmtr_name ):
        return str(regulator)
        #return self.circuitMap[prmtr_name] + "." + str(name)
    
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
                complx = partname + '_' + self.fullRegulatorName(regulator, partname)
                
                # the binding/unbinding reactions Pr + R <-> Pr_R
                self.reactions.append( Reaction([partname, regulator], [complx], "dnaBind") )
                self.reactions.append( Reaction([complx], [partname, self.fullRegulatorName(regulator, partname) ], "dnaUnbind") )
                self.species.append( partname )
                self.species.append( self.fullRegulatorName(regulator, partname) )
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


    def getMolecules(self, weaverOutput):

        for mol in weaverOutput.moleculeList:
            # Add all molecules to list then do a unique
            #for i in range( len(mol.before) ):
            #    print "scope:", str(mol), mol.before[i].scope
            #
            #for i in range( len(mol.after) ):
            #    print "scope:", str(mol), mol.after[i].scope

            print "scope:", str(mol), mol.before[0].scope.circuitName

            circuitName = mol.before[0].scope.circuitName
            spl = circuitName + "." + str(mol).split("(")[0]
            self.species.append(spl)
            print "newname:", spl 

            #print mol, mol.before, mol.after

            # up until this point we have captured everything other than additional molecule-molecule reactions
            for j in range(len(mol.before)):
                # this indicates that downstream is not a part but molecules
                if type(mol.before[j]) == type(list()):
                    #print mol.before[j][0], mol.before[j][1], mol
                    self.reactions.append( Reaction([mol.before[j][0], mol.before[j][1]], [mol], "complexAss") )
                    self.reactions.append( Reaction([mol], [mol.before[j][0], mol.before[j][1]], "complexDiss") )
                    self.reactions.append( Reaction([mol], [], "complexDeg" ) )
                    
