from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os, copy

class MassActionKineticsRNA(Aspect):
    
    def mainAspect(self):
        self.addWeaverOutput(self.getReactions)

        self.reactions = []
        self.nreactions = 0
        self.species = []
        self.nspecies = 0
        self.stoichiometry_matrix = 0

    def getReactions(self, weaverOutput):
        self.promoterMap = weaverOutput.buildPromoterMap()
        self.getReactionsMassActionRNA()

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
                    self.reactions.append( Reaction([mol.before[j][0], mol.before[j][1]], [mol], "complexAss" ) )
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

    
    def getReactionsMassActionRNA(self):
        for key in range( len(self.promoterMap) ):
            mapping = self.promoterMap[key]
            partname = str( mapping.prmtr_name )
            regulator = mapping.regulators[0]
            codings = [str(x) for x in mapping.coding ]
            #print partname, regulator, codings

            if regulator == None:
                # This is a const promoter
                # need to add:
                # pr -> mX + pr, mX -> X, mX -> 0, X -> 0 

                prods = copy.deepcopy(codings)

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
                # This is a regulated promoter
                complx = partname + '_' + str(regulator)
                
                # the binding/unbinding reactions Pr + R <-> Pr_R
                self.reactions.append( Reaction([partname, regulator], [complx], "dnaBind" ) )
                self.reactions.append( Reaction([complx], [partname, regulator], "dnaUnbind") )
                self.species.append( partname )
                self.species.append( str(regulator) )
                self.species.append( complx )

                # if positive promoter then the bound complex expresses
                if mapping.regulators[0] == 1:
                    prods = copy.deepcopy(codings)
                    
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
                    prods = copy.deepcopy(codings)
                    
                    mprods = ["m"+str(x) for x in prods]
                    mprods.append(partname)
                    self.reactions.append( Reaction([partname], mprods,"rnaTransc") )

                    for p in codings:
                        self.reactions.append( Reaction([p], [], "proteinDeg") )
                        self.reactions.append( Reaction(["m"+str(p)], [], "rnaDeg") )
                        self.reactions.append( Reaction(["m"+str(p)], [p], "proteinTransl") )
                        self.species.append( p )
                        self.species.append( "m"+str(p) )
