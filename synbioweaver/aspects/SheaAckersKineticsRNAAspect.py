from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os, copy

class SheaAckersKineticsRNA(Aspect):
    
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

        return [self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters]

    def getReactionsSheaAckersRNA(self):
        for key in range( len(self.promoterMap) ):
            mapping = self.promoterMap[key]
            partname = str( mapping.prmtr_name )

            regulator = mapping.regulators[0]
            codings = [str(x) for x in mapping.coding ]
            #print partname, regulator, codings

            located = False
            if partname in self.locatedParts:
                located = True

            if regulator == None or located == True:
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
                # This is a regulated promoter
                complx = partname + '_' + str(regulator)
                
                # if positive promoter
                if mapping.regulators[0] == 1:
                    prods = copy.deepcopy(codings)
                    mprods = ["m"+str(x) for x in prods]
                    Rprod = Reaction([], mprods, "rnaTransc")
                    Rprod.assignSA(mapping)
                    self.reactions.append(Rprod)

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
                    Rprod = Reaction([], mprods, "rnaTransc" )
                    Rprod.assignSA(mapping)
                    self.reactions.append(Rprod)

                    for p in codings:
                        self.reactions.append( Reaction([p], [], "proteinDeg") )
                        self.reactions.append( Reaction(["m"+str(p)], [], "rnaDeg") )
                        self.reactions.append( Reaction(["m"+str(p)], [p], "proteinTransl") )
                        self.species.append( p )
                        self.species.append( "m"+str(p) )


    def getMolecules(self, weaverOutput):

        for mol in weaverOutput.moleculeList:
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
                    

