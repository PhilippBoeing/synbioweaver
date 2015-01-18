from synbioweaver.core import *
import numpy, os, copy

# Each reaction is assigned one of these processes
processes = [ "dnaBind", "dnaUnbind", "rnaTransc", "rnaDeg", "proteinTransl", "proteinExp", "proteinDeg", "complexAss", "complexDiss", "complexDeg" ]

class Reaction:
    param_counter = 0

    def __init__(self, reactants, products, process):
        
        # create a reaction string for output
        react = ' + '.join(map(str, reactants))
        prod = ' + '.join(map(str, products))
        self.reac_string = react + ' -> ' + prod
        #print self.reac_string

        # assign the products and reactants
        self.reactants = map(str, reactants)
        self.products = map(str, products)
        
        # assign a process
        if process in processes:
            self.process = process
        else:
            print "Error: process undefined:", process
            exit(1)

        self.param = ""
        self.rate = ""
            
    def assignMassAction(self):
        # assign a parameter number
        Reaction.param_counter += 1
        self.param = 'p'+str(Reaction.param_counter)

        # assign a rate
        if len(self.reactants) > 0:
            rt = '*'.join(map(str, self.reactants))
            self.rate = rt+'*'+self.param
        else:
            self.rate = self.param

    def assignSA(self, prmtrMapping):
        totalregs = len( prmtrMapping.regulators )
        totalpos = 0
        Reaction.param_counter += 1
        par = "p"+str(Reaction.param_counter)

        self.param = par

        num = "( " + par;
        denom = "/( 1 + " + par 
        for i in range(0,totalregs):
            Reaction.param_counter += 1
            par = "p"+str(Reaction.param_counter)
            self.param = self.param + ", " + par

            denom = denom + " + " + par + "*" + str( prmtrMapping.regulators[i] )
                
            if prmtrMapping.polarities[i] == +1:
                totalpos += 1
                num = num + " + " + par + "*" + str( prmtrMapping.regulators[i] )

        num = num + " )"
        denom = denom + " )"
    
        self.rate = num+denom
        

# This helper class associates promoters with downstream coding regions, regulators and polarities
class promoterMapping:

    def __init__(self, prmtr_name, regulators, polarities, coding):
        self.prmtr_name = prmtr_name
        self.regulators = regulators
        self.polarities = polarities
        self.coding = coding

class ReactionNetworks(Aspect):
    
    def mainAspect(self):
        Reaction.param_counter = 0
        ReactionNetworks.builtMap = False
        # The goal of this aspect is to parse the design and generate a set of reactions, rates and parameters
        self.addTypeAdvice(PartSignature('*.RBS+'), self.notNeeded, 'generatePromoterMap')
        self.addTypeAdvice(PartSignature('*.Terminator+'), self.notNeeded, 'generatePromoterMap')
        self.addTypeAdvice(PartSignature('*.CodingRegion+'), self.notNeeded, 'generatePromoterMap')
        self.addTypeAdvice(PartSignature('*.Promoter+'), self.promoterCoding, 'generatePromoterMap')

        self.addWeaverOutput(self.buildPromoterMap)
        
        self.promoterMap = []

    def promoterCoding(self, part):
        # set up a structure to hold regulators, type, coding
        #print part
        prmtr_name = part.__class__.__name__ 
        
        if isinstance(part,NegativePromoter):
            regulator = part.getRegulatedBy()[0]
            self.promoterMap.append( promoterMapping( prmtr_name, [regulator],[-1],[] ) )
        elif isinstance(part,PositivePromoter):
            regulator = part.getRegulatedBy()[0]
            self.promoterMap.append( promoterMapping( prmtr_name, [regulator],[1],[] ) )
        else:
            self.promoterMap.append( promoterMapping( prmtr_name, [None],[],[] ) )
    
        # loop over and get all coding regions until at end of parts list or we reach another promoter        
        nextpart = part.getAfterPart()
        while nextpart != None and isinstance(nextpart,Promoter) == False:
            
            if isinstance(nextpart,CodingRegion) == True:
                #print "\t", nextpart, "Coding"
                coding_name = nextpart.getCodingFor()[0]
                self.promoterMap[-1].coding.append( coding_name )
    
            nextpart = nextpart.getAfterPart()
            
    def notNeeded(self, part):
        #print part, "not needed for modelling"
        pass

    
    def buildPromoterMap(self, weaverOutput):
        if ReactionNetworks.builtMap == False:
            #print "Building promoter map"
            for part in weaverOutput.partList:
                part.generatePromoterMap()
            
            ReactionNetworks.builtMap = True 

        return self.promoterMap
        
# general calculation of a stochiometry matrix from a list of species and a list of reactions
def stoichiometry(nspecies, nreactions, species, reactions ):
    inputs = numpy.zeros(shape=(nreactions, nspecies))
    outputs = numpy.zeros(shape=(nreactions, nspecies))

    # map the species
    species_map = {}
    for i in range(nspecies):
        #print self.species[i], i
        species_map[ species[i] ] = i
        
    # loop over and calculate the reactants
    for j in range(nreactions):
        reactants = reactions[j].reactants
        for k in reactants:
            # The += allows for multimerisation
            inputs[j, species_map[k] ] += 1
        
    # loop over and calculate the products
    for j in range(nreactions):
        products = reactions[j].products
        for k in products:
            # handle empty string
            if len(k) > 0:
                # The += allows for multimerisation
                outputs[j, species_map[k.strip()] ] += 1
                    
    return outputs - inputs

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

class SheaAckersKineticsRNA(Aspect):
    
    def mainAspect(self):
        self.addWeaverOutput(self.getReactions)

        self.reactions = []
        self.nreactions = 0
        self.species = []
        self.nspecies = 0
        self.stoichiometry_matrix = 0
        self.locatedParts = {}

    def getReactions(self, weaverOutput):
        self.promoterMap = weaverOutput.buildPromoterMap()

        if getattr(weaverOutput, "getLocatedParts", None) != None:
            self.locatedParts = weaverOutput.getLocatedParts()
            print "located:", self.locatedParts 
        else:
            self.locatedParts = []

        self.getReactionsSheaAckersRNA()
        
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

        # assign mass action rates and parameters to certain processes
        for r in self.reactions:
            if r.process == "rnaDeg" or r.process == "proteinDeg" or r.process == "proteinTransl" or r.process == "complexDiss"\
               or r.process == "complexAss" or r.process == "complexDeg":
                r.assignMassAction()

        # calculate stoichiometry
        self.stoichiometry_matrix = stoichiometry(self.nspecies, self.nreactions, self.species, self.reactions)
       
        return [self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix]

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


class PrintReactionNetwork(Aspect):
    
    def mainAspect(self):
        self.addWeaverOutput(self.printReactionNetwork)

    def printReactionNetwork(self,weaverOutput):
        self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix = weaverOutput.getReactions()

        self.printReactions()
        self.printStoichiometryMatrix()

    def printReactions(self):
        print "\n\nGenerated set of reactions:"
        for i in range(self.nreactions):
            print self.reactions[i].reac_string.rjust(45," "), "\t\t\t", self.reactions[i].rate.rjust(35," "), "\t\t\t", self.reactions[i].param.rjust(10," "), "\t\t\t", self.reactions[i].process.rjust(15," ")
        
    def printStoichiometryMatrix(self):
        print "\n\nGenerated stoichiometry matrix:"
        print "".rjust(45," "),
        for j in range(self.nspecies):
            print repr(self.species[j]).rjust(15," "),
        print ""

        for i in range(self.nreactions):
            print self.reactions[i].reac_string.rjust(45," "),
            for j in range(self.nspecies):
                print repr(self.stoichiometry_matrix[i,j]).rjust(15," "),
            print ""
