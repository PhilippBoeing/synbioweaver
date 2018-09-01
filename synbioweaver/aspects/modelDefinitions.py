# definitions for creating a mathematical model representation of the system
import numpy
from copy import *

# Each reaction is assigned one of these processes
processes = { "rnaTransc":1, "rnaDeg":2, "proteinTransl":3, "proteinExp":4, "proteinDeg":5, "dnaBind":6, "dnaUnbind":7, "complexAss":8, "complexDiss":9, "complexDeg":10, "context":11 }

species_type = { "mRNA":1, "protein":2, "promoter":3, "bindingDNA":4, "molecule":5, "general":6 }

class Species:
    def __init__(self, scope, name, type="general"):
        self.scope = scope
        self.name = name
        self.type = type

    #def __str__(self):
    #    return str(self.scope) + "::" + str(self.name)

    def __str__(self):
        return str(self.name)

def sort_react( x ):
    # get either reactant or product
    if len(x.reactants) != 0:
        tosort = x.reactants[0]
    else:
        tosort = x.products[0]
        
    return ( tosort.scope, processes[ x.process ], tosort.name )

class Model:
    def __init__(self, listOfSpecies, listOfReactions, listOfParameters):
        self.species = listOfSpecies
        self.reactions = listOfReactions
        self.parameters = listOfParameters

        self.nspecies = len(self.species)
        self.nreactions = len(self.reactions)

        self.removeDuplicateSpecies()
        self.removeZeros()
        self.orderReactions()
        self.orderSpecies()

        # Assign parameters based on current ordering

        self.stoichiometry_matrix = self.stoichiometry(self.nspecies, self.nreactions, self.species, self.reactions)

    def orderSpecies(self):
        # order alphabetically, could order by species type in the future
        self.species = sorted( self.species, key=lambda spc: ( str(spc.scope), species_type[ spc.type ], str(spc.name)) )
        #print [x.name for x in self.species]

    def orderReactions(self):
        # order by namespace, process, alphabetical
        #self.reactions = sorted( self.reactions, key=lambda reac: (processes[ reac.process ], str(reac.reactants[0])) )
        self.reactions = sorted( self.reactions, key=sort_react )

    def removeDuplicateSpecies(self):
        #print "before dup", map(str,self.species)
        seen = set() # holds the observed scope::names
        uniq = []    # holds the unique species objects
        for x in self.species:
            if str(x) not in seen:
                uniq.append(x)
                seen.add(str(x))

        self.species = copy(uniq)
        self.nspecies = len(self.species)
        #print "after dup", map(str,self.species) 

    def removeZeros(self):

        ########## species
        nonZero = []
        for x in self.species:
            if not x.name == "zero":
                nonZero.append(x)
                
        self.species = copy(nonZero)
        self.nspecies = len(self.species)
        ## print "after zero", map(str,self.species) 

        ########## reactions
        for i in range(self.nreactions):
            reacs = self.reactions[i].reactants
            newreacs = [x for x in reacs if x.name != "zero"]

            prods = self.reactions[i].products
            newprods = [x for x in prods if x.name != "zero"]

            self.reactions[i].reactants = newreacs
            self.reactions[i].products = newprods   

        return

    # general calculation of a stochiometry matrix from a list of species and a list of reactions
    def stoichiometry(self, nspecies, nreactions, species, reactions ):
        inputs = numpy.zeros(shape=(nreactions, nspecies))
        outputs = numpy.zeros(shape=(nreactions, nspecies))

        # map the species from name to ordering in the list
        # calling str( Species ) combines the name with the scope
        species_map = {}
        for i in range(nspecies):
            ## overloaded str operator will combine scope::name
            species_map[ str(species[i]) ] = i
        #print "Model::stoichiometry", species_map
        
        # loop over and calculate the reactants
        for j in range(nreactions):
            reactants = reactions[j].reactants
            for k in reactants:
                # forget anything not in the species map ie zero for example
                if str(k).strip() in species_map:
                    # The += allows for multimerisation
                    inputs[j, species_map[ str(k).strip() ] ] += 1
        
        # loop over and calculate the products
        for j in range(nreactions):
            products = reactions[j].products
            for k in products:
                # handle empty string
                if len(str(k)) > 0:
                    # forget anything not in the species map ie zero for example
                    if str(k).strip() in species_map:
                        # The += allows for multimerisation
                        outputs[j, species_map[ str(k).strip() ] ] += 1
                    
        return outputs - inputs



class Reaction:
    # a global parameter counter to number parameters
    param_counter = 0

    def __init__(self, reactants, products, process):
        
        # assign the products and reactants
        self.reactants = reactants
        self.products = products

        # assign a process
        if process in processes.keys():
            self.process = process
        else:
            print "Error: process undefined:", process
            exit(1)

        self.param = []
        self.rate = ""

    # returns a string representation of the reaction
    def reactionString(self):
        # create a reaction string for output
        react = ' + '.join(map(str, self.reactants))
        prod = ' + '.join(map(str, self.products))
        return react + ' -> ' + prod

    # returns a string representation of the parameters
    def paramString(self):
        return ', '.join(self.param)

    def assignMassAction(self):
        Reaction.param_counter += 1
        
        # here there is only one parameter per reactions
        par = 'p'+str(Reaction.param_counter)

        self.param.append( par )

        # assign a rate
        if len(self.reactants) > 0:
            rt = '*'.join(map(str, self.reactants))
            self.rate = rt+'*' + par
        else:
            self.rate = par 

    def assignMassActionExisting(self, par):
        self.param.append( par )

        # assign a rate
        if len(self.reactants) > 0:
            rt = '*'.join(map(str, self.reactants))
            self.rate = rt+'*' + par
        else:
            self.rate = par

    def assignSA(self, prmtrMapping):
        totalregs = len( prmtrMapping.getRegulators() )
        totalpos = 0
        Reaction.param_counter += 1
        par = "p"+str(Reaction.param_counter)

        self.param.append( par )

        num = "( " + par;
        denom = "/( 1 + " + par 
        for i in range(0,totalregs):
            Reaction.param_counter += 1
            par = "p"+str(Reaction.param_counter)
            self.param.append( par )

            denom = denom + " + " + par + "*" + str( prmtrMapping.getRegulators()[i] )
                
            if prmtrMapping.getPolarities()[i] == +1:
                totalpos += 1
                num = num + " + " + par + "*" + str( prmtrMapping.getRegulators()[i] )

        num = num + " )"
        denom = denom + " )"
    
        Reaction.param_counter += 1
        par = "p"+str(Reaction.param_counter)
        self.param.append( par )

        self.rate = par + "*( " +  num+denom + ")"
        


