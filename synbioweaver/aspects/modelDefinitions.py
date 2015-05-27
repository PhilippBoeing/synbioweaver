# definitions for creating a mathematical model representation of the system
import numpy
from copy import *

# Each reaction is assigned one of these processes
processes = [ "dnaBind", "dnaUnbind", "rnaTransc", "rnaDeg", "proteinTransl", "proteinExp", "proteinDeg", "complexAss", "complexDiss", "complexDeg", "context" ]

class Species:
    def __init__(self, scope, name):
        self.scope = scope
        self.name = name

    def __str__(self):
        return str(self.scope) + "::" + str(self.name)

class Model:
    def __init__(self, listOfSpecies, listOfReactions):
        self.species = listOfSpecies
        self.reactions = listOfReactions

        self.nspecies = len(self.species)
        self.nreactions = len(self.reactions)

        self.removeDuplicateSpecies()
        self.removeZeros()
        self.orderReactions()

        self.stoichiometry_matrix = self.stoichiometry(self.nspecies, self.nreactions, self.species, self.reactions)

    def orderReactions(self):
        pass

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
    param_counter = 0

    def __init__(self, reactants, products, process):
        
        # assign the products and reactants
        self.reactants = reactants
        self.products = products

        # assign a process
        if process in processes:
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
        # assign a parameter number
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
        


