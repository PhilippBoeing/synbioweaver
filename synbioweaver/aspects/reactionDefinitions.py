# This contains the definitions of classes and methods needed across the reaction networks generation
import numpy

# Each reaction is assigned one of these processes
processes = [ "dnaBind", "dnaUnbind", "rnaTransc", "rnaDeg", "proteinTransl", "proteinExp", "proteinDeg", "complexAss", "complexDiss", "complexDeg", "context" ]

class Reaction:
    param_counter = 0

    def __init__(self, reactants, products, process):
        
        # assign the products and reactants
        self.reactants = map(str, reactants)
        self.products = map(str, products)

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
        

# This helper class associates promoters with downstream coding regions, regulators and polarities
#class promoterMapping:
#
#    def __init__(self, prmtr_name, regulators, polarities, coding):
#        self.prmtr_name = prmtr_name
#        self.regulators = regulators
#        self.polarities = polarities
#        self.coding = coding

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
            # forget anything not in the species map ie zero for example
            if k.strip() in species_map:
                # The += allows for multimerisation
                inputs[j, species_map[k] ] += 1
        
    # loop over and calculate the products
    for j in range(nreactions):
        products = reactions[j].products
        for k in products:
            # handle empty string
            if len(k) > 0:
                # forget anything not in the species map ie zero for example
                if k.strip() in species_map:
                    # The += allows for multimerisation
                    outputs[j, species_map[k.strip()] ] += 1
                    
    return outputs - inputs
