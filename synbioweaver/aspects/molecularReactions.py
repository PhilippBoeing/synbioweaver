from synbioweaver.core import *
#from synbioweaver.aspects.reactionDefinitions import *
from synbioweaver.aspects.modelDefinitions import *

from synbioweaver.aspects.promoterMappingAspect import *

# This class handles the molecular interactions for the reaction abstractions
# Should probably be a proper abstract base class

class MolecularReactions:
    
    def __init__(self):
        pass

    # Helper function
    def examineContext(self, mol):
        print str(mol.scope.circuitName), str(mol)
        
        for j in range(len(mol.before)):
            if type(mol.before[j]) == type(list()):
                print "\tdetected before list:",
                for k in range( len(mol.before[j]) ):
                    print mol.before[j][k].scope.circuitName, mol.before[j][k],
                print ""
            else:
                print "\tdetected before sing:", mol.before[j].scope.circuitName, mol.before[j]

        for j in range(len(mol.after)):
            if type(mol.after[j]) == type(list()):
                print "\tdetected after list:",
                for k in range( len(mol.after[j]) ):
                    print mol.after[j][k].scope.circuitName, mol.after[j][k],
                print ""
            else:
                print "\tdetected after sing:", mol.after[j]

    # Return a species object from the molecule
    def getMoleculeId(self, mol):
        #if PromoterMapping.multiComp == True and str(mol) != "zero":
        #    return str(mol.scope.circuitName) + "-" + str(mol)
        #else:
        #    return str(mol)
        return( Species(mol.scope.circuitName, str(mol), "molecule") )

    def addMolecularReaction(self, mol):
        # BEFORE -> mol
        for j in range(len(mol.before)):
            if type(mol.before[j]) == type(list()):
                reacs = []
                for k in range( len(mol.before[j]) ):
                    #print mol.before[j][k].scope.circuitName, mol.before[j][k],
                    reacs.append( self.getMoleculeId(mol.before[j][k]) )

                self.reactions.append( Reaction(reacs, [self.getMoleculeId(mol)], "complexAss") )
            else:
                # check here that it isn't a part -> molecule reaction which we will handle separately
                if isinstance(mol.before[j], Molecule) == True:
                    #print "\tdetected before sing:", mol.before[j].scope.circuitName, mol.before[j]
                    self.reactions.append( Reaction([ self.getMoleculeId(mol.before[j]) ], [ self.getMoleculeId(mol) ], "complexAss") )

        
        # mol -> AFTER
        for j in range(len(mol.after)):
            if type(mol.after[j]) == type(list()):
                prods = []
                for k in range( len(mol.after[j]) ):
                    #print mol.after[j][k].scope.circuitName, mol.after[j][k],
                    prods.append( self.getMoleculeId(mol.after[j][k]) )
                
                self.reactions.append( Reaction([ self.getMoleculeId(mol) ], prods, "complexDiss") )
            else:
                # check here that it isn't a part -> molecule reaction which we will handle separately
                if isinstance(mol.after[j], Molecule) == True:
                    #print "\tdetected before sing:", mol.after[j].scope.circuitName, mol.after[j]
                    self.reactions.append( Reaction([ self.getMoleculeId(mol) ], [ self.getMoleculeId(mol.after[j]) ], "complexDiss") )

    def removeDuplicateReactions(self):
        
        dups = []
        for i in range(len(self.reactions)):
            for j in range(i+1,len(self.reactions)):
                #print self.reactions[i].reactionString(), self.reactions[j].reactionString()
                if self.reactions[i].reactionString() == self.reactions[j].reactionString():
                    #print j, self.reactions[j].reactionString(), " is duplicated"
                    dups.append( j )
         
        self.reactions = [ i for j,i in enumerate( self.reactions ) if j not in dups ]
           
    def getMolecules(self, weaverOutput):
            
        for i in range(len(weaverOutput.moleculeList)):
            # Add all molecules to list then do a unique
            sp = self.getMoleculeId(weaverOutput.moleculeList[i])
            # cleave off the regulated by
            # self.species.append( sp.split("(")[0] )
            self.species.append(sp)

            #self.examineContext( weaverOutput.moleculeList[i] )
            self.addMolecularReaction(weaverOutput.moleculeList[i])

        if len(weaverOutput.wovenCircuitList) > 0:
            for circ in weaverOutput.wovenCircuitList:
                for i in range(len(circ.moleculeList)):
                    self.species.append(self.getMoleculeId(circ.moleculeList[i]) )
                    #self.examineContext( circ.moleculeList[i] )
                    self.addMolecularReaction(circ.moleculeList[i])

        # transfer reactions A -> B will appear twice so may need to remove duplicate reactions
        self.removeDuplicateReactions()
         
                    
