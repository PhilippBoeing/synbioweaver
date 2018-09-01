from synbioweaver import *
import math

class MoleculeExpressionTrace():
    def __init__():
        var = 1

    def calculateTransfer(self,molecule):
        #print "calling calculateTransfer molecule", molecule
        
        result = 0
        # is this molecule expressed by a coding site?
        partCandidates = molecule.getBeforeNodes(CodingRegion)
        if (len(partCandidates) > 0):
            # ok, need to check if the coding region is expressed
            for part in partCandidates:
                codingRegion = part
                while True:
                    partCandidates = part.getBeforeNodes(Part)
                    if len(partCandidates) == 0:
                        break #the while...
                    part = partCandidates[0]
                    if isinstance(part,Terminator):
                        break #the while...
                    elif isinstance(part,PositivePromoter):
                        result = result + part.calculateTransfer()
                    elif isinstance(part,NegativePromoter):
                        # negative promoter transfer curve:
                        result = result + part.calculateTransfer()

            return result*codingRegion.codingRate

        else:
            # are there molecules connected to this molecule via import export?
            moleculeCandidates = molecule.getBeforeNodes(Molecule)
            if len(moleculeCandidates) > 0:
                result = 0
                for mol in moleculeCandidates:
                    result += mol.calculateTransfer()
                return result

            # otherwise...
            return molecule.signal

