from synbioweaver import *
import math

class MoleculeExpressionTrace(Aspect):
    def mainAspect(self):
        everyMoleculeSelector = MoleculeSignature('*.*')
        regulatedPromoterSelector = PartSignature('*.Promoter+(Molecule+)')
        codingRegionSelector = PartSignature('*.CodingRegion+(Molecule+)')
        self.addTypeAdvice(everyMoleculeSelector,0,'signal')
        self.addTypeAdvice(everyMoleculeSelector,self.calculateTransfer,'calculateTransfer')

        #self.addTypeAdvice(PartSignature('*.PositivePromoter+(aTc)'),self.transferfuncptet,'calculateTransfer')

        self.addTypeAdvice(PartSignature('*.PR+(CI)'),self.transferPR,'calculateTransfer')
        self.addTypeAdvice(PartSignature('*.PR_OR2+(CI)'),self.transferPR_OR2,'calculateTransfer')
        
        self.addTypeAdvice(PartSignature('*.CodingRegion+(CI)'),1,'codingRate') # codingRate compared to YFP
        self.addTypeAdvice(PartSignature('*.CodingRegion+(CFP)'),1,'codingRate') # codingRate compared to YFP

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

    def transferPR(self,part):
        l = part.getBeforeNodes(Molecule)[0].calculateTransfer()
        n = 2.4
        kd = 55
        beta = 220
        result = beta/(1 + math.pow(l/kd,n))
        return result

    def transferPR_OR2(self,part):
        l = part.getBeforeNodes(Molecule)[0].calculateTransfer()
        n = 1.7
        kd = 120
        beta = 255
        result = beta/(1 + math.pow(l/kd,n))
        return result

   
