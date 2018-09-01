from synbioweaver import *
from moleculeExpressionTrace import *
import math

class LibraryRosenfeld2005(Aspect, MoleculeExpressionTrace):
    def mainAspect(self):
        everyMoleculeSelector = MoleculeSignature('*.*')
        regulatedPromoterSelector = PartSignature('*.Promoter+(Molecule+)')
        codingRegionSelector = PartSignature('*.CodingRegion+(Molecule+)')
        self.addTypeAdvice(everyMoleculeSelector,0,'signal')
        self.addTypeAdvice(everyMoleculeSelector,self.calculateTransfer,'calculateTransfer')

        self.addTypeAdvice(PartSignature('*.PR+(CI)'),self.transferPR,'calculateTransfer')
        self.addTypeAdvice(PartSignature('*.PR_OR2+(CI)'),self.transferPR_OR2,'calculateTransfer')
        
        self.addTypeAdvice(PartSignature('*.CodingRegion+(CI)'),1,'codingRate') # codingRate compared to YFP
        self.addTypeAdvice(PartSignature('*.CodingRegion+(CFP)'),1,'codingRate') # codingRate compared to YFP

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

   
