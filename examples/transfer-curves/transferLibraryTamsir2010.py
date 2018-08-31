from synbioweaver import *
import math

class MoleculeExpressionTrace(Aspect):
    def mainAspect(self):
        everyMoleculeSelector = MoleculeSignature('*.*')
        regulatedPromoterSelector = PartSignature('*.Promoter+(Molecule+)')
        codingRegionSelector = PartSignature('*.CodingRegion+(Molecule+)')
        self.addTypeAdvice(everyMoleculeSelector,0,'signal')
        self.addTypeAdvice(everyMoleculeSelector,self.calculateTransfer,'calculateTransfer')

        self.addTypeAdvice(PartSignature('*.PositivePromoter+(Ara)'),self.transferfuncpbad,'calculateTransfer')
        self.addTypeAdvice(PartSignature('*.PositivePromoter+(aTc)'),self.transferfuncptet,'calculateTransfer')
        self.addTypeAdvice(PartSignature('*.PositivePromoter+(LasI)'),self.transferfuncplas,'calculateTransfer')
        self.addTypeAdvice(PartSignature('*.PositivePromoter+(RhlI)'),self.transferfuncprhli,'calculateTransfer')

        self.addTypeAdvice(PartSignature('*.NegativePromoter+(Cl)'),self.transferfuncpcl,'calculateTransfer')

        self.addTypeAdvice(PartSignature('*.CodingRegion+(LasI)'),0.0005,'codingRate') # codingRate compared to YFP
        self.addTypeAdvice(PartSignature('*.CodingRegion+(RhlI)'),20,'codingRate') # codingRate compared to YFP
        self.addTypeAdvice(PartSignature('*.CodingRegion+(Cl)'),9,'codingRate') # codingRate compared to YFP
        self.addTypeAdvice(PartSignature('*.CodingRegion+(YFP)'),1,'codingRate') # codingRate compared to YFP

    def calculateTransfer(self,molecule):
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

    def transferfuncplas(self,part):
        l = part.getBeforeNodes(Molecule)[0].calculateTransfer()
        n = 1.4
        kd = 0.2 #uM
        k1 = 0.002
        k2c0 = 100
        cpart = math.pow(l,n)/(math.pow(kd,n)+math.pow(l,n))
        k2c = k2c0*cpart
        result = (k1 + k2c)/(1+ k1 + k2c)
        return result*690

    def transferfuncptet(self,part):
        l = part.getBeforeNodes(Molecule)[0].calculateTransfer()
        n = 1
        kd = 1.7 #ng/mL
        k1 = 350
        k2c0 = 160
        cpart = math.pow(l,n)/(math.pow(kd,n)+math.pow(l,n))
        cfpart = 1-cpart
        k2cf = k2c0*cfpart
        result = k1/(1+ k1 + 2 * k2cf + pow(k2cf,2))
        return result*3000

    def transferfuncprhli(self,part):
        l = part.getBeforeNodes(Molecule)[0].calculateTransfer()
        # no characterization, just buffer...
        return l

    def transferfuncpbad(self,part):
        l = part.getBeforeNodes(Molecule)[0].calculateTransfer()
        n = 2.8
        kd = 0.09 #mM
        k1 = 0.009
        k2c0 = 37.5
        k3c0 = 3.4
        cpart = math.pow(l,n)/(math.pow(kd,n)+math.pow(l,n))
        cfpart = 1-cpart
        result = (k1+k2c0*cpart)/(1+ k1 + 2 * k2c0*cpart + k3c0*cfpart)
        return result*7650

    def transferfuncpcl(self,part):
        r = part.getBeforeNodes(Molecule)[0].calculateTransfer()
        k1 = 350
        k2 = 0.015 #au
        k3 = 0.5 #au
        k0 = 0.18
        result = (k1)/(1+ k1 + k2*r + k3*r + (k2*k3*k0*r*r))
        return result*181
