from synbioweaver import *


class BooleanLogicTrace(Aspect):
    def mainAspect(self):
        everyMoleculeSelector = MoleculeSignature('*.*')
        regulatedPromoterSelector = PartSignature('*.Promoter+(Molecule+)')
        self.addTypeAdvice(everyMoleculeSelector,False,'moleculeExpressed')
        self.addTypeAdvice(everyMoleculeSelector,self.calculateExpressed,'calculateExpressed')
        self.addTypeAdvice(regulatedPromoterSelector,self.regulatedPromoterActivated,'isActivated')


    def calculateExpressed(self,molecule):
        moleculeCandidates = molecule.getBeforeNodes(Molecule)
        for mol in moleculeCandidates:
            if mol.calculateExpressed() == True:
                return True

        # is this molecule expressed by a coding site?
        partCandidates = molecule.getBeforeNodes(CodingRegion)
        if (len(partCandidates) > 0):
            # ok, need to check if the coding region is expressed
            for part in partCandidates:
                while True:
                    partCandidates = part.getBeforeNodes(Part)
                    if len(partCandidates) == 0:
                        break #the while...
                    part = partCandidates[0]
                    if isinstance(part,Terminator):
                        break #the while...
                    elif isinstance(part,ConstitutivePromoter):
                        return True
                    elif isinstance(part,PositivePromoter):
                        if part.isActivated():
                            return True
                            # however, if not, we continue until we hit a terminator
                            # the coding region might be activated by another promoter!
                    elif isinstance(part,NegativePromoter):
                        if not part.isActivated():
                            return True
            return False

        else:
            return molecule.moleculeExpressed

    def regulatedPromoterActivated(self,part):
        mol = part.getBeforeNodes(Molecule)[0]
        return mol.calculateExpressed()