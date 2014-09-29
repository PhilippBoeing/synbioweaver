from aosb import *

class TwoPromoters(Circuit):
    def mainCircuit(self):
        self.addPart(Promoter)
        self.addPart(NegativePromoter(Protein))
        
class PromoterTypeAdvice(Aspect):
    def mainAspect(self):
        # define three signatures
        nonRegulatedPromoterSignature = PartSignature('*.Promoter+()')
        regulatedPromoterSignature = PartSignature('*.Promoter+(Molecule+)')
        anyPromoterSignature = PartSignature('*.Promoter+')
        
        # add Boolean Attribute "isRegulated" to two different types of Promoters
        self.addTypeAdvice(nonRegulatedPromoterSignature,False,"isRegulated")
        self.addTypeAdvice(regulatedPromoterSignature,True,"isRegulated")
        
        # add a type advice method
        self.addTypeAdvice(anyPromoterSignature,self.printIsRegulated,"printIsRegulated")
        
    def printIsRegulated(self,part):
        if part.isRegulated == False:
            print str(part)+" is not regulated by a Molecule."
        else:
            print str(part)+" is regulated by a Molecule."

compiledDesign = Weaver(TwoPromoters,PromoterTypeAdvice).output()
for part in compiledDesign.partList:
    part.printIsRegulated()