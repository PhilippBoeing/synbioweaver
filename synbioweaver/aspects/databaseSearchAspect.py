from synbioweaver.core import *
from synbioweaver.parts import *

class PromoterSearch(Aspect):
    
    def mainAspect(self):
        self.addWeaverOutput(self.getLocatedParts)
        self.locatedParts = {}

    def getLocatedParts(self,weaverOutput):
        self.promoterMap = weaverOutput.buildPromoterMap()
        
        for prm in self.promoterMap:
            if str(prm.prmtr_name) in promoters:
                #print "located", str(prm.prmtr_name)
                
                class_name = "db" + str(prm.prmtr_name)
                module = __import__("parts")
                class_ = getattr(module, class_name)
                instance = class_()
                #print instance.transfer_function()
                
                # add the name of the part and an instance to the dictionary
                self.locatedParts[str(prm.prmtr_name)] = instance
        
        return self.locatedParts
    
