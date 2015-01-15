from synbioweaver import *
from synbioweaver.designaspects import DesignRules

class CodingGFP(Circuit):
    def mainCircuit(self):
        declareNewMolecule('GFP') 
        self.addPart(Promoter)
        self.addPart(CodingRegion(GFP))
        
class DesignRules(Aspect):
    def mainAspect(self):
        anyCodingRegion = PartSignature('*.CodingRegion+')
        
        beforeCodingRegion = PointCut(anyCodingRegion,PointCut.BEFORE)
        afterCodingRegion = PointCut(anyCodingRegion,PointCut.AFTER)
        
        self.addAdvice(beforeCodingRegion,self.insertRBS)
        self.addAdvice(afterCodingRegion,self.insertTerminator)
    
    def insertRBS(self,context):
        self.addPart(RBS)
    
    def insertTerminator(self,context):
        self.addPart(Terminator)
        
class PrintStack(Aspect):
    def mainAspect(self):
        self.addAdvice(PointCut('*.*',PointCut.AFTER),self.printStack)
    
    def printStack(self,context):
        print ','.join(str(context.within[i]) for i in range(len(context.within))) + " add: " + str(context.part)
    

compiledDesign = Weaver(CodingGFP,DesignRules,PrintStack).output()
print compiledDesign