from synbioweaver import *

declareNewMolecule('GFP')

class CodingGFP(Circuit):
    def mainCircuit(self):
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

compiledDesign = Weaver(CodingGFP,DesignRules).output()
print compiledDesign