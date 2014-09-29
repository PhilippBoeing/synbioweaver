from aosb import *

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
        
class NumberOfPartsAspect(object):
    def mainAspect(self):
        self.addWeaverOutput(self.printNumberOfParts)
        
    def printNumberOfParts(self,weaverOutput):
        print "The design has "+str(len(weaverOutput.partList))+" parts."

compiledDesign = Weaver(CodingGFP,DesignRules,NumberOfPartsAspect).output()
compiledDesign.printNumberOfParts()