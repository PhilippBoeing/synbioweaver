from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import DesignRules

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

        self.addTypeAdvice(PartSignature('*.Terminator'),self.terminatorRuleString,'terminatorRuleString')
        
        self.addAdvice(beforeCodingRegion,self.insertRBS)
        self.addAdvice(afterCodingRegion,self.insertTerminator)

    def terminatorRuleString(self,part):
        print 'hello\n'
        for before in part.before:
            print before

    def insertRBS(self,context):
        self.addPart(RBS)
    
    def insertTerminator(self,context):
        self.addPart(Terminator)
        
class NumberOfPartsAspect(Aspect):
    def mainAspect(self):
        self.addWeaverOutput(self.printNumberOfParts)
        
    def printNumberOfParts(self,weaverOutput):
        print "The design has "+str(len(weaverOutput.partList))+" parts."

compiledDesign = Weaver(CodingGFP,DesignRules,NumberOfPartsAspect).output()
compiledDesign.printNumberOfParts()
