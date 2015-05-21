from synbioweaver import *

declareNewMolecule('GFP')

class CodingGFP(Circuit):
    def mainCircuit(self):

        #self.addPart(Promoter)
        #self.addPart(RBS)
        #self.addPart(CodingRegion(GFP))
        #self.addPart(Terminator)

        #self.addPart(Promoter)
        ##self.addPart(RBS)
        #self.addPart(CodingRegion(GFP))
        #self.addPart(Terminator)

        self.addPart(Promoter)
        self.addPart(RBS)
        self.addPart(CodingRegion(GFP))
        self.addPart(Terminator)

        #self.addPart(Promoter)
        ##self.addPart(RBS)
        #self.addPart(CodingRegion(GFP))
        ##self.addPart(Terminator)
        
class DesignRulesRobust(Aspect):
    def mainAspect(self):
        anyCodingRegion = PartSignature('*.CodingRegion+')
        #terminatedCodingRegion =  PartSignature('*.CodingRegion+') % PartSignature('*.Terminator+')
        #unTerminatedCodingRegion =  PartSignature('*.CodingRegion+') % PointCutExpressionNot( PartSignature('*.Terminator+') )
        
        beforeCodingRegion = PointCut(anyCodingRegion,PointCut.BEFORE)
        afterCodingRegion = PointCut(anyCodingRegion,PointCut.AFTER)

        self.addAdvice(beforeCodingRegion,self.insertRBS)
        self.addAdvice(afterCodingRegion,self.insertTerminator)
        
    def printInfo(self,context):
        print "\tPRINT INFO::", context.part
        
    
    def insertRBS(self,context):
        # Check to see if any of the previous parts are RBS. If not add an RBS
        if all( isinstance(x,RBS) == False for x in context.part.before):
            self.addPart(RBS)
    
    def insertTerminator(self,context):
        print "\tPRINT INFO::", context
        print "\tPRINT INFO::", context.part
        print "\tPRINT INFO::", context.part.getAfterPart()
        print "\tPRINT INFO::", context.part.after
        
        #if all( isinstance(x,Terminator) == False for x in context.part.after):
        #    self.addPart(CodingRegion(GFP))
        #    self.addPart(Terminator)

compiledDesign = Weaver(CodingGFP,DesignRulesRobust).output()
print compiledDesign
