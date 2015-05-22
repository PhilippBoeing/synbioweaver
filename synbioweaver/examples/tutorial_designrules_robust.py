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

        self.addPart(Promoter)
        ##self.addPart(RBS)
        #self.addPart(CodingRegion(GFP))
        ##self.addPart(Terminator)
        
class DesignRulesRobust(Aspect):
    def mainAspect(self):
        anyCodingRegion = PartSignature('*.CodingRegion+')
        anyTerminator = PartSignature('*.Terminator+')
        #unTerminatedCodingRegion =  PartSignature('*.CodingRegion+') % PointCutExpressionNot( PartSignature('*.Terminator+') )


        beforeCodingRegion = PointCut(anyCodingRegion,PointCut.BEFORE)
        afterCodingRegion = PointCut(anyCodingRegion,PointCut.AFTER)
        insteadOfTerminator = PointCut(PartSignature("!DesignRulesAspect.Terminator+"),PointCut.REPLACE)

        self.addAdvice(beforeCodingRegion,self.insertRBS)
        self.addAdvice(afterCodingRegion,self.insertTerminator)
        self.addAdvice(insteadOfTerminator,self.replaceTerminator)

    def printInfo(self,context):
        print "\tPRINT INFO::", context.part


    def insertRBS(self,context):
        # Check to see if any of the previous parts are RBS. If not add an RBS
        if all( isinstance(x,RBS) == False for x in context.part.before):
            self.addPart(RBS)

    def insertTerminator(self,context):
        self.addPart(Terminator)

    def replaceTerminator(self,context):
        if all( isinstance(x,Terminator) == False for x in context.part.before):
            return False #in this case, don't replace the terminator
        else:
            return True #in this case, replace the part with nothing

compiledDesign = Weaver(CodingGFP,DesignRulesRobust).output()
print compiledDesign
