from synbioweaver import *

declareNewMolecule('GFP')

class CodingGFP(Circuit):
    def mainCircuit(self):
        self.addPart(Promoter)
        #self.addPart(RBS)
        self.addPart(CodingRegion(GFP))
        #self.addPart(Terminator)
        
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


class CheckParts(Aspect):
    def mainAspect(self):
        # match any terminator that has been added by DesignRules
        addedTerminator = PartSignature("DesignRules.Terminator+")

        insteadOfTerminator = PointCut(addedTerminator,PointCut.REPLACE)

        self.addAdvice(insteadOfTerminator,self.replaceTerminator)

    def replaceTerminator(self,context):
        print "HERE:"
        if all( isinstance(x,Terminator) == False for x in context.part.before):
            return False #in this case, don't replace the terminator
        else:
            return True #in this case, replace the part with nothing

compiledDesign = Weaver(CodingGFP,DesignRules,CheckParts).output()
print compiledDesign
