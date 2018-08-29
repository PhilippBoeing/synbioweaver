from synbioweaver import *

declareNewMolecule('GFP')
declareNewMolecule('B1')
declareNewMolecule('B2')
declareNewMolecule('A1')
declareNewMolecule('A2')
declareNewMolecule('R1')

class CodingGFP(Circuit):
    def mainCircuit(self):
        self.addPart(Promoter)
        self.addPart(RBS)
        self.addPart(CodingRegion(GFP))
        self.addPart(Terminator)

class PrecedenceRules(Aspect):
    def mainAspect(self):
        CodingRegion = PartSignature('CodingGFP.CodingRegion+(GFP)')
        
        beforeCodingRegion1 = PointCut(CodingRegion,PointCut.BEFORE)
        beforeCodingRegion2 = PointCut(CodingRegion,PointCut.BEFORE)

        afterCodingRegion1 = PointCut(CodingRegion,PointCut.AFTER)
        afterCodingRegion2 = PointCut(CodingRegion,PointCut.AFTER)

        atCodingRegion1 = PointCut(CodingRegion,PointCut.REPLACE)

        # The resultant circuit depends on the relative precedence of B2, A2 and R1
        # if B2,A2 > R1 then we end up with B1-B2-R1-A2-A1
        # if B2,A2 < R1 then we end up with B1-R1-A1
        
        self.addAdvice(beforeCodingRegion1,self.insertCodingB1,100)
        self.addAdvice(beforeCodingRegion2,self.insertCodingB2,50)

        self.addAdvice(afterCodingRegion1,self.insertCodingA1,100)
        self.addAdvice(afterCodingRegion2,self.insertCodingA2,50)

        self.addAdvice(atCodingRegion1,self.insertCodingR1,60)
    
    def insertCodingB1(self,context):
        self.addPart(CodingRegion(B1))
    
    def insertCodingB2(self,context):
        self.addPart(CodingRegion(B2))

    def insertCodingA1(self,context):
        self.addPart(CodingRegion(A1))

    def insertCodingA2(self,context):
        self.addPart(CodingRegion(A2))

    def insertCodingR1(self,context):
        self.addPart(CodingRegion(R1))
        return True

compiledDesign = Weaver(CodingGFP,PrecedenceRules).output()
print compiledDesign
