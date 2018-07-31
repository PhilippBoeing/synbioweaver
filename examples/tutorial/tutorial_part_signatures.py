# This example runs the examples from figure 5 of the manuscript

from synbioweaver import *

class BBa_B0030(RBS):
    def __init__(self):
        super(BBa_B0030, self).__init__()
        self.type = "biobrick"
        self.sequence = "attaaagaggagaaa"

declareNewMolecule('GFP')
declareNewMolecule('TetR')

declareNewPart('pConst', ConstitutivePromoter, [])
declareNewPart('cTetR', CodingRegion, moleculesAfter=[TetR])
declareNewPart('ter',Terminator)
declareNewPart('pNeg', NegativePromoter, [TetR])
declareNewPart('rbs',RBS)

class SimpleCircuit(Circuit):
    def mainCircuit(self):
        self.createMolecule(TetR)

        self.addPart(pConst)
        self.addPart(BBa_B0030)
        self.addPart(cTetR)
        self.addPart(Terminator)
        self.addPart(pNeg)
        self.addPart(RBS)
        self.addPart(CodingRegion(GFP))
        self.addPart(Terminator)

class ExaminePartSignatures(Aspect):
    def __init__(self, example = "A"):
        super(ExaminePartSignatures,self).__init__()
        self.example = example

    def mainAspect(self):
        #anyCodingRegion = PartSignature('*.CodingRegion+')
        
        #beforeCodingRegion =
        #afterCodingRegion = PointCut(anyCodingRegion,PointCut.AFTER)
        
        #self.id = "A"
        print "Printing context for advice:", self.example
        if self.example == "A":
            self.addAdvice(PointCut('SimpleCircuit.BBa_B0030',PointCut.BEFORE) ,self.printContext, 100)
        if self.example == "B":
            self.addAdvice(PointCut('Simple*.Promoter+',PointCut.BEFORE) ,self.printContext)
        if self.example == "C":
            self.addAdvice(PointCut('Simple*.Promoter+(TetR)',PointCut.BEFORE) ,self.printContext)
        if self.example == "D":
            self.addAdvice(PointCut('Simple*.*(TetR)',PointCut.BEFORE) ,self.printContext,1)
        if self.example == "E":
            self.addAdvice(PointCut('Simple*.*(Protein+)',PointCut.BEFORE) ,self.printContext,1)
        if self.example == "F":
            self.addAdvice(PointCut('Simple*.Promoter+()',PointCut.BEFORE) ,self.printContext,1)
        if self.example == "G":
            self.addAdvice(PointCut('Simple*.!Terminator()',PointCut.BEFORE) ,self.printContext,1)
        if self.example == "H":
            self.addAdvice(PointCut(PartSignature('*.BBa_B0030') % PartSignature('*.CodingRegion+') ,PointCut.AFTER) ,self.printContext)

    def printContext(self,context):
        
        print context.part
    

compiledDesign = Weaver(SimpleCircuit,ExaminePartSignatures("A")).output()
compiledDesign = Weaver(SimpleCircuit,ExaminePartSignatures("B")).output()
compiledDesign = Weaver(SimpleCircuit,ExaminePartSignatures("C")).output()
compiledDesign = Weaver(SimpleCircuit,ExaminePartSignatures("D")).output()
compiledDesign = Weaver(SimpleCircuit,ExaminePartSignatures("E")).output()
compiledDesign = Weaver(SimpleCircuit,ExaminePartSignatures("F")).output()
compiledDesign = Weaver(SimpleCircuit,ExaminePartSignatures("G")).output()
compiledDesign = Weaver(SimpleCircuit,ExaminePartSignatures("H")).output()
#print compiledDesign
