from synbioweaver import *

declareNewMolecule('GFP')

class CodingGFP(Circuit):
    def __init__(self, regulator=None):
        super(CodingGFP,self).__init__()
        self.regulator = regulator

    def mainCircuit(self):
        if self.regulator is None:
            self.addPart(Promoter)
        else:
            self.createMolecule(self.regulator)
            self.addPart(PositivePromoter(self.regulator))
        self.addPart(CodingRegion(GFP))


class DesignRules(Aspect):
    def __init__(self, numberOfTerminators = 1):
        super(DesignRules,self).__init__()
        self.numberOfTerminators = numberOfTerminators

    def mainAspect(self):
        anyCodingRegion = PartSignature('*.CodingRegion+')

        beforeCodingRegion = PointCut(anyCodingRegion,PointCut.BEFORE)
        afterCodingRegion = PointCut(anyCodingRegion,PointCut.AFTER)

        self.addAdvice(beforeCodingRegion,self.insertRBS)
        self.addAdvice(afterCodingRegion,self.insertTerminators)

    def insertRBS(self,context):
        self.addPart(RBS)

    def insertTerminators(self,context):
        for _ in range(self.numberOfTerminators):
            self.addPart(Terminator)

compiledDesign = Weaver(CodingGFP,DesignRules).output()
print compiledDesign

declareNewMolecule("regulatorMolecule")

# now we're going to change the design using parameters to the circuits / designs
compiledDesign = Weaver(CodingGFP(regulatorMolecule),DesignRules(numberOfTerminators=3)).output()
print compiledDesign
