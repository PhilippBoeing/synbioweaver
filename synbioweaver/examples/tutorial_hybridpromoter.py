from synbioweaver import *

# follows tutorial_designrules.py but uses a HybridPromoter

class CodingGFP(Circuit):
    def mainCircuit(self):
        declareNewMolecule('GFP')
        declareNewMolecule('MoleculeA')
        declareNewMolecule('MoleculeB')
        declareNewMolecule('MoleculeC')
        self.createMolecule(MoleculeA)
        self.createMolecule(MoleculeB)
        self.createMolecule(MoleculeC)

        # self.addPart(HybridPromoter([MoleculeA,MoleculeB,MoleculeC],{MoleculeA:True,MoleculeB:False,MoleculeC:True}))
        declareNewPart("aVeryHybridPromoter",HybridPromoter,[MoleculeA,MoleculeB,MoleculeC],
                            regulatorInfoMap={MoleculeA:True,MoleculeB:False,MoleculeC:True})
        self.addPart(aVeryHybridPromoter)
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