from aosb import *
import abcOutput

class SimpleSwitch(Circuit):
    def mainCircuit(self):

        declareNewMolecule('TetR')
        declareNewMolecule('LacI')
        declareNewMolecule('mCherry')
        declareNewMolecule('GFP')
        declareNewMolecule('TetR_TetR')
        declareNewMolecule('LacI_LacI')
        declareNewPart('PLtetO', NegativePromoter, [TetR_TetR])
        declareNewPart('Ptrc_2', NegativePromoter, [LacI_LacI])
        declareNewMolecule('PLtetO_TetR_TetR')
        declareNewMolecule('Ptrc_2_LacI_LacI')

        self.createMolecule(TetR)
        self.createMolecule(LacI)
        self.createMolecule(mCherry)
        self.createMolecule(GFP)
        self.createMolecule(Ptrc_2)
        self.createMolecule(PLtetO)
        self.createMolecule(TetR_TetR)
        self.createMolecule(LacI_LacI)
        self.createMolecule(PLtetO_TetR_TetR)
        self.createMolecule(Ptrc_2_LacI_LacI)

        self.addPart(Ptrc_2)
        self.addPart(CodingRegion(TetR))
        self.addPart(Ptrc_2)
        self.addPart(CodingRegion(mCherry))
        self.addPart(PLtetO)
        self.addPart(CodingRegion(LacI))
        self.addPart(PLtetO)
        self.addPart(CodingRegion(GFP))

        self.reactionFrom(PLtetO,TetR) >> self.reactionTo(PLtetO_TetR_TetR)
        self.reactionFrom(Ptrc_2,LacI) >> self.reactionTo(Ptrc_2_LacI_LacI)
        self.reactionFrom(TetR,TetR) >> self.reactionTo(TetR_TetR)
        self.reactionFrom(LacI,LacI) >> self.reactionTo(LacI_LacI)

compiledDesign = Weaver(SimpleSwitch, abcOutput.CudaOutput).output()
print compiledDesign.cudaOutput()
