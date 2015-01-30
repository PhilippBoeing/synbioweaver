from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.runCudaSim import *
from synbioweaver.aspects.writeCudaFileODE import *
from synbioweaver.aspects.writeABCInputFileODE import *

# this example is based on the toggle switch implemented in Litcofsky et al. Nature Methods (2012)

class SimpleSwitch(Circuit):
    def mainCircuit(self):

        declareNewMolecule('TetR')
        declareNewMolecule('LacI')
        declareNewMolecule('mCherry')
        declareNewMolecule('GFP')
        declareNewMolecule('TetRn2')
        declareNewMolecule('LacIn2')
        declareNewMolecule('LacIn4')
        declareNewPart('PLtetO', NegativePromoter, [TetRn2])
        declareNewPart('Ptrc_2', NegativePromoter, [LacIn4])
        declareNewMolecule('PLtetO_TetRn2')
        declareNewMolecule('Ptrc_2_LacIn4')
        declareNewMolecule('aTc')
        declareNewMolecule('TetR_aTc')
        declareNewMolecule('IPTG')
        declareNewMolecule('LacIn4_IPTG')

        self.createMolecule(TetR)
        self.createMolecule(LacI)
        self.createMolecule(mCherry)
        self.createMolecule(GFP)
        self.createMolecule(TetRn2)
        self.createMolecule(LacIn2)
        self.createMolecule(LacIn4)
        self.createMolecule(Ptrc_2)
        self.createMolecule(PLtetO)
        self.createMolecule(PLtetO_TetRn2)
        self.createMolecule(Ptrc_2_LacIn4)
        self.createMolecule(aTc)
        self.createMolecule(TetR_aTc)
        self.createMolecule(IPTG)
        self.createMolecule(LacIn4_IPTG)


        self.addPart(Ptrc_2)
        self.addPart(CodingRegion(TetR))
        self.addPart(Ptrc_2)
        self.addPart(CodingRegion(mCherry))
        self.addPart(PLtetO)
        self.addPart(CodingRegion(LacI))
        self.addPart(PLtetO)
        self.addPart(CodingRegion(GFP))

        # these each add 2X -> X2, X2 -> 2X and X2 -> 0
        self.reactionFrom(TetR, TetR) >> self.reactionTo( TetRn2 )
        self.reactionFrom(LacI, LacI) >> self.reactionTo( LacIn2 )
        self.reactionFrom(LacIn2, LacIn2) >> self.reactionTo( LacIn4 )

        self.reactionFrom(aTc, TetR) >> self.reactionTo( TetR_aTc )
        self.reactionFrom(IPTG, LacIn4) >> self.reactionTo( LacIn4_IPTG )


compiledDesign5 = Weaver(SimpleSwitch, DesignRules, PromoterMapping, MassActionKineticsProtein, WriteABCInputFileODE, WriteCudaFileODE, RunCudaSim).output()

compiledDesign5.writeABCInputFileODE()
compiledDesign5.writeCudaFileODE()
#compiledDesign5.runCudaSim()