from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.MassActionKineticsRNAAspect import *
from synbioweaver.aspects.SheaAckersKineticsRNAAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *
from synbioweaver.aspects.runCudaSim import *
from synbioweaver.aspects.writeCudaFileODE import *


declareNewMolecule('GFP')
declareNewMolecule('LacI')
declareNewMolecule('AraC')
declareNewMolecule('LacIn4')
declareNewMolecule('AraCn2')
declareNewMolecule('zero')

declareNewPart("P1", HybridPromoter, [LacIn4, AraCn2], regulatorInfoMap={LacIn4:False, AraCn2:True})
declareNewPart("P2", HybridPromoter, [LacIn4, AraCn2], regulatorInfoMap={LacIn4:False, AraCn2:True})
declareNewPart("P3", HybridPromoter, [LacIn4, AraCn2], regulatorInfoMap={LacIn4:False, AraCn2:True})

declareNewMolecule('A')
declareNewMolecule('B')
declareNewMolecule('An2')
declareNewMolecule('Bn2')
declareNewMolecule('In')

declareNewPart('Pa', NegativePromoter, [An2])
declareNewPart('Pb', NegativePromoter, [Bn2])

class S(Circuit):
    def mainCircuit(self):

        self.createMolecule(GFP)
        self.createMolecule(LacI)
        self.createMolecule(AraC)
        self.createMolecule(LacIn4)
        self.createMolecule(AraCn2)

        self.addPart(P1)
        self.addPart(CodingRegion(AraC))

        self.addPart(P2)
        self.addPart(CodingRegion(LacI))

        self.addPart(P3)
        self.addPart(CodingRegion(GFP))

        self.reactionFrom(LacI, LacI, LacI, LacI) >> self.reactionTo( LacIn4 )
        self.reactionFrom(LacIn4) >> self.reactionTo( zero )
        self.reactionFrom(AraC, AraC) >> self.reactionTo( AraCn2 )
        self.reactionFrom(AraCn2) >> self.reactionTo( zero )

class T(Circuit):
    def mainCircuit(self):
        
        self.createMolecule(A)
        self.createMolecule(An2)
        self.createMolecule(B)
        self.createMolecule(Bn2)
        self.createMolecule(In)
        
        self.addPart(Pa)
        self.addPart(CodingRegion(B))
        self.addPart(Pb)
        self.addPart(CodingRegion(A))
        
        self.reactionFrom(A, A) >> self.reactionTo( An2 )
        self.reactionFrom(B, B) >> self.reactionTo( Bn2 )
        self.reactionFrom(An2) >> self.reactionTo( zero )
        self.reactionFrom(Bn2) >> self.reactionTo( zero )

        self.reactionFrom(zero) >> self.reactionTo( In )
        self.reactionFrom(In,Bn2) >> self.reactionTo( In )

class toggleOsc(Circuit):
    def mainCircuit(self):
        self.addCircuit(S)
        self.addCircuit(T)

        
compiledDesign = Weaver(toggleOsc, DesignRules, PromoterMapping, SheaAckersKineticsRNA, PrintReactionNetwork, WriteCudaFileODE
                        #RunCudaSim(name="growth",vpars=sim_prior, inits=sim_inits, times=times, obs=sim_obs) 
                        ).output()

compiledDesign.printReactionNetwork()
cudaStr = compiledDesign.writeCudaFileODE()
cudaFile = open("toggleOsc_combined.cu","w")
print >>cudaFile, cudaStr
cudaFile.close()
