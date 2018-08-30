from synbioweaver.core import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsProteinAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

declareNewMolecule('GFP')
declareNewMolecule('HK')
declareNewMolecule('HKS')
declareNewMolecule('HKp')
declareNewMolecule('RR')
declareNewMolecule('HKpRR')
declareNewMolecule('RRp')
declareNewMolecule('S')

declareNewPart('Pc', ConstitutivePromoter)
declareNewPart('Psignal', PositivePromoter, [RRp] )

class biosensorTCS(Circuit):
    def mainCircuit(self):
        self.createMolecule(S)
        self.createMolecule(HKp)
        self.createMolecule(HKS)
        self.createMolecule(RRp)

        # Operon for the two component signalling system
        self.addPart(Pc)
        self.addPart(RBS)
        self.addPart( CodingRegion(HK) )
        self.addPart(RBS)
        self.addPart( CodingRegion(RR) )
        self.addPart(Terminator)

        # Readout of the signal
        self.addPart(Psignal)
        self.addPart(RBS)
        self.addPart( CodingRegion(GFP) )
        self.addPart(Terminator)
        
        # Signal, S, binding to histidine kinase
        self.reactionFrom(S,HK) >> self.reactionTo(HKS)
        self.reactionFrom(HKS) >> self.reactionTo(HKp,S)

        # Phosphorylation of the response regulator
        self.reactionFrom(HKp,RR) >> self.reactionTo(HKpRR)
        self.reactionFrom(HKpRR) >> self.reactionTo(HKp,RRp)

        # Dephosphorylation reactions
        self.reactionFrom(HKp) >> self.reactionTo(HK)
        self.reactionFrom(RRp) >> self.reactionTo(RR)

        
compiledDesign = Weaver(biosensorTCS, PromoterMapping, MassActionKineticsProtein, PrintReactionNetwork).output()
compiledDesign.printReactionNetwork()

