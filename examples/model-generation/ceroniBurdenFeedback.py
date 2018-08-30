from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from synbioweaver.aspects.promoterMappingAspect import *
from synbioweaver.aspects.MassActionKineticsRNAAspect import *
from synbioweaver.aspects.printReactionNetworkAspect import *

#class RNACodingRegion(Part):
#    def __init__(self, codesFor):
#        super(RNACodingRegion, self).__init__()
#        self.precompileMoleculesAfter.append(checkAndSetMolecule(codesFor))
#
#    def getCodingFor(self):
#        result = self.getAfterNodes(Molecule)
#        if len(result) > 0:
#            return result

#        return self.precompileMoleculesAfter
        
declareNewMolecule('sgRNA')
declareNewMolecule('dCas9')
declareNewMolecule('sgRNAvdCas9')
declareNewMolecule('Bdn')
declareNewMolecule('mCherry')

declareNewPart('Pc', ConstitutivePromoter)
declareNewPart('PhtpG1', PositivePromoter, [Bdn])
declareNewPart('P_BAD', NegativePromoter, [sgRNAvdCas9])

class circuit(Circuit):
    def mainCircuit(self):
        self.createMolecule(Bdn)
        self.createMolecule(sgRNAvdCas9)
        
        self.addPart(Pc)
        self.addPart(RBS)
        self.addPart(CodingRegion(dCas9))
        self.addPart(Terminator)

        self.addPart(PhtpG1)
        self.addPart(RNACodingRegion(sgRNA))
        self.addPart(Terminator)

        self.addPart(P_BAD)
        self.addPart(RBS)
        self.addPart(RNACodingRegion(mCherry))
        self.addPart(Terminator)
        
        self.reactionFrom(sgRNA,dCas9) >> self.reactionTo(sgRNAvdCas9)
        self.reactionFrom(mCherry) >> self.reactionTo(Bdn, mCherry)
        self.reactionFrom(dCas9) >> self.reactionTo(Bdn, dCas9)
        

compiledDesign = Weaver(circuit, PromoterMapping, MassActionKineticsRNA, PrintReactionNetwork).output()
compiledDesign.printReactionNetwork()



