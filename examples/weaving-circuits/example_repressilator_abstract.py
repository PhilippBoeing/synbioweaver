from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *

# GFP circuit definitions
declareNewMolecule('AbstractMolecule')
declareNewMolecule('externalInducer')
declareNewMolecule('GFP')

# Abstract repressilator definitions
declareNewMolecule('A',AbstractMolecule)
declareNewMolecule('B',AbstractMolecule)
declareNewMolecule('C',AbstractMolecule)

# AND gate design based on Wang et al, 2011, 
# Engineering modular and orthogonal genetic logic gates for robust digital-like synthetic biology
# Nature communications
declareNewMolecule('HrpR')
declareNewMolecule('HrpS')
declareNewMolecule('HrpSandR') # not ideal...
declareNewPart('hrpL', PositivePromoter, [HrpSandR] ) # not ideal...

class InduceGFP(Circuit):
    # code GFP dependent on an external inducer
    def mainCircuit(self):
        self.createMolecule(externalInducer)
        self.addPart(PositivePromoter(externalInducer))
        self.addPart(CodingRegion(GFP))
        
class CodingGFP(Circuit):
    # constitutively code gfp
    def mainCircuit(self):
        self.addPart(ConstitutivePromoter)
        self.addPart(CodingRegion(GFP))
        
class AbstractRepressilatorComposite(Circuit):
    # an abstract repressilator circuit
    #       A(C), B(A), C(B)
           
    def declareRepressilatorProteins(self):
        self.repressilationProteins = []
        self.repressilationProteins.append(A)
        self.repressilationProteins.append(B)
        self.repressilationProteins.append(C)
        
    def getSignal(self):
        return A
        
    def mainCircuit(self):
        self.createMolecule(A)
        self.createMolecule(B)
        self.createMolecule(C)
        self.declareRepressilatorProteins()
        lastprotein = self.repressilationProteins[-1]
        
        for protein in self.repressilationProteins:
            self.addPart(NegativePromoter(lastprotein))
            self.addPart(CodingRegion(protein))
            lastprotein = protein

class InstantiateAbstractMolecules(Aspect):
    # from a list of non-cross-talking transcription factors
    # instantiate actual molecules in an abstract gene regulatory network
    def mainAspect(self):
        self.independentMolecules = [declareNewMolecule('TetR'),declareNewMolecule('LacL'),declareNewMolecule('cl')]
        
        self.abstractConcreteMapping = {}
        self.addAdvice(PointCut('*.*(AbstractMolecule+)',PointCut.BEFORE),self.makeConcrete,30)

    def makeConcrete(self,context):
        #print "makeConcrete::", context.part, context.part.precompileMoleculesBefore, context.part.precompileMoleculesAfter

        for i, m in enumerate(context.part.precompileMoleculesBefore):
            if not m in self.abstractConcreteMapping:
                if len(self.independentMolecules) > 0: 
                    self.abstractConcreteMapping[m] = self.independentMolecules.pop()
                    context.within[0].createMolecule( self.abstractConcreteMapping[m] )
                else:                                                                                                                                                                                                     
                    raise Exception('No more independent molecules')  
                    
        for i, m in enumerate(context.part.precompileMoleculesAfter):
            if not m in self.abstractConcreteMapping:
                if len(self.independentMolecules) > 0:
                    self.abstractConcreteMapping[m] = self.independentMolecules.pop()
                    context.within[0].createMolecule( self.abstractConcreteMapping[m] )
                else:
                    raise Exception('No more independent molecules')

        for i, m in enumerate(context.part.precompileMoleculesAfter):
            context.part.precompileMoleculesAfter[i] = self.abstractConcreteMapping[m]
        for i, m in enumerate(context.part.precompileMoleculesBefore):
            context.part.precompileMoleculesBefore[i] = self.abstractConcreteMapping[m]
            
class AndGate(Circuit):
    moleculeInputOne = None
    moleculeInputTwo = None
    promoterOnePositive = True
    promoterTwoPositive = True
    moleculeOutput = None
    
    def mainCircuit(self):
        self.createMolecule(HrpR)
        self.createMolecule(HrpS)

        if self.promoterOnePositive:
            self.addPart(PositivePromoter(self.moleculeInputOne))
        else:
            self.addPart(NegativePromoter(self.moleculeInputOne))
        
        self.addPart(CodingRegion(HrpR))     
        
        if self.promoterTwoPositive:
            self.addPart(PositivePromoter(self.moleculeInputTwo))
        else:
            self.addPart(NegativePromoter(self.moleculeInputTwo))
        
        self.addPart(CodingRegion(HrpS))
        
        self.addPart(hrpL)
        self.addPart(CodingRegion(self.moleculeOutput))

class Repressilation(Aspect):    
    # weave repressilator into a one-protein coding circuit
    def mainAspect(self):
        self.repressilatorComposite = AbstractRepressilatorComposite()
        
        # promoter is unregulated
        # (not from this aspect or the composite)
        unregulatedPromoterReplace = PointCut(
            PartSignature('!AbstractRepressilatorComposite.Promoter+()') &
            PartSignature('!Repressilation.Promoter+()'),PointCut.REPLACE)
        self.addAdvice(unregulatedPromoterReplace, self.insertRegulatedBySignalPromoter, 10)
        ## regulated case
        #regulatedCodingExpressionReplace = PointCut(
        #(PartSignature('!AbstractRepressilatorComposite.Promoter+(Molecule+)') &
        #PartSignature('!Aspect+.Promoter+(Molecule+)') &
        #PartSignature('!AndGate.Promoter+(Molecule+)')) % 
        #PartSignature('*.CodingRegion+(Molecule+)'), PointCut.REPLACE)
        #self.addAdvice(regulatedCodingExpressionReplace, self.insertBooleanAndSignal, 10)

        # add repressilatior composite
        #promoterAndCodingRegion = PartSignature('Repressilation.Promoter+') % PartSignature('*.RBS+') % PartSignature('*.CodingRegion')
        #afterMainCoding = PointCut(promoterAndCodingRegion | PartSignature('Repressilation.AndGate'),PointCut.AFTER)
        promoterAndCodingRegion = PartSignature('Repressilation.Promoter+') % PartSignature('*.RBS+') % PartSignature('*.CodingRegion')
        afterMainCoding = PointCut(promoterAndCodingRegion,PointCut.AFTER)        
        self.addAdvice(afterMainCoding, self.insertComposite, 10)
        
    def insertRegulatedBySignalPromoter(self,context):
        signal = self.repressilatorComposite.getSignal()
        context.within[0].createMolecule( signal )
        self.addPart(PositivePromoter(signal))

    def insertBooleanAndSignal(self,context):        
        originalpromoter = context.part.before
        originalsignal = originalpromoter.moleculeConnection
        originaloutput = context.part.moleculeConnection
        andGate = AndGate()
        andGate.moleculeInputOne = originalsignal
        if isinstance(context.part,PositivePromoter):
            andGate.setInputOneRepressing = False
        else:
            andGate.setInputOneRepressing = True
        andGate.moleculeInputTwo = self.repressilatorComposite.getSignal()
        andGate.setInputOneRepressing = True
        andGate.moleculeOutput = originaloutput
        self.addPart(andGate)
        
    def insertComposite(self,context):
        print "inserting composite"
        #self.addPart(AbstractRepressilatorComposite) # problem adding this 
        self.addPart(Promoter()) # even this 

# Testing instantiation : all working
#print Weaver(InduceGFP).output()
#print Weaver(CodingGFP).output()
#print Weaver(AbstractRepressilatorComposite).output()
#print Weaver(AbstractRepressilatorComposite,InstantiateAbstractMolecules).output()
#print Weaver(AbstractRepressilatorComposite,DesignRules,InstantiateAbstractMolecules).output()

# repressilator with GFP output (currently not working)
print Weaver(CodingGFP,DesignRules,Repressilation).output()

# repressilator with inducible output
#print Weaver(InduceGFP,Repressilation).output()

# repressilator with GFP output
#print Weaver(CodingGFP,Repressilation,InstantiateAbstractMolecules).output()

# repressilator with inducible output
# print Weaver(InduceGFP,Repressilation,InstantiateAbstractMolecules).output()



