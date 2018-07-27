from synbioweaver.core import *
from synbioweaver.aspects.pigeonOutputAspect import *

# GFP circuit definitions
declareNewMolecule('AbstractMolecule')
declareNewMolecule('externalInducer')
declareNewMolecule('GFP')
declareNewPart('cGFP', CodingRegion, moleculesAfter=[GFP] )

class InduceGFP(Circuit):
    # code GFP dependent on an external inducer
    def mainCircuit(self):
        self.createMolecule(externalInducer)
        self.addPart(PositivePromoter(externalInducer))
        self.addPart(cGFP)
        
class CodingGFP(Circuit):
    # constitutively code gfp
    def mainCircuit(self):
        self.addPart(ConstitutivePromoter)
        self.addPart(cGFP)
        
class InstantiateAbstractMolecules(Aspect):
    # from a list of non-cross-talking transcription factors
    # instantiate abstract molecules
    def mainAspect(self):
        self.independentMolecules = [declareNewMolecule('TetR'),declareNewMolecule('LacL'),declareNewMolecule('cl')]
        
        self.abstractConcreteMapping = {}
        self.addAdvice(PointCut('*.*(AbstractMolecule+)',PointCut.BEFORE),self.makeConcrete,30)
        
    def makeConcrete(self,context):
        if not context.part.moleculeConnection in self.abstractConcreteMapping:
            if len(self.independentMolecules) > 0:
                self.abstractConcreteMapping[context.part.moleculeConnection] = self.independentMolecules.pop()
            else:
                raise Exception('No more independent molecules')
        
        context.part.moleculeConnection = self.abstractConcreteMapping[context.part.moleculeConnection]
        
class AbstractRepressilatorComposite(Circuit):
    declareNewMolecule('A')
    declareNewMolecule('B')
    declareNewMolecule('C')

    def getSignal(self):
        return A
        
    def mainCircuit(self):
        declareNewPart('p1', NegativePromoter, [C] )
        declareNewPart('p2', NegativePromoter, [A] )
        declareNewPart('p3', NegativePromoter, [B] )
        declareNewPart('cA', CodingRegion, moleculesAfter=[A] )
        declareNewPart('cB', CodingRegion, moleculesAfter=[B] )
        declareNewPart('cC', CodingRegion, moleculesAfter=[C] )
        self.addPart( p1 )
        self.addPart( cA )
        self.addPart( p2 )
        self.addPart( cB )
        self.addPart( p3 )
        self.addPart( cC )
   
class AndGate(Circuit):
    moleculeInputOne = None
    moleculeInputTwo = None
    promoterOnePositive = True
    promoterTwoPositive = True
    moleculeOutput = None
    
    def mainCircuit(self):
        # declare internal andGate parts / molecules - doing this here has advantage:
        # will throw warning if already used, i.e. if cross-talk potential
        
        # AND gate design based on Wang et al, 2011, 
        # Engineering modular and orthogonal genetic logic gates for robust digital-like synthetic biology
        # Nature communications

        declareNewMolecule('HrpR')
        declareNewMolecule('HrpS')
        declareNewMolecule('HrpSandR') # not ideal...
        declareNewPart('hrpL', PositivePromoter, [HrpSandR] ) # not ideal...

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
        # regulated case
        #regulatedCodingExpressionReplace = PointCut(
        #(PartSignature('!AbstractRepressilatorComposite.Promoter+(Molecule+)') &
        #PartSignature('!Aspect+.Promoter+(Molecule+)') & 
        #PartSignature('!AndGate.Promoter+(Molecule+)')) % PartSignature('*.CodingRegion+(Molecule+)'), PointCut.REPLACE)
        #self.addAdvice(regulatedCodingExpressionReplace, self.insertBooleanAndSignal, 10)

        # add repressilatior composite
        #promoterAndCodingRegion = PartSignature('Repressilation.Promoter+') % PartSignature('*.RBS+') % PartSignature('*.CodingRegion')
        #afterMainCoding = PointCut(promoterAndCodingRegion | PartSignature('Repressilation.AndGate'),PointCut.AFTER)

        promoterAndCodingRegion = PartSignature('Repressilation.Promoter+') % PartSignature('*.RBS+') % PartSignature('*.CodingRegion+')
        afterMainCoding = PointCut(promoterAndCodingRegion,PointCut.AFTER)

        self.addAdvice(afterMainCoding, self.insertComposite, 10)
        
    def insertRegulatedBySignalPromoter(self,context):
        signal = self.repressilatorComposite.getSignal()
        context.within[0].createMolecule( signal )
        declareNewPart('ppA', PositivePromoter, [signal] )
        self.addPart(ppA)
        return True

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
        return True

    def insertComposite(self,context):
        print "insertComposite"
        print context, context.part
        #self.addPart(self.repressilatorComposite)
        context.within[0].createMolecule( A )
        context.within[0].createMolecule( B )
        context.within[0].createMolecule( C )
        self.addPart(AbstractRepressilatorComposite)
        return True
        
class DesignRules(Aspect):
    beforeCodingRegion = PointCut('*.CodingRegion+',PointCut.BEFORE)
    afterCodingRegion = PointCut('*.CodingRegion+',PointCut.AFTER)
    
    def insertRBS(self,context):
        self.addPart(RBS)
    
    def insertTerminator(self,context):
        self.addPart(Terminator)
        
    def mainAspect(self):
        self.addAdvice(self.beforeCodingRegion,self.insertRBS,1)
        self.addAdvice(self.afterCodingRegion,self.insertTerminator,1)
        
# repressilator with GFP output
print "\n\nRUNNING"
#print Weaver(CodingGFP,DesignRules).output()
out = Weaver(CodingGFP,Repressilation,DesignRules,PigeonOutput).output() #.printPigeonOutput()
print out
print out.printPigeonOutput()

# repressilator with inducible output
#print Weaver(InduceGFP,Repressilation,DesignRules).output()

# repressilator with GFP output
#print Weaver(CodingGFP,Repressilation,DesignRules,InstantiateAbstractMolecules).output()

# repressilator with inducible output
#print Weaver(InduceGFP,Repressilation,DesignRules,InstantiateAbstractMolecules).output()
