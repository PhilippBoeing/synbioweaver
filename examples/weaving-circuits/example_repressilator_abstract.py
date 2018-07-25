from synbioweaver.core import *

# GFP circuit definitions
declareNewMolecule('AbstractMolecule')
declareNewMolecule('externalInducer')
declareNewMolecule('GFP')
declareNewPart('P1', PositivePromoter, [externalInducer])
declareNewPart('Pc', ConstitutivePromoter)

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

class PrintStack(Aspect):
    def mainAspect(self):
        self.addAdvice(PointCut('*.*',PointCut.AFTER),self.printStack)
    
    def printStack(self,context):
        print ','.join(str(context.within[i]) for i in range(len(context.within))) + " add: " + str(context.part)

class InduceGFP(Circuit):
    # code GFP dependent on an external inducer
    def mainCircuit(self):
        self.createMolecule(externalInducer)
        self.createMolecule(GFP)
        self.addPart(P1)
        self.addPart(CodingRegion(GFP))
        
class CodingGFP(Circuit):
    # constitutively code gfp
    def mainCircuit(self):
        self.createMolecule(GFP)
        self.addPart(Pc)
        self.addPart(CodingRegion(GFP))
        
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
    # an abstract repressilator circuit
   
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

        #self.exportMolecule(A)
            
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
        # regulated case
        regulatedCodingExpressionReplace = PointCut(
        (PartSignature('!AbstractRepressilatorComposite.Promoter+(Molecule+)') &
        PartSignature('!Aspect+.Promoter+(Molecule+)') &
        PartSignature('!AndGate.Promoter+(Molecule+)')) % 
        PartSignature('*.CodingRegion+(Molecule+)'), PointCut.REPLACE)
        self.addAdvice(regulatedCodingExpressionReplace, self.insertBooleanAndSignal, 10)

        # add repressilatior composite
        promoterAndCodingRegion = PartSignature('Repressilation.Promoter+') % PartSignature('*.RBS+') % PartSignature('*.CodingRegion')
        afterMainCoding = PointCut(promoterAndCodingRegion | PartSignature('Repressilation.AndGate'),PointCut.AFTER)
        self.addAdvice(afterMainCoding, self.insertComposite, 10)
        
    def insertRegulatedBySignalPromoter(self,context):
        declareNewPart('Pg', PositivePromoter, [ self.repressilatorComposite.getSignal() ])
        self.addPart(Pg)
        
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
        self.addPart(AbstractRepressilatorComposite)
        
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

#Weaver(InduceGFP,PrintStack)
#Weaver(CodingGFP,PrintStack)
#Weaver(AbstractRepressilatorComposite,PrintStack)
#Weaver(AndGate,PrintStack)
        
# repressilator with GFP output
#print Weaver(CodingGFP,Repressilation).output()

# repressilator with inducible output
#print Weaver(InduceGFP,Repressilation,DesignRules).output()

# repressilator with GFP output
#print Weaver(CodingGFP,Repressilation,DesignRules,InstantiateAbstractMolecules).output()

# repressilator with inducible output
# print Weaver(InduceGFP,Repressilation,DesignRules,InstantiateAbstractMolecules).output()



