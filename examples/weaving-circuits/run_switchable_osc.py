from synbioweaver.core import *
from synbioweaver.aspects.pigeonOutputAspect import *
from synbioweaver.aspects.printStackAspect import *
from synbioweaver.aspects.designRulesAspect import *

# GFP circuit definitions
declareNewMolecule('exIn')
declareNewMolecule('GFP')
declareNewPart('cGFP', CodingRegion, moleculesAfter=[GFP] )

class CodingGFP(Circuit):
    # constitutively code gfp
    def mainCircuit(self):
        declareNewPart('Pconst',ConstitutivePromoter)
        self.addPart(Pconst)
        self.addPart(cGFP)

class InduceGFP(Circuit):
    # code GFP dependent on an external inducer
    
    def getSignal(self):
        return exIn

    def getOutputName():
        return cGFP

    def mainCircuit(self):
        self.createMolecule(exIn)
        declareNewPart('Pin', PositivePromoter, [exIn] )
        self.addPart(Pin)
        self.addPart(cGFP)
                
class AbstractRepressilatorComposite(Circuit):
    def getSignal(self):
        return A
        
    def mainCircuit(self):
        declareNewMolecule('A')
        declareNewMolecule('B')
        declareNewMolecule('C')
        self.createMolecule(A)
        self.createMolecule(B)
        self.createMolecule(C)

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
        declareNewPart('hrpL', HybridPromoter, [HrpS, HrpR], regulatorInfoMap={HrpS:True,HrpR:True} )

        if self.promoterOnePositive:
            declareNewPart('pos1', PositivePromoter, [self.moleculeInputOne] )
            self.addPart(pos1)
        else:
            declareNewPart('neg1', NegativePromoter, [self.moleculeInputOne] )
            self.addPart(neg1)
       
        declareNewPart('cHrpR', CodingRegion, moleculesAfter=[HrpR] ) 
        self.addPart(cHrpR) 
        
        if self.promoterTwoPositive:
            declareNewPart('pos2', PositivePromoter, [self.moleculeInputTwo] )
            self.addPart(pos2)
        else:
            declareNewPart('neg2', NegativePromoter, [self.moleculeInputTwo] )
            self.addPart(neg2)

        declareNewPart('cHrpS', CodingRegion, moleculesAfter=[HrpS] )
        self.addPart(cHrpS)

        self.addPart(hrpL)
        self.addPart(cGFP)

class ConcreteAndGate(Circuit):
    def mainCircuit(self):
        # same AND gate but a concrete version for plotting

        declareNewMolecule('HrpR')
        declareNewMolecule('HrpS')
        declareNewPart('hrpL', HybridPromoter, [HrpS, HrpR], regulatorInfoMap={HrpS:True,HrpR:True} )

        declareNewMolecule('input1')
        declareNewMolecule('input2')
        self.createMolecule(input1)
        self.createMolecule(input2)

        declareNewPart('pos1', PositivePromoter, [input1] )
        self.addPart(pos1)
        declareNewPart('cHrpR', CodingRegion, moleculesAfter=[HrpR] ) 
        self.addPart(cHrpR) 
        
        declareNewPart('pos2', PositivePromoter, [input2] )
        self.addPart(pos2)
        declareNewPart('cHrpS', CodingRegion, moleculesAfter=[HrpS] )
        self.addPart(cHrpS)

        self.addPart(hrpL)
        self.addPart(cGFP)

class Repressilation(Aspect):    
    # weave repressilator into a one-protein coding circuit
    def mainAspect(self):
        self.repressilatorComposite = AbstractRepressilatorComposite()
        
        # promoter is unregulated (not from this aspect or the composite)
        unregulatedPromoterReplace = PointCut(
            PartSignature('!AbstractRepressilatorComposite.Promoter+()') &
            PartSignature('!Repressilation.Promoter+()'),PointCut.REPLACE)
        self.addAdvice(unregulatedPromoterReplace, self.insertRegulatedBySignalPromoter, 10)
        
        # regulated case
        regulatedCodingExpressionReplace = PointCut(
            (PartSignature('InduceGFP.PositivePromoter+') &
             PartSignature('!AbstractRepressilatorComposite.Promoter+(Molecule+)') & 
             PartSignature('!Aspect+.Promoter+(Molecule+)') & 
             PartSignature('!AndGate.Promoter+(Molecule+)')) % PartSignature('*.CodingRegion+(Molecule+)'), 
            PointCut.REPLACE)
        self.addAdvice(regulatedCodingExpressionReplace, self.insertBooleanAndSignal, 10)

        # add repressilatior composite
        promoterAndCodingRegion = PartSignature('Repressilation.Promoter+') % PartSignature('*.RBS+') % PartSignature('*.CodingRegion+') % PartSignature('*.Terminator+') 
        andGateRegion = PartSignature('Repressilation.AndGate')
        afterCircuit = PointCut(andGateRegion | promoterAndCodingRegion, PointCut.AFTER)
        self.addAdvice(afterCircuit, self.insertComposite, 10)

    def insertRegulatedBySignalPromoter(self,context):
        #print "Fired! Simple system"
        
        signal = self.repressilatorComposite.getSignal()
        context.within[0].createMolecule( signal )
        declareNewPart('ppA', PositivePromoter, [signal] )
        
        self.addPart(ppA)
        return True

    def insertBooleanAndSignal(self,context):        
        #print "Fired! Boolean And Signal"
        #print "insertBooleanAndSignal:", context.part, context.part.before,  context.part.before[0].getRegulatedBy(), context.part.getCodingFor()
        
        originalpromoter = context.part.before[0]
        originalsignal = context.part.before[0].getRegulatedBy()[0]
        originaloutput = context.part.getCodingFor()[0]

        andGate = AndGate()
        
        # Creat first input for AND gate    
        context.within[0].createMolecule( exIn )
        andGate.moleculeInputOne = exIn 

        if isinstance(originalpromoter,PositivePromoter):
            andGate.setInputOneRepressing = False
        else:
            andGate.setInputOneRepressing = True

        # Create second input for AND gate
        in2 = self.repressilatorComposite.getSignal()
        context.within[0].createMolecule( in2 )
        andGate.moleculeInputTwo = in2
       
        andGate.setInputOneRepressing = True
        andGate.moleculeOutput = originaloutput

        self.addPart(andGate)
        return True

    def insertComposite(self,context):
        #print "insertComposite:", context, context.part
        context.within[0].createMolecule( A )
        context.within[0].createMolecule( B )
        context.within[0].createMolecule( C )
        self.addPart(self.repressilatorComposite)
        return True
                
print "\n\n"
out = Weaver(CodingGFP,DesignRules,PigeonOutput).output()
print out
print out.printPigeonOutput()

out = Weaver(InduceGFP,DesignRules,PigeonOutput).output()
print out
print out.printPigeonOutput()

out = Weaver(AbstractRepressilatorComposite,DesignRules,PigeonOutput).output()
print out
print out.printPigeonOutput()

#out = Weaver(ConcreteAndGate,DesignRules,PigeonOutput).output()
#print out
#print out.printPigeonOutput()


print "\n\nRUNNING UNINDUCED"
out = Weaver(CodingGFP,DesignRules,Repressilation,PrintStack,PigeonOutput).output()
print out
print out.printPigeonOutput()

# repressilator with inducible output
print "\n\nRUNNING INDUCED"
out = Weaver(InduceGFP,DesignRules,Repressilation,PrintStack,PigeonOutput).output()
print out
print out.printPigeonOutput()
