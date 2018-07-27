from synbioweaver.core import *
from synbioweaver.aspects.designRulesAspect import *
from instantiateAbstractMoleculesAspect import *
from synbioweaver.aspects.pigeonOutputAspect import *

declareNewMolecule('AHL')
declareNewMolecule('GFP')
declareNewMolecule('mCherry')
declareNewMolecule('CFP')
declareNewMolecule('LacI')
declareNewMolecule('IPTG')
declareNewMolecule('LacI_IPTG')
declareNewMolecule('zero')

declareNewPart('t1',Terminator)
declareNewPart('t2',Terminator)
declareNewPart('t3',Terminator)
declareNewPart('r1',RBS )
declareNewPart('r2',RBS )
declareNewPart('r3',RBS )

declareNewPart('P_lux', PositivePromoter, [AHL] )
declareNewPart('cGFP', CodingRegion, moleculesAfter=[GFP] )

declareNewPart('P_lac', NegativePromoter, [LacI] )
declareNewPart('cLacI', CodingRegion, moleculesAfter=[LacI] )
declareNewPart('Pc', ConstitutivePromoter )
declareNewPart('cmCherry', CodingRegion, moleculesAfter=[mCherry] )

class AHLReceiver(Circuit):
    def mainCircuit(self):
        self.createMolecule(AHL)

        self.addPart(P_lux)
        self.addPart(r1)
        self.addPart(cGFP)
        self.addPart(t1)

class Inverter(Circuit):

    def mainCircuit(self):
        self.createMolecule(IPTG)

        #self.addPart(Pc)
        self.addPart(r2)
        self.addPart(cLacI)
        self.addPart(t2)
        self.addPart(P_lac)
        #self.addPart(r3)
        #self.addPart(cmCherry)
        #self.addPart(t3)
        
        #inducer must come first to display in pigeon properly
        self.reactionFrom(IPTG, LacI) >> self.reactionTo( LacI_IPTG )
        self.reactionFrom(LacI_IPTG) >> self.reactionTo( zero )

class InverterFull(Circuit):

    def mainCircuit(self):
        self.createMolecule(IPTG)

        self.addPart(Pc)
        self.addPart(r2)
        self.addPart(cLacI)
        self.addPart(t2)
        self.addPart(P_lac)
        self.addPart(r3)
        self.addPart(cmCherry)
        self.addPart(t3)
        
        #inducer must come first to display in pigeon properly
        self.reactionFrom(IPTG, LacI) >> self.reactionTo( LacI_IPTG )
        self.reactionFrom(LacI_IPTG) >> self.reactionTo( zero )

class Inversion(Aspect):    
    
    def mainAspect(self):
        self.Inverter = Inverter()

        PosPromoter = PartSignature('AHLReceiver.PositivePromoter+')
        afterPosPromoter = PointCut(PosPromoter,PointCut.AFTER)
        self.addAdvice(afterPosPromoter,self.insertInverter)

    def insertInverter(self,context):
        #print context, context.part
        self.addPart(Inverter())
        

class SwapReporter(Aspect):
    def mainAspect(self):
        Reporter = PartSignature('*.CodingRegion+(GFP)')
        atReporter = PointCut(Reporter,PointCut.REPLACE)
        self.addAdvice(atReporter,self.replaceWithCFP)
        #self.addAdvice(atReporter,self.replaceWithCircuit)

    def replaceWithCFP(self,context):
        self.createMolecule(CFP)
        declareNewPart('cCFP', CodingRegion,  moleculesAfter=[CFP])
        self.addPart(cCFP)
        return True
        
    #def replaceWithCircuit(self,context):
    #    self.addPart(InverterFull())
    #    return True
        



print "AHL receiver"
d1 = Weaver(AHLReceiver, PigeonOutput).output()
print d1.printPigeonOutput()

print "Inverter"
d2 = Weaver(Inverter, PigeonOutput).output()
print d2.printPigeonOutput()

#print "InverterFull"
#d = Weaver(InverterFull, PigeonOutput).output()
#print d.printPigeonOutput()



# Here we use the Inverter circuit as a new part to insert into the AHL receiver
print "Inverted AHL receiver"
d3 = Weaver(AHLReceiver, Inversion, PigeonOutput).output()
print d3.printPigeonOutput()

print "Inverted AHL receiver + GFP -> CFP"
d4 = Weaver(AHLReceiver, Inversion, SwapReporter, PigeonOutput).output()
print d4.printPigeonOutput()

