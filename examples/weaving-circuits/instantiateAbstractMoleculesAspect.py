from synbioweaver.core import *

import copy

class InstantiateAbstractMolecules(Aspect):
    # from a list of non-cross-talking transcription factors
    # instantiate actual molecules in an abstract gene regulatory network

    def __init__(self, repressors):
        super(InstantiateAbstractMolecules,self).__init__()
        self.independentMolecules = copy.copy(repressors)

    def mainAspect(self):
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
