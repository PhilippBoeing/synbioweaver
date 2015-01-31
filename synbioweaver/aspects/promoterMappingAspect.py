from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os, copy

# This helper class associates promoters with downstream coding regions, regulators and polarities
class promoterMapping:

    def __init__(self, prmtr_name, prmtr_scope, regulators, polarities, coding):
        self.prmtr_name = prmtr_name
        self.prmtr_scope = prmtr_scope
        self.regulators = regulators
        self.polarities = polarities
        self.coding = coding

    def getId(self):
        if PromoterMapping.multiComp == True:
            return str(self.prmtr_scope) + "." + str(self.prmtr_name)
        else:
            return str(self.prmtr_name)

    def getRegulators(self):
        if PromoterMapping.multiComp == True:
            return [ str(self.prmtr_scope) + "." + str(x) for x in self.regulators ]
        else:
            return [ str(x) for x in self.regulators ]
            
    def getCodings(self):
        if PromoterMapping.multiComp == True:
            return [ str(self.prmtr_scope) + "." + str(x) for x in self.coding ]
        else:
            return [ str(x) for x in self.coding ]

    def getPolarities(self):
        return self.polarities

class PromoterMapping(Aspect):
    
    def mainAspect(self):
        Reaction.param_counter = 0
        PromoterMapping.builtMap = False
        PromoterMapping.multiComp = False
       
        self.addTypeAdvice(PartSignature('*.*'), self.isPromoter, 'isPromoter')
        self.addTypeAdvice(PartSignature('*.Promoter+'), self.promoterCoding, 'generatePromoterMap')
        self.addWeaverOutput(self.buildPromoterMap)
        
        self.promoterMap = []
   
    def promoterCoding(self, part):
        # set up a structure to hold regulators, type, coding
        prmtr_name = part.__class__.__name__ 
        prmtr_scope = part.scope.circuitName
        
        if not isinstance(part,ConstitutivePromoter):
            regulators = part.getRegulatedBy()
            pols = []
            for regulator in regulators:
                if isinstance(part,NegativePromoter) or (isinstance(part,HybridPromoter) and regulator in part.getRepressors()) :
                    pols.append(-1)
                if isinstance(part,PositivePromoter) or (isinstance(part,HybridPromoter) and regulator in part.getInducers()) :
                    pols.append(+1)
                
            # print prmtr_name, regulator, pol
            self.promoterMap.append( promoterMapping( prmtr_name, prmtr_scope, regulators, pols, [] ) )
        else:
            self.promoterMap.append( promoterMapping( prmtr_name, prmtr_scope, [],[],[] ) )

        # loop over and get all coding regions until at end of parts list or we reach another promoter        
        nextpart = part.getAfterPart()
        while nextpart != None and isinstance(nextpart,Promoter) == False:
            
            if isinstance(nextpart,CodingRegion) == True:
                coding_name = nextpart.getCodingFor()[0]
                self.promoterMap[-1].coding.append( coding_name )
    
            nextpart = nextpart.getAfterPart()
            
    def isPromoter(self, part):
        #print "part", part
        if isinstance(part,Promoter):
            return True

    def buildPromoterMap(self, weaverOutput):
        if PromoterMapping.builtMap == False:
           
            # main circuit
            for part in weaverOutput.partList:
                if part.isPromoter() == True:
                    part.generatePromoterMap()
                    
            # sub circuits
            if len(weaverOutput.wovenCircuitList) > 0:
                PromoterMapping.multiComp = True
                for circ in weaverOutput.wovenCircuitList:
                    for part in circ.partList:
                        if part.isPromoter() == True:
                            part.generatePromoterMap()
            
            #for pr in self.promoterMap:
            #    print pr.getId(), pr.getRegulators(), pr.getCodings(), pr.getPolarities()
                            
            PromoterMapping.builtMap = True 

        return self.promoterMap
        
