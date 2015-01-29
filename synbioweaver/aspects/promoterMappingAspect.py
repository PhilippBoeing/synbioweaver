from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os, copy

# This helper class associates promoters with downstream coding regions, regulators and polarities
class promoterMapping:

    def __init__(self, prmtr_name, regulators, polarities, coding):
        self.prmtr_name = prmtr_name
        self.regulators = regulators
        self.polarities = polarities
        self.coding = coding

class PromoterMapping(Aspect):
    
    def mainAspect(self):
        Reaction.param_counter = 0
        PromoterMapping.builtMap = False
        # The goal of this aspect is to parse the design and generate a set of reactions, rates and parameters
        self.addTypeAdvice(PartSignature('*.*'), self.isPromoter, 'isPromoter')
        self.addTypeAdvice(PartSignature('*.Promoter+'), self.promoterCoding, 'generatePromoterMap')

        self.addWeaverOutput(self.buildPromoterMap)
        
        self.promoterMap = []
        self.circuitMap = {}

    def promoterCoding(self, part):
        # set up a structure to hold regulators, type, coding
        #print part
        prmtr_name = part.__class__.__name__ 
        
        # if more than one circuit, append the cicuit name
        if len(self.circuitMap) > 0:
            prmtr_name = self.circuitMap[prmtr_name] + "." + prmtr_name

        if isinstance(part,NegativePromoter):
            regulator = part.getRegulatedBy()[0]
            self.promoterMap.append( promoterMapping( prmtr_name, [regulator],[-1],[] ) )
        elif isinstance(part,PositivePromoter):
            regulator = part.getRegulatedBy()[0]
            self.promoterMap.append( promoterMapping( prmtr_name, [regulator],[1],[] ) )
        else:
            self.promoterMap.append( promoterMapping( prmtr_name, [None],[],[] ) )
    
        # loop over and get all coding regions until at end of parts list or we reach another promoter        
        nextpart = part.getAfterPart()
        while nextpart != None and isinstance(nextpart,Promoter) == False:
            
            if isinstance(nextpart,CodingRegion) == True:
                #print "\t", nextpart, "Coding"
                coding_name = nextpart.getCodingFor()[0]
                
                if len( self.circuitMap ) > 0:
                    full_coding_name = self.circuitMap[ part.__class__.__name__  ] + "." + str(coding_name)
                    self.promoterMap[-1].coding.append( full_coding_name )
                else:
                    self.promoterMap[-1].coding.append( coding_name )
    
            nextpart = nextpart.getAfterPart()
            
    def isPromoter(self, part):
        #print "part", part
        if isinstance(part,Promoter):
            return True

    def buildPromoterMap(self, weaverOutput):
        if PromoterMapping.builtMap == False:
           
            # if this is zero we are dealing with one compartment
            if len(weaverOutput.wovenCircuitList) == 0:
                for part in weaverOutput.partList:
                    if part.isPromoter() == True:
                        #self.circuitMap[ part.__class__.__name__ ] = "main"
                        part.generatePromoterMap()
            
            else:
                for circ in weaverOutput.wovenCircuitList:
                    for part in circ.partList:
                        if part.isPromoter() == True:
                            self.circuitMap[ part.__class__.__name__ ] = circ.circuitName
                            part.generatePromoterMap()
            
            #print self.circuitMap
            #print self.promoterMap
            PromoterMapping.builtMap = True 

        return [self.promoterMap, self.circuitMap]
        
