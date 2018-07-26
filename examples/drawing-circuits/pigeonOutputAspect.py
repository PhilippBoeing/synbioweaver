from synbioweaver.core import *

class PigeonOutput(Aspect):

    def mainAspect(self):
        self.addWeaverOutput(self.printPigeonOutput)

    def printPigeonOutput(self,weaverOutput):

        print "\nPigeon Output:\n"

        for p in weaverOutput.partList:
            name = p.__class__.__name__
            if( isinstance(p,Promoter) ):
                print 'p', name, '1'
            elif( isinstance(p,RBS) ):
                print 'r', name, '1'
            elif( isinstance(p,CodingRegion) ):
                print 'c', name, '1'
            elif( isinstance(p,Terminator) ):
                print 't', name, '1'

        print "# Arcs"
        
        for pr in weaverOutput.partList:
            name = pr.__class__.__name__

            if( isinstance(pr,NegativePromoter) ):
                for neg in pr.getRegulatedBy():
                    if len(neg.getBeforeNodes(Part)) > 0:
                        print neg.getBeforeNodes(Part)[0].__class__.__name__,  "rep", name   
                    else:
                        print neg.__class__.__name__, "rep", name

            if( isinstance(pr,PositivePromoter) ):
                for pos in pr.getRegulatedBy():
                    if len(pos.getBeforeNodes(Part)) > 0:
                        print pos.getBeforeNodes(Part)[0].__class__.__name__,  "ind", name
                    else:
                        print pos.__class__.__name__, "ind", name

            if( isinstance(pr,HybridPromoter) ):
                for neg in pr.getRepressors():
                    if len(neg.getBeforeNodes(Part)) > 0:
                        print neg.getBeforeNodes(Part)[0].__class__.__name__,  "rep", name
                    else:
                        print neg.__class__.__name__, "rep", name

                for pos in pr.getInducers():
                    if len(pos.getBeforeNodes(Part)) > 0:
                        print pos.getBeforeNodes(Part)[0].__class__.__name__,  "ind", name
                    else:
                        print pos.__class__.__name__, "ind", name

        print "\n"

