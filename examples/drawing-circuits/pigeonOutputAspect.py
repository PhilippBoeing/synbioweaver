from synbioweaver.core import *

#http://code.activestate.com/recipes/578948-flattening-an-arbitrarily-nested-list-in-python/
def flatten(lis):
    """Given a list, possibly nested to any level, return it flattened."""
    new_lis = []
    for item in lis:
        if type(item) == type([]):
            new_lis.extend(flatten(item))
        else:
            new_lis.append(item)
    return new_lis

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
        # this will store a list of arrays [mol, arcfrom, arcto]
        arcs = []
        for pr in weaverOutput.partList:
            name = pr.__class__.__name__

            if( isinstance(pr,NegativePromoter) ):
                for neg in pr.getRegulatedBy():
                    if len(neg.getBeforeNodes(Part)) > 0:
                        print neg.getBeforeNodes(Part)[0].__class__.__name__,  "rep", name   
                        arcs.append( [neg.__class__.__name__, neg.getBeforeNodes(Part)[0].__class__.__name__, name] )
                    else:
                        print neg.__class__.__name__, "rep", name

            if( isinstance(pr,PositivePromoter) ):
                for pos in pr.getRegulatedBy():
                    if len(pos.getBeforeNodes(Part)) > 0:
                        print pos.getBeforeNodes(Part)[0].__class__.__name__,  "ind", name
                        arcs.append( [pos.__class__.__name__,pos.getBeforeNodes(Part)[0].__class__.__name__, name] )
                    else:
                        print pos.__class__.__name__, "ind", name

            if( isinstance(pr,HybridPromoter) ):
                for neg in pr.getRepressors():
                    if len(neg.getBeforeNodes(Part)) > 0:
                        print neg.getBeforeNodes(Part)[0].__class__.__name__,  "rep", name
                        arcs.append( [neg.__class__.__name__,neg.getBeforeNodes(Part)[0].__class__.__name__, name] )
                    else:
                        print neg.__class__.__name__, "rep", name

                for pos in pr.getInducers():
                    if len(pos.getBeforeNodes(Part)) > 0:
                        print pos.getBeforeNodes(Part)[0].__class__.__name__,  "ind", name
                        arcs.append( [pos.__class__.__name__,pos.getBeforeNodes(Part)[0].__class__.__name__, name] )
                    else:
                        print pos.__class__.__name__, "ind", name

        # Look for inducer molecule regulations
        #print "arcs:", arcs

        for mol in weaverOutput.moleculeList:
            #print mol, mol.after, mol.before
            
            # skip reactions involving promoters
            if len(mol.after) > 0:
                if isinstance(mol.after[0],Promoter):
                    continue
            
            flat_list = [mol] + flatten(mol.after) + flatten(mol.before)
            flat_list = [str(x) for x in flat_list]
            #print flat_list

            for a in arcs:
                if a[0] in flat_list:  
                    #print "\tfound:", a, flat_list
                    if 'zero' in flat_list:
                        print flat_list[2], 'rep', a[1]+"-"+a[2]
                    else:
                        print flat_list[1], 'ind', a[1]+"-"+a[2]

        print "\n"

