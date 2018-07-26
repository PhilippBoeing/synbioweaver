from synbioweaver.core import *

class PrintStack(Aspect):
    def mainAspect(self):
        self.addAdvice(PointCut('*.*',PointCut.AFTER),self.printStack)
    
    def printStack(self,context):
        print ','.join(str(context.within[i]) for i in range(len(context.within))) + " add: " + str(context.part)
