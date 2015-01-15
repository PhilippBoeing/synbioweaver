from core import Aspect,RBS,Terminator,PointCut


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