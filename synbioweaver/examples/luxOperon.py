from aosb import *
from reactionNetworksAspect import *

# this example is based on the inducible lux quorum sensing system




class LuxOperon(Circuit):
    def mainCircuit(self):

        # constitutive promoter
        #declareNewMolecule('GFP')
        #declareNewPart('Pc', ConstitutivePromoter)
        #self.addPart(Pc)
        #self.addPart(CodingRegion(GFP))

        # positve feedback
        declareNewMolecule('AHL')
        declareNewMolecule('AHLn2')
        declareNewMolecule('AHLn4')
        declareNewMolecule('R')
        declareNewMolecule('R_AHLn4')
        declareNewMolecule('GFP')

        self.createMolecule(AHL) 
        self.createMolecule(AHLn2) 
        self.createMolecule(AHLn4) 
        self.createMolecule(R)
        self.createMolecule(R_AHLn4)
        
        # positive promoters
        declareNewPart('Plux', PositivePromoter, [R_AHLn4])
        self.addPart(Plux)
        self.addPart(RBS)
        self.addPart(CodingRegion(R))
        self.addPart(CodingRegion(AHL))
        self.addPart(Terminator)
        
        declareNewPart('Pc', NegativePromoter, [R])
        self.addPart(Pc)
        self.addPart(CodingRegion(GFP))
        self.addPart(Terminator)

        self.reactionFrom(AHL,AHL) >> self.reactionTo( AHLn2 )
        self.reactionFrom(AHLn2,AHLn2) >> self.reactionTo( AHLn4 )
        self.reactionFrom(R,AHLn4) >> self.reactionTo( R_AHLn4 )


# Take the compiled design and add Type Advice that generates a mass action model
compiledDesign1 = Weaver(LuxOperon, ReactionNetworks, MassActionKineticsRNA, PrintReactionNetwork).output()
compiledDesign1.printReactionNetwork()

compiledDesign1 = Weaver(LuxOperon, ReactionNetworks, SheaAckersKineticsRNA, PrintReactionNetwork).output()
compiledDesign1.printReactionNetwork()
