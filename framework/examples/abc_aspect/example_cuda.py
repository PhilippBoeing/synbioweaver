from aosb import *
import abcOutput

"""
STILL TO DO:
There can only be one reaction per species now

I dont know how to include the reacted molecules in the reactions list. so must be added manually.
This is not good as if the reacted molecule is named differently by the user it will not

Positive promoter
Degradation reactions?

"""

fit = 'species7 species8'
epsilon = 1
particles = 100
alpha = 0.9
#simtype = 'Gillespie'
prior_distribution = ['constant', 'constant', 'constant', 'constant', 'constant', 'constant', 'constant', 'constant', 'constant', 'constant']
init_cond__distribution = ['constant', 'constant', 'constant', 'constant', 'constant', 'constant', 'uniform', 'uniform', 'constant', 'constant']
#The same order  as the reactions capable from each part. Dimerisation last. First in compartment so always constant 1
priors = [0.4, 0.05, 0.4, 0.05, 0.4, 0.05, 0.4, 0.05, 0.005, 0.005]
#The same order as the order of self.createMolecule
initial_conditions = [0, 0, 0, 0, 1, 1, [0, 10], [0, 10], 0, 0]
#abcOutput.setInputVariables(epsilon, fit, particles, alpha, prior_distribution, priors)


class SimpleSwitch(Circuit):
    def mainCircuit(self):

        declareNewMolecule('TetR')
        declareNewMolecule('LacI')
        declareNewMolecule('mCherry')
        declareNewMolecule('GFP')
        declareNewMolecule('TetR_TetR')
        declareNewMolecule('LacI_LacI')
        declareNewPart('PLtetO', NegativePromoter, [TetR_TetR])
        declareNewPart('Ptrc_2', NegativePromoter, [LacI_LacI])
        declareNewMolecule('PLtetO_TetR_TetR')
        declareNewMolecule('Ptrc_2_LacI_LacI')

        self.createMolecule(TetR)
        self.createMolecule(LacI)
        self.createMolecule(mCherry)
        self.createMolecule(GFP)
        self.createMolecule(Ptrc_2)
        self.createMolecule(PLtetO)
        self.createMolecule(TetR_TetR)
        self.createMolecule(LacI_LacI)
        self.createMolecule(PLtetO_TetR_TetR)
        self.createMolecule(Ptrc_2_LacI_LacI)

        self.addPart(Ptrc_2)
        self.addPart(CodingRegion(TetR))
        self.addPart(Ptrc_2)
        self.addPart(CodingRegion(mCherry))
        self.addPart(PLtetO)
        self.addPart(CodingRegion(LacI))
        self.addPart(PLtetO)
        self.addPart(CodingRegion(GFP))

        self.reactionFrom(PLtetO,TetR) >> self.reactionTo(PLtetO_TetR_TetR)
        self.reactionFrom(Ptrc_2,LacI) >> self.reactionTo(Ptrc_2_LacI_LacI)
        self.reactionFrom(TetR,TetR) >> self.reactionTo(TetR_TetR)
        self.reactionFrom(LacI,LacI) >> self.reactionTo(LacI_LacI)

compiledDesign = Weaver(SimpleSwitch, abcOutput.CudaOutput).output()
print compiledDesign.cudaOutput()
