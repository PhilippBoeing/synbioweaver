from synbioweaver.core import *
import pysb
from pysb import kappa
import matplotlib.pyplot as plt

class PYSBmodel(Aspect):
    observableProteins = []

    def __init__(self,rnapCount = 700,ribosomeCount = 18000,operonCount = 1,
                 rnapBindingHigh = 0.0007,rnapBindingLow=7e-07,
                 transcriptionRate = 10,transcriptionFactorBinding = 0.01,
                 transcriptionFactorUnbinding = 0.09,transcriptionFactorUnbindingSole=2.24,
                 transcriptionFactorFactorDegradation=0.00115,rbsBinding = 0.000166,
                 translationInitiationRate=0.167,translationRate=10,ribosomeFalloff = 0.01,
                 rnapFalloffRate=1,rnaDegradation=0.0058):
        self.rnapCount = rnapCount
        self.ribosomeCount = ribosomeCount
        self.operonCount = operonCount
        self.rnapBindingHigh = rnapBindingHigh
        self.rnapBindingLow = rnapBindingLow
        self.transcriptionRate = transcriptionRate
        self.transcriptionFactorBinding = transcriptionFactorBinding
    def mainAspect(self):
        self.partIdCounter = 0

        self.addWeaverOutput(self.runKappaSimulation)

        self.addTypeAdvice(PartSignature('*.NegativePromoter+'),self.pysbPartRulesNegPromoter,'pysbPartRules')
        self.addTypeAdvice(PartSignature('*.RBS+'),self.pysbPartRulesRBS,'pysbPartRules')
        self.addTypeAdvice(PartSignature('*.CodingRegion+'),self.pysbPartRulesCodingRegion,'pysbPartRules')
        self.addTypeAdvice(PartSignature('*.Terminator+'),self.pysbPartRulesTerminator,'pysbPartRules')
        self.addTypeAdvice(PartSignature('*.*'),self.initForPysb,'initForPysb')
        self.addTypeAdvice(MoleculeSignature('*.*'),self.initMoleculeForPysb,'initForPysb')



        # setup model
        self.model = pysb.Model()
        pysb.SelfExporter.target_globals = self.__dict__;

        # Transcription = Transcription and Translation...

        # using parameters from kappa repressilator model

        # Transcription = Transcription and Translation...
        pysb.Monomer('RNAP',['dna','rna'])
        pysb.Monomer('Ribosome',['rna'])

        pysb.Parameter('RNAPcount',700)
        pysb.Parameter('RibosomeCount',18000)
        pysb.Initial(self.RNAP(dna = None,rna = None),self.RNAPcount)
        pysb.Initial(self.Ribosome(rna = None),self.RibosomeCount)
        pysb.Parameter('RNAPBindingHigh', 0.0007)
        pysb.Parameter('RNAPBindingLow', 7e-07)
        pysb.Parameter('OperonCount',1)
        pysb.Parameter('TranscriptionRate',10)
        pysb.Parameter('TranscriptionFactorBinding',0.01)
        pysb.Parameter('TranscriptionFactorUnbinding',0.09)
        pysb.Parameter('TranscriptionFactorUnbindingSole',2.24)
        pysb.Parameter('TranscriptionFactorDegradation',0.00115)
        pysb.Parameter('RBSBinding',0.000166)
        pysb.Parameter('TranslationInitiationRate',0.167)
        pysb.Parameter('TranslationRate',10)
        pysb.Parameter('RibosomeFalloff',0,01)
        pysb.Parameter('RNAPfalloffrate',1)
        pysb.Parameter('RNAdegradation',0.0058)

        self.dnaPartTypeList = []
        self.rnaPartTypeList = []




    def initForPysb(self,part):
        part.pysbname = part.__class__.__name__
        # dna part type
        if not isinstance(part, Promoter):
            self.dnaPartTypeList.append(part.pysbname)
        else:
            # a promoter has 4 dna subsections
            for i in range(1,5):
                self.dnaPartTypeList.append(part.pysbname+'p'+str(i))

        # rna part type
        if not isinstance(part,Terminator):
            self.rnaPartTypeList.append(part.pysbname)

    def initMoleculeForPysb(self,molecule):
        molecule.pysbmonomer = pysb.Monomer(str(molecule),['dna'])



    def pysbPartRulesNegPromoter(self,part):
        regulatorMol = part.getBeforeNodes(Molecule)[0]

        # transcription factor binding
        pysb.Rule('TF_Binding_to_' + part.pysbname + 'p2' + '_noTF',
                self.DNA(binding = None,partType = part.pysbname + 'p3',upstream = 2) +
                regulatorMol.pysbmonomer(dna = None) +
                self.DNA(downstream = 2,binding = None,partType = part.pysbname + 'p2')
                <>
                self.DNA(binding = None,partType = part.pysbname + 'p3',upstream = 3) +
                regulatorMol.pysbmonomer(dna = 1) +
                self.DNA(downstream = 3,binding = 1,partType = part.pysbname + 'p2'),
                self.TranscriptionFactorBinding, self.TranscriptionFactorUnbindingSole)

        pysb.Rule('TF_Binding_to_' + part.pysbname + 'p2' + '_TF',
                self.DNA(binding = 1,partType = part.pysbname + 'p3',upstream = 2) +
                regulatorMol.pysbmonomer(dna = 1) +
                self.DNA(downstream = 2,binding = None,partType = part.pysbname + 'p2') +
                regulatorMol.pysbmonomer(dna = None)
                <>
                self.DNA(binding = 2,partType = part.pysbname + 'p3',upstream = 3) +
                regulatorMol.pysbmonomer(dna = 2) +
                self.DNA(downstream = 3,binding = 1,partType = part.pysbname + 'p2') +
                regulatorMol.pysbmonomer(dna = 1),
                self.TranscriptionFactorBinding,self.TranscriptionFactorUnbinding)

        pysb.Rule('TF_Binding_to_' + part.pysbname + 'p3' + '_noTF',
                self.DNA(binding = None,partType = part.pysbname + 'p3',upstream = 2) +
                regulatorMol.pysbmonomer(dna = None) +
                self.DNA(downstream = 2,binding = None,partType = part.pysbname + 'p2')
                <>
               self.DNA(binding = 1,partType = part.pysbname + 'p3',upstream = 3) +
                regulatorMol.pysbmonomer(dna = 1) +
                self.DNA(downstream = 3,binding = None,partType = part.pysbname + 'p2'),
                self.TranscriptionFactorBinding,self.TranscriptionFactorUnbindingSole)

        pysb.Rule('TF_Binding_to_' + part.pysbname + 'p3' + '_TF',
                self.DNA(binding = None,partType = part.pysbname + 'p3',upstream = 2) +
                regulatorMol.pysbmonomer(dna = 1) +
                self.DNA(downstream = 2,binding = 1,partType = part.pysbname + 'p2') +
                regulatorMol.pysbmonomer(dna = None)
                <>
                self.DNA(binding = 1,partType = part.pysbname + 'p3',upstream = 3) +
                regulatorMol.pysbmonomer(dna = 2) +
                self.DNA(downstream =3,binding = 2,partType = part.pysbname + 'p2') +
                regulatorMol.pysbmonomer(dna = 1),
                self.TranscriptionFactorBinding,self.TranscriptionFactorUnbinding)

        # RNAP binding
        pysb.Rule('RNAP_Binding_To_' + part.pysbname+'_no_TF',
                self.DNA(binding = None,partType = part.pysbname + 'p3',upstream = 2, downstream = 1) +
                self.DNA(binding = None,partType = part.pysbname + 'p4',upstream = 1) +
                self.RNAP(dna = None,rna = None) +
                self.DNA(binding = None,partType = part.pysbname + 'p2',downstream = 2)
                >>
               self.DNA(binding = None,partType = part.pysbname + 'p3',upstream = 3, downstream = 1) +
                self.DNA(binding = 2,partType = part.pysbname + 'p4',upstream = 1) +
                self.RNAP(dna = 2,rna = None) +
                self.DNA(binding = None,partType = part.pysbname + 'p2',downstream = 3),
                self.RNAPBindingHigh)

        pysb.Rule('RNAP_Binding_To_' + part.pysbname+'_TF_on_p2',
                self.DNA(binding = None,partType = part.pysbname + 'p3',upstream = 2,downstream = 1) +
                self.DNA(binding = None,partType = part.pysbname + 'p4',upstream = 1) +
                self.RNAP(dna = None,rna = None) +
                self.DNA(binding = 3,partType = part.pysbname + 'p2',downstream = 2) +
                regulatorMol.pysbmonomer(dna = 3)
                >>
                self.DNA(binding = None,partType = part.pysbname + 'p3',upstream = 3,downstream = 1) +
                self.DNA(binding = 2,partType = part.pysbname + 'p4',upstream = 1) +
                self.RNAP(dna = 2,rna = None) +
                self.DNA(binding = 4,partType = part.pysbname + 'p2',downstream = 3) +
                regulatorMol.pysbmonomer(dna = 4),
                self.RNAPBindingLow)

        pysb.Rule('RNAP_Binding_To_' + part.pysbname+'_TF_on_p3',
                self.DNA(binding = 3,partType = part.pysbname + 'p3',upstream = 2,downstream = 1) +
                self.DNA(binding = None,partType = part.pysbname + 'p4',upstream = 1) +
                self.RNAP(dna = None,rna = None) +
                self.DNA(binding = None,partType = part.pysbname + 'p2',downstream = 2) +
                regulatorMol.pysbmonomer(dna = 3)
                >>
                self.DNA(binding = 4,partType = part.pysbname + 'p3',upstream = 3,downstream = 1) +
                self.DNA(binding = 2,partType = part.pysbname + 'p4',upstream = 1) +
                self.RNAP(dna = 2,rna = None) +
                self.DNA(binding = None,partType = part.pysbname + 'p2',downstream = 3) +
                regulatorMol.pysbmonomer(dna = 4),
                self.RNAPBindingLow)

        pysb.Rule('RNAP_Binding_To_' + part.pysbname+'_TF_on_p2_and_p3',
                self.DNA(binding = 3,partType = part.pysbname + 'p3',upstream = 2,downstream = 1) +
                self.DNA(binding = None,partType = part.pysbname + 'p4',upstream = 1) +
                self.RNAP(dna = None,rna = None) +
                self.DNA(binding = 4,partType = part.pysbname + 'p2',downstream = 2) +
                regulatorMol.pysbmonomer(dna = 3) +
                regulatorMol.pysbmonomer(dna = 4)
                >>
                self.DNA(binding = 4,partType = part.pysbname + 'p3',upstream = 3,downstream = 1) +
                self.DNA(binding = 2,partType = part.pysbname + 'p4',upstream = 1) +
                self.RNAP(dna = 2,rna = None) +
                self.DNA(binding = 5,partType = part.pysbname + 'p2',downstream = 3) +
                regulatorMol.pysbmonomer(dna = 4) +
                regulatorMol.pysbmonomer(dna = 5),
                self.RNAPBindingLow)

        # Transcription
        pysb.Rule('Transcription_Initiation_Of_'+part.pysbname,
                self.DNA(binding = 1,partType = part.pysbname + 'p4',downstream = 2) +
                self.RNAP(dna = 1,rna = None) +
                self.DNA(upstream = 2,binding = None)
                >>
                self.DNA(binding = None,partType = part.pysbname + 'p4',downstream = 3) +
                self.RNAP(dna = 1,rna = 2) +
                self.DNA(upstream = 3,binding = 1) +
                self.RNA(binding = None,upstream = None,downstream = 2,partType = part.pysbname),
                self.TranscriptionRate)



    def pysbPartRulesRBS(self,part):
        try:
            # Transcripton
            pysb.Rule('RBS_'+ part.pysbname + '_Transcription',
                          self.DNA(binding = 1,downstream = 2,partType = part.pysbname) +
                          self.RNAP(dna = 1,rna = 3) +
                          self.DNA(upstream = 2,binding = None) +
                          self.RNA(downstream = 3)
                          >>
                          self.DNA(binding = None,downstream = 2,partType = part.pysbname) +
                          self.RNAP(dna = 1,rna = 3) +
                          self.DNA(upstream = 2,binding = 1) +
                          self.RNA(downstream = 4) +
                          self.RNA(binding = None,downstream = 3,upstream = 4,partType = part.pysbname)
                          ,self.TranscriptionRate)

            # Translation
            pysb.Rule('RBS_' + part.pysbname + '_Ribosome_Binding',
                        self.RNA(binding = None,partType = part.pysbname) +
                        self.Ribosome(rna = None)
                        >>
                        self.RNA(binding = 1,partType = part.pysbname) +
                        self.Ribosome(rna = 1),
                        self.RBSBinding)
        except:
            pass
            # because we only need this rule once, but the method will be called per RBS







    def pysbPartRulesTerminator(self,part):
        try:
            pysb.Rule('Transcription_Termination_'+ part.pysbname,
                          self.DNA(binding = 1, partType = part.pysbname) +
                          self.RNAP(dna = 1,rna = 2) +
                          self.RNA(downstream = 2)
                          >>
                          self.DNA(binding = None, partType = part.pysbname) +
                          self.RNAP(dna = None,rna = None) +
                          self.RNA(downstream = None),
                          self.TranscriptionRate)
        except:
            pass
            # because we only need this rule once, but the method will be called per terminator



    def pysbPartRulesCodingRegion(self,part):
        codingfor = part.getAfterNodes(Molecule)[0]

        # Transcription

        pysb.Rule('Transcription_'+ part.pysbname,
                      self.DNA(binding = 1,downstream = 2,partType = part.pysbname) +
                      self.RNAP(dna = 1,rna = 3) +
                      self.DNA(upstream = 2,binding = None) +
                      self.RNA(downstream = 3)
                      >>
                      self.DNA(binding = None,downstream = 2,partType = part.pysbname) +
                      self.RNAP(dna = 1,rna = 3) +
                      self.DNA(upstream = 2,binding = 1) +
                      self.RNA(downstream = 4) +
                      self.RNA(binding = None, downstream = 3,upstream = 4,partType = part.pysbname),
                      self.TranscriptionRate)

        # Translation
        pysb.Rule('Translation_Initiation_'+part.pysbname,
                    self.RNA(binding = 2, downstream = 1) +
                    self.RNA(binding = None,upstream = 1,partType = part.pysbname) +
                    self.Ribosome(rna = 2)
                    >>
                    self.RNA(binding = None, downstream = 1) +
                    self.RNA(binding = 2,upstream = 1,partType = part.pysbname) +
                    self.Ribosome(rna = 2),
                    self.TranslationInitiationRate)

        pysb.Rule('Translation_'+part.pysbname,
                    self.RNA(binding = 1,partType = part.pysbname) +
                    self.Ribosome(rna = 1)
                    >>
                    self.RNA(binding = None,partType = part.pysbname) +
                    self.Ribosome(rna = None) +
                    codingfor.pysbmonomer(dna = None),
                    self.TranslationRate)

        # codingFor Degradation
        pysb.Rule('Degradation_of_'+str(codingfor),
                    codingfor.pysbmonomer(dna = None) >>
                    None,
                    self.TranscriptionFactorDegradation)


    def runKappaSimulation(self,weaverOutput):
        #generate populate dna / rna part lists
        for part in weaverOutput.partList:
            part.initForPysb()

        # create molecule monomers.
        for molecule in weaverOutput.moleculeList:
            molecule.initForPysb()

        pysb.Monomer('DNA',['binding','downstream','upstream','partType'],{'partType':self.dnaPartTypeList})
        pysb.Monomer('RNA',['binding','downstream','upstream','partType'],{'partType':self.rnaPartTypeList})

        circuitinitial = []

        # create part initials
        bindingIndex = 0
        for part in self.dnaPartTypeList[:-1]:
            if bindingIndex == 0:
                circuitinitial.append(pysb.MonomerPattern(self.DNA,{'partType' : part,'binding' : None,'upstream' : None,'downstream' : bindingIndex+1},None))
            else:
                circuitinitial.append(pysb.MonomerPattern(self.DNA,{'partType' : part,'binding' : None,'upstream' : bindingIndex,'downstream' : bindingIndex+1},None))
            bindingIndex += 1

        circuitinitial.append(pysb.MonomerPattern(self.DNA,{'partType' : self.dnaPartTypeList[-1],'binding' : None,'upstream' : bindingIndex,'downstream' : None},None))
        pattern = pysb.ComplexPattern(circuitinitial,None)
        pysb.Initial(pysb.ComplexPattern(circuitinitial,None),self.OperonCount)

        #generate rules
        for part in weaverOutput.partList:
            part.pysbPartRules()

        # additional rules
        pysb.Rule('RNAP_falloff',
                self.DNA(binding = 1,downstream = 3) +
                self.RNAP(dna = 1,rna = 2) +
                self.RNA(downstream = 2) +
                self.DNA(upstream = 3,binding=pysb.ANY)
                >>
                self.DNA(binding = None,downstream = 1) +
                self.RNAP(dna = None,rna = None)+
                self.RNA(downstream = None) +
                self.DNA(upstream = 1,binding=pysb.ANY),
                self.RNAPfalloffrate)

        pysb.Rule('Ribosome_falloff',
                self.Ribosome(rna = 1) %
                self.RNA(binding = 1)
                >>
                self.Ribosome(rna = None) +
                self.RNA(binding = None),
                self.RibosomeFalloff)

        pysb.Rule('RNA_degradation',
                self.RNA(binding = None,downstream = None)
                >>
                None,
                self.RNAdegradation)


        for obs in weaverOutput.moleculeList:
            pysb.Observable('obs'+str(obs),obs.pysbmonomer())


        x = kappa.run_simulation(self.model, time=25000)
        legend = []


        for obs in weaverOutput.moleculeList:
            plt.plot(x['time'],x['obs'+str(obs)])
            legend.append(str(obs))

        plt.legend(legend,loc='upper left')
        plt.show()
