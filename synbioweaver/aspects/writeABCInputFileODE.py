from synbioweaver.core import *
from synbioweaver.aspects.printReactionNetworkAspect import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy as np
import readData
class InputFile:
    def __init__(self, input_file_name, params, molecule_list, epsilon, fit, particles, alpha, prior_distribution, init_cond_distribution, priors, initial_conditions, data):
        self.input_file_name = input_file_name
        self.cuda_file_name = 'model'
        self.params = params
        self.molecule_list = molecule_list
        self.epsilon = epsilon
        self.fit = fit
        self.particles = particles
        self.alpha = alpha
        self.simtype = 'ODE'
        self.prior_distribution = prior_distribution
        self.init_cond_distribution = init_cond_distribution
        self.priors = priors
        self.initial_conditions = initial_conditions
        self.data = data

    def writeInput(self):
        out_file = open(self.input_file_name, 'w')
        out_file.write('<input>' + '\n')
        out_file.write('\n')
        out_file.write('# Number of models for which details are described in this input file' + '\n')
        out_file.write('<modelnumber> 1 </modelnumber>' + '\n')
        out_file.write('# Restart from previous (pickled) population?' + '\n')
        out_file.write('<restart> False </restart>' + '\n')
        out_file.write('# Automatic epsilon schedule. Provide a vector of final epsilons and the alpha (defaults to 0.9)' + '\n')
        out_file.write('<autoepsilon>' + '\n')
        out_file.write('<finalepsilon> ' + str(self.epsilon) + '</finalepsilon>' + '\n')
        out_file.write('<alpha> ' + str(self.alpha) + '</alpha>' + '\n')
        out_file.write('</autoepsilon>' + '\n')
        out_file.write('<particles> ' + str(self.particles) + '</particles>' + '\n')
        out_file.write('# Beta is the number of times to simulate each sampled parameter set.' + '\n')
        out_file.write('# This is only applicable for models simulated using Gillespie and SDE' + '\n')
        out_file.write('<beta> 1 </beta>' + '\n')
        out_file.write('# Internal timestep for solver.' + '\n')
        out_file.write('# Make this small for a stiff model.' + '\n')
        out_file.write('<dt> 0.01 </dt>' + '\n')
        out_file.write('# The pertubation kernels are computed with respect to the previous parameter distribution'+ '\n')
        out_file.write('# Currently uniform and normal are implemented' + '\n')
        out_file.write('<kernel> uniform </kernel>' + '\n')
        out_file.write('# Probability of perturbing the sampled model (ignored when modelnumber = 1)' + '\n')
        out_file.write('<modelkernel> 0.7 </modelkernel>' + '\n')
        out_file.write('# rtol and atol can be specified here. If the model is stiff then setting these to small might help the simulation to run' + '\n')
        out_file.write('#<rtol> </rtol>' + '\n')
        out_file.write('#<atol> </atol>' + '\n')
        out_file.write('\n')
        out_file.write('<data>' + '\n')
        out_file.write('# times: For ABC SMC, times must be a whitespace delimited list'+ '\n')
        out_file.write('# In simulation mode these are the timepoints for which the simulations will be output' + '\n')
        out_file.write('<times> ')
        for i in self.data:
             out_file.write(str(i[0]) + ' ')
        out_file.write('</times>' + '\n')
        out_file.write('# variables: For ABC SMC, whitespace delimited lists of concentrations (ODE or SDE) or molecule numbers (Gillespie)' + '\n')
        out_file.write('# Denote your data via tags <v1> </v1> or <var1> </var1> or <v2> </v2> etc. The tags are ignored and the data read in order' + '\n')
        out_file.write('# For simulation these data are ignored'+ '\n')
        out_file.write('# See fitting instruction below if the dimensionality of your data sets differ from the dimensionality of your model' + '\n')
        out_file.write('<variables>' + '\n')
        for i in range(1, self.data.shape[1], 2):
            out_file.write(' <var' + str(i)+'> ')
            for k in self.data:
                out_file.write(str(k[i]) + ' ')
            out_file.write('</var' + str(i)+'> ' + '\n')
        out_file.write('</variables>' + '\n')
        out_file.write('</data>' + '\n')
        out_file.write('<models>' + '\n')
        out_file.write('<model1>' + '\n')
        out_file.write('\n')
        out_file.write('<name> ' + self.cuda_file_name + ' </name>' + '\n')
        out_file.write('<source> ' + self.cuda_file_name + ' </source>' + '\n')
        out_file.write('# type: the method used to simulate your model. ODE, SDE or Gillespie.' + '\n')
        out_file.write('<type> ' + self.simtype + ' </type>' + '\n')
        out_file.write('# Fitting information. If fit is None, all species in the model are fitted to the data in the order they are listed in the model.' + '\n')
        out_file.write('# Otherwise, give a whitespace delimited list of fitting instructions the same length as the dimensions of your data.' + '\n')
        out_file.write('# Use speciesN to denote the Nth species in your model. Simple arithmetic operations can be performed on the species from your model.' + '\n')
        out_file.write('# For example, to fit the sum of the first two species in your model to your first variable, write fit: species1+species2' + '\n')
        out_file.write('<fit> ' + self.fit + ' </fit>' + '\n')
        out_file.write('# Priors on initial conditions and parameters:' + '\n')
        out_file.write('# one of ' + '\n')
        out_file.write('#       constant, value ' + '\n')
        out_file.write('#       normal, mean, variance ' + '\n')
        out_file.write('#       uniform, lower, upper' + '\n')
        out_file.write('#       lognormal, mean, variance ' + '\n')
        out_file.write('\n')
        out_file.write('<initial>' + '\n')
        init_cond_count = -1
        for species in self.molecule_list:
            init_cond_count += 1
            if '(' in str(species):
                species, regulator = str(species).split('(')
            if len(str(self.initial_conditions[init_cond_count]).split()) > 1:
                out_file.write('<' + str(species) + '> ' + self.init_cond_distribution[init_cond_count] + ' ' + str(self.initial_conditions[init_cond_count][0]) + ' ' + str(self.initial_conditions[init_cond_count][1]) + ' </' + str(species) + '>' + '\n')
            elif len(str(self.initial_conditions[init_cond_count]).split()) == 1:
                out_file.write('<' + str(species) + '> ' + self.init_cond_distribution[init_cond_count] + ' ' + str(self.initial_conditions[init_cond_count]) + ' </' + str(species) + '>' + '\n')
        out_file.write('</initial>' + '\n')
        out_file.write('\n')
        out_file.write('<parameters>' + '\n')
        prior_count = -1
        out_file.write('<p0> constant 1 </p0>' + '\n')

        for param in self.params:
            prior_count += 1
            if len(str(self.priors[prior_count]).split()) > 1:
                out_file.write('<' + str(param) + '> ' + self.prior_distribution[prior_count] + ' ' + str(self.priors[prior_count][0]) + ' ' + str(self.priors[prior_count][1])  + ' </' + str(param) + '>' + '\n')
            elif len(str(self.priors[prior_count]).split()) == 1:
                out_file.write('<' + str(param) + '> ' + self.prior_distribution[prior_count] + ' ' + str(self.priors[prior_count]) + ' </' + str(param) + '>' + '\n')
        out_file.write('</parameters>' + '\n')
        out_file.write('\n')
        out_file.write('</model1>' + '\n')
        out_file.write('</models>' + '\n')
        out_file.write('</input>' + '\n')
        out_file.close()
        return self.input_file_name


class WriteABCInputFileODE(Aspect):

    def mainAspect(self):
        self.addWeaverOutput(self.writeABCInputFileODE)

    def writeABCInputFileODE(self, weaverOutput):
        if getattr(weaverOutput, "getContext", None) != None:
            self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getContext()
        else:
            if getattr(weaverOutput, "getReactions", None) != None:
                self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getReactions()
        self.reaction_list = []
        self.initial_conditions = []
        self.init_cond_distribution = []
        self.prior_distribution = []
        self.priors = []
        self.rates = []
        for i in range(self.nreactions):
            self.rates.append(self.reactions[i].rate)
            self.reaction_list.append(self.reactions[i])

        for mol in weaverOutput.moleculeList:
            self.init_cond_distribution.append('uniform')
            tmp2 = '0 1'
            tmpinp2 = tmp2.split(" ")
            self.initial_conditions.append([tmpinp2[0], tmpinp2[1]])

        for i in range(len(self.reaction_list)):

            if self.reactions[i].process == "dnaBind":
                tmp = '1 2'
            elif self.reactions[i].process == "dnaUnbind":
                tmp = '2 3'
            elif self.reactions[i].process == "rnaTransc":
                tmp = '4 5'
            elif self.reactions[i].process == "rnaDeg":
                tmp = '6 7'
            elif self.reactions[i].process == "proteinTransl":
                tmp = '8 9'
            elif self.reactions[i].process == "proteinExp":
                tmp = '10 11'
            elif self.reactions[i].process == "proteinDeg":
                tmp = '12 13'
            elif self.reactions[i].process == "complexAss":
                tmp = '14 15'
            elif self.reactions[i].process == "complexDiss":
                tmp = '16 17'
            elif self.reactions[i].process == "complexDeg":
                tmp = '18 19'
            elif self.reactions[i].process == "context":
                tmp = '20 21'

            self.prior_distribution.append('uniform')
            tmpinp = tmp.split(" ")
            self.priors.append([tmpinp[0], tmpinp[1]])

        self.epsilon = 1.0
        fit_tmp = 'GFP'
        fit_tmpsp = fit_tmp.split(" ")
        self.particles = 10
        self.alpha = 0.1
        mol_strings = []

        for i in weaverOutput.moleculeList:
            mol_strings.append(str(i))

        if len(fit_tmpsp) == 1:
            self.fit = 'species' + str(mol_strings.index(str(fit_tmp)))
        elif len(fit_tmpsp) > 1:
            self.fit = ""
            for i in range(len(fit_tmpsp)):
                self.fit += 'species' + str(mol_strings.index(str(fit_tmpsp[i]))) + " "

        dat = readData.DataMatrix('data.txt')
        self.data = dat.readData()
        switch_input_file = InputFile('input_file.xml', self.parameters, weaverOutput.moleculeList, self.epsilon, self.fit, self.particles, self.alpha, self.prior_distribution, self.init_cond_distribution, self.priors, self.initial_conditions, self.data)
        switch_input_file.writeInput()
