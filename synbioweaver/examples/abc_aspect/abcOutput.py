from synbioweaver import *
import numpy
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import stats

class Reaction:
    param_counter = 0

    def __init__(self, reactants, products):
        react = ' + '.join(map(str, reactants))
        prod = ' + '.join(map(str, products))
        self.reaction = react + ' -> ' + prod

    def rates(self, reactants):
        Reaction.param_counter += 1
        rt = '*'.join(map(str, reactants))
        self.param = 'p'+str(Reaction.param_counter)
        self.rate = rt+'*'+'p'+str((Reaction.param_counter))
        return [self.rate, self.param]


class DistanceFunction:
    def __init__(self, distance_file):
        self.distance_file = distance_file
        self.number_of_distances = 3
        self.targets_in_order = [0, 20, 0]
        self.times_for_targets = [[0, 20], [21, 70], [71, 100]]

    def writeDistanceFile(self):
        distance_name=[]
        main_file = open(self.distance_file, 'w')
        main_file.write('import numpy as np' + '\n')
        main_file.write('\n')
        main_file.write('\n')
        main_file.write('def distance(data1, data2, parameters, model):' + '\n')
        main_file.write('\n')
        main_file.write('    #data1 is simulated' + '\n')
        main_file.write('    #data2 is real' + '\n')
        main_file.write('\n')
        segment_counter=-1
        for g in range(self.number_of_distances):
            segment_counter += 1
            distance_name.append('d' + str(segment_counter))
            main_file.write('    d' + str(segment_counter) + ' = 0' + '\n')
            main_file.write('    for i in range(' + str(self.times_for_targets[segment_counter][0]) + ', ' + str(self.times_for_targets[segment_counter][1]) + '):' + '\n')
            main_file.write('        d' + str(segment_counter) + ' += (data1[i,0] - ' + str(self.targets_in_order[segment_counter]) + ')*(' + 'data1[i,0] - ' + str(self.targets_in_order[segment_counter]) + ')' + '\n')
            main_file.write('    d' + str(segment_counter) + ' = np.sqrt('+'d'+ str(segment_counter) + '/(' + str(self.times_for_targets[segment_counter][1]) + ' - ' + str(self.times_for_targets[segment_counter][0]) + '))'+ '\n')
        temp_var = ','.join(distance_name)
        main_file.write('\n')
        main_file.write('    return [' + temp_var + ']')
        main_file.close()
        return self.distance_file


class InputFile:
    def __init__(self, input_file_name, cuda_file_name, param_list, molecule_list, epsilon, fit, particles, alpha, prior_distribution, init_cond__distribution, priors, initial_conditions):
        self.input_file_name = input_file_name
        self.cuda_file_name = cuda_file_name
        self.param_list = param_list
        self.molecule_list = molecule_list
        self.epsilon = epsilon
        self.fit = fit
        self.particles = particles
        self.alpha = alpha
        self.simtype = 'Gillespie'
        self.prior_distribution = prior_distribution
        self.init_cond__distribution = init_cond__distribution
        self.priors = priors
        self.initial_conditions = initial_conditions

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
        out_file.write('<times>   0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100  </times>' + '\n')
        out_file.write('# variables: For ABC SMC, whitespace delimited lists of concentrations (ODE or SDE) or molecule numbers (Gillespie)' + '\n')
        out_file.write('# Denote your data via tags <v1> </v1> or <var1> </var1> or <v2> </v2> etc. The tags are ignored and the data read in order' + '\n')
        out_file.write('# For simulation these data are ignored'+ '\n')
        out_file.write('# See fitting instruction below if the dimensionality of your data sets differ from the dimensionality of your model' + '\n')
        out_file.write('<variables>' + '\n')
        out_file.write(' <var1> </var1>' + '\n')
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
        out_file.write('# Otherwise, give a whitespace delimited list of fitting instrictions the same length as the dimensions of your data.' + '\n')
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
                out_file.write('<' + str(species) + '> ' + self.init_cond__distribution[init_cond_count] + ' ' + str(self.initial_conditions[init_cond_count][0]) + ' ' + str(self.initial_conditions[init_cond_count][1]) + ' </' + str(species) + '>' + '\n')
            elif len(str(self.initial_conditions[init_cond_count]).split()) == 1:
                out_file.write('<' + str(species) + '> ' + self.init_cond__distribution[init_cond_count] + ' ' + str(self.initial_conditions[init_cond_count]) + ' </' + str(species) + '>' + '\n')
        out_file.write('</initial>' + '\n')
        out_file.write('\n')
        out_file.write('<parameters>' + '\n')
        prior_count = -1
        out_file.write('<p0> constant 1 </p0>' + '\n')

        for param in self.param_list:
            prior_count += 1
            out_file.write('<' + str(param) + '> ' + self.prior_distribution[prior_count] + ' ' + str(self.priors[prior_count]) + ' </' + str(param) + '>' + '\n')
        out_file.write('</parameters>' + '\n')
        out_file.write('\n')
        out_file.write('</model1>' + '\n')
        out_file.write('</models>' + '\n')
        out_file.write('</input>' + '\n')
        out_file.close()
        return self.input_file_name


class CudaOutput(Aspect):

    def mainAspect(self):
        self.addTypeAdvice(PartSignature('*.PositivePromoter+'), self.posPromoterCuda, 'cudacode')
        self.addTypeAdvice(PartSignature('*.NegativePromoter+'), self.negPromoterCuda, 'cudacode')
        self.addTypeAdvice(PartSignature('*.CodingRegion+'), self.reactionsCuda, 'cudacode')
        self.addWeaverOutput(self.cudaOutput)

    def negPromoterCuda(self, part):
        partname = part.__class__.__name__
        regulator = part.getRegulatedBy()[0]
        complx = str(partname) + '_' + str(regulator)
        reaction = Reaction([partname, regulator], [complx])
        rate, param = reaction.rates([partname, regulator])
        return [reaction.reaction, rate, param, partname]

    def posPromoterCuda(self, part):
        #TO DO
        pass

    def reactionsCuda(self,part):
        partname = part.getCodingFor()[0]
        before_part = part.getBeforePart()
        before_name = before_part.__class__.__name__
        reaction = Reaction([before_name], [partname, before_name])
        degradation = Reaction([partname], [])
        rated, paramd = degradation.rates([partname])
        rate, param = reaction.rates([before_name])
        return [[reaction.reaction, degradation.reaction], [rate, rated], [param, paramd], partname]

    def stoichiometry(self, moleculeList, reaction_list):
        #Build the stoichiometry matrix
        self.stoichiometry_matrix = numpy.zeros(shape=(len(reaction_list), len(moleculeList)))
        molecule_counter = 0
        #For each molecule, go through each reaction
        for molecule in moleculeList:
            if '(' in str(molecule):
                molecule, regulator = str(molecule).split('(')
            reaction_counter = -1
            for rction in reaction_list:
                reaction_counter += 1
                if str(molecule) in rction:
                    #If the species was found in the reaction, split it to reactants and products. Then count number of times you found it (remove the ones that are in complex)
                    #If in reactants (-) if in products (+)
                    splitreactprod = rction.split('->')
                    counterreac = splitreactprod[0].count(str(molecule)) - splitreactprod[0].count('_'+str(molecule))
                    if counterreac > 0:
                        self.stoichiometry_matrix[reaction_counter, molecule_counter] = self.stoichiometry_matrix[reaction_counter, molecule_counter] - counterreac
                    counterprod = splitreactprod[1].count(str(molecule)) + splitreactprod[1].count('_'+str(molecule)) - splitreactprod[1].count(str(molecule) + '_')
                    if counterprod > 0:
                        self.stoichiometry_matrix[reaction_counter, molecule_counter] = self.stoichiometry_matrix[reaction_counter, molecule_counter] + counterprod
                if reaction_counter == len(reaction_list)-1:
                    molecule_counter += 1
                    break
        return self.stoichiometry_matrix

    def writeCuda(self, stoichiometry_matrix, reaction_list, molecule_list, rates_list, param_list):

        self.cuda_file = str('CudaOutput.cu')
        self.name = 'CudaOutput'
        out_file = open(self.cuda_file,'w')
        out_file.write("#define NSPECIES " + str(len(molecule_list)) + "\n")
        out_file.write("#define NPARAM " + str(len(rates_list)+1) + "\n")
        out_file.write("#define NREACT " + str(len(reaction_list)) + "\n")
        out_file.write("\n")
        out_file.write('#define leq(a,b) a<=b' + '\n')
        out_file.write('#define neq(a,b) a!=b' + '\n')
        out_file.write('#define geq(a,b) a>=b' + '\n')
        out_file.write('#define lt(a,b) a<b' + '\n')
        out_file.write('#define gt(a,b) a>b' + '\n')
        out_file.write('#define eq(a,b) a==b' + '\n')
        out_file.write('#define and_(a,b) a&&b' + '\n')
        out_file.write('#define or_(a,b) a||b' + '\n')
        out_file.write("\n")
        for i in range(len(molecule_list)):
            if '(' in str(molecule_list[i]):
                name, regulator = str(molecule_list[i]).split('(')
                out_file.write("#define " + str(name) + " y[" + str(i) + "]" + "\n")
            else:
                out_file.write("#define " + str(molecule_list[i]) + " y[" + str(i) + "]" + "\n")
        out_file.write("\n")
        out_file.write("#define p0 tex2D(param_tex,0,tid)" + "\n")
        for i in range(len(param_list)):
            out_file.write("#define " + str(param_list[i]) + " tex2D(param_tex," + str(i+1) + ",tid)" + "\n")
        out_file.write("\n")
        out_file.write("__constant__ int smatrix[]={"+"\n")

        numpy.savetxt(out_file, stoichiometry_matrix, fmt='%s,', )

        out_file.write("};"+"\n"+"\n"+"\n")

        out_file.write("__device__ void hazards(int *y, float *h, float t, int tid){ " + "\n")
        eq_count = -1
        for i in rates_list:
            if i is not None:
                eq_count += 1
                out_file.write("\t" + "h[" + str(eq_count) + "] = p0*" + str(i) + ';' +"\n")
        out_file.write("}")
        out_file.close()
        return self.name

    def equations(self, stoichiometry_matrix, rates_list):
        self.equns = [None]*len(rates_list)
        self.lenSRow = len(stoichiometry_matrix[0])
        self.lenS = len(stoichiometry_matrix)
        self.lenh = len(rates_list)
        for rowS in range(self.lenS):
            non_zero_indeces = []
            for elementRow in range(self.lenSRow):
                if stoichiometry_matrix[rowS][elementRow] != 0:
                    non_zero_indeces.append(elementRow)
            for non_zero in range(len(non_zero_indeces)):
                if non_zero == 0:
                    self.equns[rowS] = str(stoichiometry_matrix[rowS][non_zero_indeces[non_zero]]) + ' * ' + str(rates_list[non_zero_indeces[non_zero]])
                else:
                    if stoichiometry_matrix[rowS][non_zero_indeces[non_zero]] >=0:
                        self.equns[rowS] += ' + '+str(stoichiometry_matrix[rowS][non_zero_indeces[non_zero]]) + ' * ' + str(rates_list[non_zero_indeces[non_zero]])
                    elif stoichiometry_matrix[rowS][non_zero_indeces[non_zero]] < 0:
                        self.equns[rowS] += str(stoichiometry_matrix[rowS][non_zero_indeces[non_zero]]) + ' * ' + str(rates_list[non_zero_indeces[non_zero]])
                if non_zero == len(non_zero_indeces):
                    self.equns[rowS] += '\n'
        return self.equns

    def runCudaSim(self):
        runfile=open('out.sh', 'w')
        runfile.write('export PYTHONPATH=$PYTHONPATH:/home/ucbtle1/abc-sysbio-code/:/home/ucbtle1/cuda-sim-code/' + '\n')
        runfile.write('exe=/home/ucbtle1/abc-sysbio-code/scripts/run-abc-sysbio'+ '\n')
        runfile.write('\n')
        runfile.write('python -u $exe --localcode --infile input_file.xml -f -of=results_simulation --cuda  --custd -S')
        runfile.close()
        os.system('./out.sh')
        return runfile

    def plot_simulation(self):
        l = np.loadtxt('results_simulation/trajectories.txt')
        self.b = scipy.delete(l, 0, 1)
        self.b = scipy.delete(self.b, 0, 1)
        self.b = scipy.delete(self.b, 0, 1)
        self.b = scipy.delete(self.b, 0, 1)
        nt = len(self.b)
        self.iA = []
        self.iB = []
        self.vA = []
        self.vB = []
        for i in range(0, nt, 2):
            self.iA.append(i)
        for i in range(1, nt, 2):
            self.iB.append(i)
        for i in self.iA:
            self.vA.append(self.b[i][-1])
        for i in self.iB:
            self.vB.append(self.b[i][-1])

        self.xmax = max(self.vA)
        self.ymax = max(self.vB)
        x = np.linspace(0, self.xmax, 1000)
        y = np.linspace(0, self.xmax, 1000)
        X, Y = np.meshgrid(x, y)
        self.positions = np.vstack([X.ravel(), Y.ravel()])
        self.values = np.vstack([self.vB, self.vA])
        self.kernel = stats.gaussian_kde(self.values)
        self.Z = np.reshape(self.kernel(self.positions).T, X.shape)
        cmap = plt.cm.autumn
        #norm = plt.cm.colors.Normalize(vmax=abs(Z).max(), vmin=-abs(Z).max())
        extent = [0, self.xmax, 0, self.ymax]
        self.fig = plt.figure()
        ax = self.fig.add_subplot(111)
        #Uncomment line below if you want filled contours
        #ax.imshow(np.rot90(Z), extent=extent, cmap=cmap, norm=norm)
        v = plt.axis()
        ax.contour(np.rot90(self.Z), hold='on', cmap=cmap, origin='upper', extent=extent, linewidth=.5)
        plt.axis(v)
        ax.set_xlim([0, self.xmax])
        ax.set_ylim([0, self.ymax])
        plt.xlabel('N(AA)')
        plt.ylabel('N(BB)')
        self.title = 'results_simulation/simulation_plot.pdf'
        plt.savefig(self.title)
        return self.title

    def cudaOutput(self, weaverOutput):
        reaction_list = []
        rates_list = []
        param_list = []
        self.initial_conditions = []
        self.init_cond__distribution = []
        self.prior_distribution = []
        self.priors = []
        self.weaverOutput = weaverOutput

        """
        For each part, go to 'mainAspect' and check for what kind of part it is. Then for each one get the reaction and the rate and then calculate the stoichiometry matrix
        If its a negative promoter, get the regulator and make the repression reaction. Regulator + Promoter -> Regulator_promoter complex. *!*manually set the _ in complex*!*
        If its a coding region, get the promoter (before part) and the name of what the coding region codes for. Promoter -> CodesFor + Promoter
        """

        for part in weaverOutput.partList:
            reac, rat, params = part.cudacode()[0:3]
            if type(params) == type(list()):
                for i in range(0, 2):
                    reaction_list.append(reac[i])
                    rates_list.append(rat[i])
                    param_list.append(params[i])
            else:
                reaction_list.append(reac)
                rates_list.append(rat)
                param_list.append(params)

        for mol in weaverOutput.moleculeList:
            tmp = raw_input("Enter initial condition value for " + str(mol) + " (distribution start end): ")
            tmpinp = tmp.split(" ")
            if len(tmpinp) > 2:
                self.init_cond__distribution.append(tmpinp[0])
                self.initial_conditions.append([tmpinp[1], tmpinp[2]])
            elif len(tmpinp) == 2:
                self.init_cond__distribution.append(tmpinp[0])
                self.initial_conditions.append(tmpinp[1])

        for reac in reaction_list:
            tmp = raw_input("Enter value of parameter for " + str(reac) + ": ")
            tmpinp = tmp.split(" ")
            if len(tmpinp) > 2:
                self.prior_distribution.append(tmpinp[0])
                self.priors.append([tmpinp[1], tmpinp[2]])
            elif len(tmpinp) == 2:
                self.prior_distribution.append(tmpinp[0])
                self.priors.append(tmpinp[1])

        self.epsilon = raw_input("Enter value for epsilon: ")
        fit_tmp = str(raw_input("Enter species to fit: "))
        fit_tmpsp = fit_tmp.split(" ")

        mol_strings = []

        for i in weaverOutput.moleculeList:
            mol_strings.append(str(i))

        if len(fit_tmpsp) == 1:
            self.fit = 'species' + str(mol_strings.index(str(fit_tmp)))
        elif len(fit_tmpsp) > 1:
            self.fit = ""
            for i in range(len(fit_tmpsp)):
                self.fit += 'species' + str(mol_strings.index(str(fit_tmpsp[i]))) + " "
        self.particles = raw_input("Enter value for particles: ")
        self.alpha = raw_input("Enter value for alpha: ")

        #Dimerisation: if a molecule is not found in the products list, then you have to create a reaction for it.
        products_list = []
        for reaction in reaction_list:
            splitreactprod = reaction.split('->')
            products_list_a = splitreactprod[1].split('+')
            for a in products_list_a:
                a = str(a.replace(" ", ""))
                products_list.append(a)
        for mol in weaverOutput.moleculeList:
            if '(' in str(mol):
                mol, regulator = str(mol).split('(')
            if str(mol) not in products_list and '_' in str(mol):
                reactants = str(mol).split('_')
                extra_reaction = Reaction(reactants, [mol])
                reaction_list.append(extra_reaction.reaction)
                rate, param = extra_reaction.rates(reactants)
                rates_list.append(rate)
                param_list.append(param)

        stoichiometry_matrix = self.stoichiometry(weaverOutput.moleculeList, reaction_list)
        cuda_file = self.writeCuda(stoichiometry_matrix, reaction_list, weaverOutput.moleculeList, rates_list, param_list)
        switch_input_file = InputFile('input_file.xml', str(cuda_file), param_list, weaverOutput.moleculeList, self.epsilon, self.fit, self.particles, self.alpha, self.prior_distribution, self.init_cond__distribution, self.priors, self.initial_conditions)
        switch_input_file.writeInput()
        switch_distance_file = DistanceFunction('customABC.py')
        switch_distance_file.writeDistanceFile()
        simulation = self.runCudaSim()
        plot_sim = self.plot_simulation()
        return cuda_file, switch_input_file.input_file_name, switch_distance_file.distance_file


