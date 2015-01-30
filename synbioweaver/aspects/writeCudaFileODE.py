from synbioweaver.core import *
import numpy


class WriteCudaFileODE(Aspect):

    def mainAspect(self):
        self.addWeaverOutput(self.writeCudaFileODE)

    def writeCuda(self, reaction_list, molecule_list, params, nspecies, nreactions):
        self.cuda_file = str('model.cu')
        self.name = 'model'
        out_file = open(self.cuda_file,'w')
        out_file.write("#define NSPECIES " + str(nspecies) + "\n")
        out_file.write("#define NPARAM " + str(nreactions+1) + "\n")
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
        for i in range(len(params)):
            out_file.write("#define " + str(params[i]) + " tex2D(param_tex," + str(i+1) + ",tid)" + "\n")
        out_file.write("\n")
        out_file.write("__device__ void operator()(int *neq, double *t, double *y, double *ydot/*, void *otherData*/){ " + "\n")

        for i in range(len(self.odes)):
            out_file.write('ydot[' + str(i) + ']=(' + self.odes[i] + ')' + '\n')
        out_file.write("};")
        out_file.write("struct myJex{" + "\n")
        out_file.write("__device__ void operator()(int *neq, double *t, double *y, int ml, int mu, double *pd, int nrowpd/*, void *otherData*/){"+"\n")
        out_file.write("return; "+"\n"+"}"+"\n"+"};")
        out_file.close()
        return self.name

    def calculateODEs(self):
        self.odes = []
        for sp in range(len(self.species)):
            ode = ''
            for rt in range(len(self.rates)):
                if self.stoichiometry_matrix.T[sp][rt] > 0:
                    ode +='+'
                    ode += str(self.stoichiometry_matrix.T[sp][rt])
                    ode +=' * '
                    ode += str(self.rates[rt])
                elif self.stoichiometry_matrix.T[sp][rt] < 0:
                    ode +=' '
                    ode += str(self.stoichiometry_matrix.T[sp][rt])
                    ode +='*'
                    ode += str(self.rates[rt])
                elif self.stoichiometry_matrix.T[sp][rt] == 0:
                    pass
            self.odes.append(ode)
        return self.odes

    def writeCudaFileODE(self, weaverOutput):
        if getattr(weaverOutput, "getContext", None) != None:
            self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getContext()
        else:
            if getattr(weaverOutput, "getReactions", None) != None:
                self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getReactions()
        self.rates = []
        for i in range(self.nreactions):
            self.rates.append(self.reactions[i].rate)
        self.calculateODEs()
        self.writeCuda(self.reactions, self.species, self.parameters, self.nspecies, self.nreactions)


