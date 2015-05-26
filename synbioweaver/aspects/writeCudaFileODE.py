from synbioweaver.core import *
import numpy, StringIO


class WriteCudaFileODE(Aspect):

    def mainAspect(self):
        self.addWeaverOutput(self.writeCudaFileODE)

    def writeCuda(self, reaction_list, molecule_list, params, nspecies, nreactions):
        out_file = StringIO.StringIO()

        out_file.write("#define NSPECIES " + str(nspecies) + "\n")
        out_file.write("#define NPARAM " + str(len(params)+1) + "\n")
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
        out_file.write("struct myFex{" + '\n')
        out_file.write("__device__ void operator()(int *neq, double *t, double *y, double *ydot/*, void *otherData*/){ " + "\n")
        out_file.write("int tid = blockDim.x * blockIdx.x + threadIdx.x;" + "\n")
        for i in range(len(self.odes)):
            out_file.write('ydot[' + str(i) + "] = p0*(" + self.odes[i] + ');' + '\n')
        out_file.write("}" + "\n")
        out_file.write("};" + "\n")
        out_file.write("struct myJex{" + "\n")
        out_file.write("__device__ void operator()(int *neq, double *t, double *y, int ml, int mu, double *pd, int nrowpd/*, void *otherData*/){"+"\n")
        out_file.write("return; "+"\n"+"}"+"\n"+"};"+"\n")
        
        return out_file.getvalue()

    def calculateODEs(self):
        self.odes = []
        for sp in range(self.nspecies):
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
       
        # We are expecting either a set of reactions or a context

        if getattr(weaverOutput, "getContext", None) != None:
            self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getContext()
            print "HERE:", self.species, self.nspecies
        else:
            if getattr(weaverOutput, "getReactions", None) != None:
                self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getReactions()
            else:
                print " writeCudaFileODE : Neither getContext() or getReactions() is available. Quitting"
                exit()

        self.rates = []
        for i in range(self.nreactions):
            self.rates.append(self.reactions[i].rate)
        self.calculateODEs()

        return self.writeCuda(self.reactions, self.species, self.parameters, self.nspecies, self.nreactions)


