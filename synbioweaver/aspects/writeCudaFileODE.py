from synbioweaver.core import *
from synbioweaver.aspects.modelDefinitions import *
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
        out_file.write('#define zero 1' + '\n') # hack
        
        for i in range(len(molecule_list)):
            if '(' in str(molecule_list[i]):
                name, regulator = str(molecule_list[i]).split('(')
                out_file.write("#define " + str(name) + " y[" + str(i) + "]" + "\n")
            else:
                out_file.write("#define " + str(molecule_list[i]) + " y[" + str(i) + "]" + "\n")
        out_file.write("\n")
        out_file.write("#define p0 tex2D(param_tex,0,tid)" + "\n")

        # write parameters out in order of reactions
        parcount = 0
        for i in range(len(reaction_list)):
            for j in range(len(reaction_list[i].param)):
                out_file.write("#define " + str(reaction_list[i].param[j]) + " tex2D(param_tex," + str(parcount+1) + ",tid)" + "\n")
                parcount += 1
                
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
        for sp in range(self.model.nspecies):
            ode = ''
            for rt in range(len(self.rates)):
                if self.model.stoichiometry_matrix.T[sp][rt] > 0:
                    ode +='+'
                    ode += str(self.model.stoichiometry_matrix.T[sp][rt])
                    ode +=' * '
                    ode += str(self.rates[rt])
                elif self.model.stoichiometry_matrix.T[sp][rt] < 0:
                    ode +=' '
                    ode += str(self.model.stoichiometry_matrix.T[sp][rt])
                    ode +='*'
                    ode += str(self.rates[rt])
                elif self.model.stoichiometry_matrix.T[sp][rt] == 0:
                    pass
            self.odes.append(ode)
        return self.odes

    def writeCudaFileODE(self, weaverOutput):
       
        # We are expecting either a set of reactions or a context

        if getattr(weaverOutput, "getContext", None) != None:
            self.model = weaverOutput.getContext()
        else:
            if getattr(weaverOutput, "getReactions", None) != None:
                self.model = weaverOutput.getReactions()
            else:
                print " writeCudaFileODE : Neither getContext() or getReactions() is available. Quitting"
                exit()

        self.rates = []
        for i in range(self.model.nreactions):
            self.rates.append(self.model.reactions[i].rate)
        self.calculateODEs()

        #print "writeCudaFileODE:", self.model.parameters

        return self.writeCuda(self.model.reactions, self.model.species, self.model.parameters, self.model.nspecies, self.model.nreactions)


