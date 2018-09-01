from synbioweaver.core import *
from synbioweaver.aspects.modelDefinitions import *
import numpy, StringIO

class WriteMultipleCudaFileODE(Aspect):

    def mainAspect(self):
        self.addWeaverOutput(self.writeMultipleCudaFileODE)

    def writeCuda(self, reaction_list, molecule_list, params, nspecies, nreactions, odes):
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
        out_file.write('#define zero 1.0' + '\n') # hack
        
        for i in range(len(molecule_list)):
            if '(' in str(molecule_list[i]):
                name, regulator = str(molecule_list[i]).split('(')
                out_file.write("#define " + str(name) + " y[" + str(i) + "]" + "\n")
            else:
                out_file.write("#define " + str(molecule_list[i]) + " y[" + str(i) + "]" + "\n")
        out_file.write("\n")
        out_file.write("#define p0 tex2D(param_tex,0,tid)" + "\n")

        # write parameters out in order of reactions
        #parcount = 0
        #for i in range(len(reaction_list)):
        #    for j in range(len(reaction_list[i].param)):
        #        out_file.write("#define " + str(reaction_list[i].param[j]) + " tex2D(param_tex," + str(parcount+1) + ",tid)" + "\n")
        #        parcount += 1

        # write parameters out in order of params list
        #for i in range(len(params)):
        #    out_file.write("#define " + str(params[i]) + " tex2D(param_tex," + str(i+1) + ",tid)" + "\n")

        # write out parameters in numerical order
        for i in range(len(params)):
            out_file.write("#define " + "p" + str(i+1) + " tex2D(param_tex," + str(i+1) + ",tid)" + "\n")
                
        out_file.write("\n")
        out_file.write("struct myFex{" + '\n')
        out_file.write("__device__ void operator()(int *neq, double *t, double *y, double *ydot/*, void *otherData*/){ " + "\n")
        out_file.write("int tid = blockDim.x * blockIdx.x + threadIdx.x;" + "\n")
        for i in range(len(odes)):
            out_file.write('ydot[' + str(i) + "] = p0*(" + odes[i] + ');' + '\n')
        out_file.write("}" + "\n")
        out_file.write("};" + "\n")
        out_file.write("struct myJex{" + "\n")
        out_file.write("__device__ void operator()(int *neq, double *t, double *y, int ml, int mu, double *pd, int nrowpd/*, void *otherData*/){"+"\n")
        out_file.write("return; "+"\n"+"}"+"\n"+"};"+"\n")
        
        return out_file.getvalue()

    def calculateODEs(self,model, rates):
        odes = []
        for sp in range(model.nspecies):
            ode = ''
            for rt in range(len(rates)):
                if model.stoichiometry_matrix.T[sp][rt] > 0:
                    ode +='+'
                    ode += str(model.stoichiometry_matrix.T[sp][rt])
                    ode +=' * '
                    ode += str(rates[rt])
                elif model.stoichiometry_matrix.T[sp][rt] < 0:
                    ode +=' '
                    ode += str(model.stoichiometry_matrix.T[sp][rt])
                    ode +='*'
                    ode += str(rates[rt])
                elif model.stoichiometry_matrix.T[sp][rt] == 0:
                    pass
            odes.append(ode)
        return odes

    def writeMultipleCudaFileODE(self, weaverOutput):
       
        # We are expecting a context returning a list of files

        if getattr(weaverOutput, "getContext", None) != None:
            self.models = weaverOutput.getContext()
        else:
            print " writeMultipleCudaFileODE : getContext() not available. Quitting"
            exit()

        modelStrings = []
        for m in self.models:
            rates = []
            for i in range(m.nreactions):
                rates.append(m.reactions[i].rate)
            odes = self.calculateODEs(m,rates)
            modelStrings.append( self.writeCuda(m.reactions, m.species, m.freeparams, m.nspecies, m.nreactions, odes) )

        return modelStrings


