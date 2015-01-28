from synbioweaver.core import *
import numpy
#from synbioweaver.aspects.reactionDefinitions import *
from synbioweaver.aspects.printReactionNetworkAspect import *
from synbioweaver.aspects.WriteInputFile import *


class WriteCudaFile(Aspect):

    def mainAspect(self):
        self.addWeaverOutput(self.writeCudaFile)

    def writeCuda(self, stoichiometry_matrix, reaction_list, molecule_list, rates_list, params, nspecies, nreactions):
        self.cuda_file = str('model.cu')
        self.name = 'model'
        out_file = open(self.cuda_file,'w')
        out_file.write("#define NSPECIES " + str(nspecies) + "\n")
        out_file.write("#define NPARAM " + str(nreactions) + "\n")
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

    def writeCudaFile(self, weaverOutput):
        self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getReactions()

        #self.reaction_list = []
        self.rates = []
        #self.params = []

        for i in range(0, self.nreactions):
            self.rates.append(self.reactions[i].rate)
            #self.params.append(self.reactions[i].param)
            #self.reaction_list.append(self.reactions[i])
        cuda_file = self.writeCuda(self.stoichiometry_matrix, self.reactions, self.species, self.rates, self.params, self.nspecies, self.nreactions)
        #return cuda_file
