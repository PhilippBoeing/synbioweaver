
from synbioweaver.core import *
import os


class RunCudaSimEmerald(Aspect):

    def mainAspect(self):
        self.addWeaverOutput(self.runCudaSim)

    def writeRunFile(self):
        runfile = open('out.sh', 'w')

        runfile.write('module load cuda')
        runfile.write('module load pycuda')
        runfile.write('module load cuda-sim')
        runfile.write('module load libsbml/5.8')
        runfile.write('#BSUB -o log.synbioweaver')
        runfile.write('#BSUB -e err.synbioweaver')
        runfile.write('#BSUB -W 30:00')
        runfile.write('export PYTHONPATH=$PYTHONPATH:/home/ucl/eisuc058/code/abc-sysbio' + '\n')
        runfile.write('export PYTHONPATH=$PYTHONPATH:/home/ucl/eisuc058/code/cuda-sim' + '\n')

        runfile.write('exe=/home/ucl/eisuc058/code/abc-sysbio/scripts/run-abc-sysbio' + '\n')
        runfile.write('\n')
        runfile.write('python -u $exe --localcode  -S --infile input_file.xml -f -of=results_simulation --cuda')
        runfile.close()


    def runCudaSimEmerald(self, weaverOutput):
        self.writeRunFile()
        os.chmod('out.sh', 0777)
        os.system('bsub out.sh')
