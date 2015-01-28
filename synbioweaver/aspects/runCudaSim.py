
from synbioweaver.core import *
import os


class RunCudaSim(Aspect):

    def mainAspect(self):
        self.addWeaverOutput(self.runCudaSim)

    def writeRunFile(self):
        runfile = open('out.sh', 'w')
        runfile.write('export PYTHONPATH=$PYTHONPATH:/home/ucbtle1/abc-sysbio-code/:/home/ucbtle1/cuda-sim-code/' + '\n')
        runfile.write('exe=/home/ucbtle1/abc-sysbio-code/scripts/run-abc-sysbio' + '\n')
        runfile.write('\n')
        runfile.write('python -u $exe --localcode  -S --infile input_file.xml -f -of=results_simulation --cuda')
        runfile.close()


    def runCudaSim(self, weaverOutput):
        self.writeRunFile()
        os.chmod('out.sh', 0777)
        os.system('nohup ./out.sh &')