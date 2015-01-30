
from synbioweaver.core import *
import os


class RunABC(Aspect):

    def mainAspect(self):
        self.addWeaverOutput(self.runABC)

    def writeRunFile(self):
        runfile = open('out.sh', 'w')
        runfile.write('export PYTHONPATH=$PYTHONPATH:/home/ucbtle1/abc-sysbio-code/:/home/ucbtle1/cuda-sim-code/' + '\n')
        runfile.write('exe=/home/ucbtle1/abc-sysbio-code/scripts/run-abc-sysbio' + '\n')
        runfile.write('\n')
        runfile.write('python -u $exe --localcode --infile input_file.xml -f -of=results --cuda')
        runfile.close()


    def runABC(self, weaverOutput):
        self.writeRunFile()
        os.chmod('out.sh', 0777)
        os.system('nohup ./out.sh &')