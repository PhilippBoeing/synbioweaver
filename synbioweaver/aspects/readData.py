import numpy as np

class DataMatrix:

    def __init__(self, data_file):
        self.data_file = data_file

    def readData(self):
        data = np.genfromtxt(self.data_file, comments='#', delimiter=' ')
        #time species1 species1sigma species2 species2sigma
        return data