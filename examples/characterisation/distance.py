from numpy import *
from numpy.linalg import *

from abcsysbio import euclidian

# read in the data, should only be done once when imported
# we extract only the median of the two observables
data = genfromtxt("data-res-growth.txt", delimiter=" ", skip_header=1)[:,(2,5)]


def distance(data1, data2, parameters, model):
    # data1 is simulated, and has shape npoints x beta
    # data2 is real
    # model is the model number

    d =  euclidian.euclidianDistance(data1, data, parameters, model)
    return d
