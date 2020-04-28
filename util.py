from __future__ import division
from random import random as rnd, randint as rndi
import numpy as np
from math import log, exp, sqrt
import cPickle as pickle


def RANDOM(lb, ub):
    return lb + (ub - lb) * rnd()


def PROD(x):
    y = 1.0
    for xi in x: y *= xi
    return y


def TO_STRING(x):
    return "".join(map(str, x))



def CALC_RELATIVE_ERROR(arr):
    mean = np.mean(arr)
    variance = np.var(arr)
    return sqrt(variance) / mean


def AVERAGE_EVERY(arr, k):
    return [np.mean(arr[i*k:(i+1)*k]) for i in range(int(len(arr)/k))]


def SAVE_OBJECT_BINARY(obj, filename):
	with open(filename, "wb") as output:
		pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def LOAD_OBJECT_BINARY(filename):
	with open(filename, "rb") as input:
		obj = pickle.load(input)
	return obj

