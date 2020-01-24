# ref: http://www.math-cs.gordon.edu/courses/mat342/handouts/gauss.pdf

import random
import copy
from fractions import Fraction
from time import perf_counter

from matrix import *


def generatematrix(size, max):
    nm = np.zeros((size, size), Fraction)
    nx = np.zeros(size, Fraction)
    for i in range(0, size):
        r = random.randrange(-max, max - 1)
        wsp = r / max
        nx[i] = Fraction(wsp)
        for j in range(0, size):
            v = random.randrange(-max, max - 1)
            wsp = v / max
            nm[i, j] = Fraction(wsp)
    nb = nm @ nx
    return nm, nx, nb


def finishmatrixfloat32(nm1, nx1, nb1):
    nm1 = nm1.astype(np.float32)
    nx1 = nx1.astype(np.float32)
    nb1 = nb1.astype(np.float32)
    nz = Matrix(nm1, nx1, nb1)
    return nz


def finishmatrixfloat64(nm1, nx1, nb1):
    nm1 = nm1.astype(np.float64)
    nx1 = nx1.astype(np.float64)
    nb1 = nb1.astype(np.float64)
    nz = Matrix(nm1, nx1, nb1)
    return nz


def finishmatrixfraction(nm1, nx1, nb1):
    nz = Matrix(nm1, nx1, nb1)
    return nz


def relerror(a1, a2):
    re = (abs(a1 - a2)) / a1
    return re


def getnormofvector(v1):
    sum = 0
    lv1 = len(v1)
    for i in range(0, lv1):
        sum = sum + (v1[i] ** 2)
    sum = sum ** (1 / 2.0)
    return sum


def testfunction(tm):
    t0 = [0, 0, 0]
    e0 = [0, 0, 0]
    tempmatrix = [None] * 3
    for i in range(0, 3):
        tempmatrix[i] = copy.deepcopy(tm)
    t0[0] = perf_counter()
    tempmatrix[0].gauss_nopivoting()
    t0[0] = perf_counter() - t0[0]
    e0[0] = relerror(getnormofvector(tempmatrix[0].oldx), getnormofvector(tempmatrix[0].x))
    print(tempmatrix[0].x)
    t0[1] = perf_counter()
    tempmatrix[1].gauss_partialpivoting()
    t0[1] = perf_counter() - t0[1]
    e0[1] = relerror(getnormofvector(tempmatrix[1].oldx), getnormofvector(tempmatrix[1].x))
    print(tempmatrix[1].x)
    t0[2] = perf_counter()
    tempmatrix[2].gauss_completepivoting()
    t0[2] = perf_counter() - t0[2]
    e0[2] = relerror(getnormofvector(tempmatrix[2].oldx), getnormofvector(tempmatrix[2].x))
    print(tempmatrix[2].x)
    return t0, e0


sizes = [100, 200, 300]

# for i in range(100, 301,5):
#     sizes.append(i)

a = len(sizes)

for i in range(0, a):
    tempma, tempx, tempb = generatematrix(sizes[i], 65536)
    print("DLA ROZMIARU %d X %d" % (sizes[i], sizes[i]))
    # print("")
    # print(".....")
    # print("DLA FLOAT32")
    # tempmat = finishmatrixfloat32(tempma, tempx, tempb)
    # time1, error1 = testfunction(tempmat)
    # print("CZAS" + str(time1))
    # print("BLAD" + str(error1))
    # print(".....")
    print("DLA FLOAT64")
    tempmat = finishmatrixfloat64(tempma, tempx, tempb)
    time1, error1 = testfunction(tempmat)
    print("CZAS" + str(time1))
    print("BLAD" + str(error1))
    print(".....")
    # print("DLA FRACTIONS")
    # tempmat = finishmatrixfraction(tempma, tempx, tempb)
    # time1, error1 = testfunction(tempmat)
    # print("CZAS" + str(time1))
    # print("BLAD" + str(error1))
    # print("====================================================")

