from __future__ import division
from math import pi, exp
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

def gaussian(alpha, coeff=1.0, nuclei_R=0.0):
    def fct(ele_r):
        function = (2 * alpha / pi) ** (3. / 4) * exp(-alpha * (abs(ele_r - nuclei_R)) ** 2)
        function *= coeff
        return function
    return fct

def polynomial(*functions):
    def poly(ele_r):
        fct = 0
        for fct_i in functions:
            fct += fct_i(ele_r)
        return fct
    return poly

def slater(zeta, coeff=1.0, nuclei_R=0.0):
    def fct(ele_r):
        function = (zeta ** 3 / pi) ** (1. / 2) * exp(-zeta * abs(ele_r - nuclei_R))
        function *= coeff
        return function
    return fct

def overlap_integral(func1, func2):
    func = lambda x: 4 * pi * x**2 * func1(x) * func2(x)
    overlap = integrate.quad(func, 0., np.inf)
    return overlap

if __name__ == '__main__':
    gaus1 = gaussian(0.109818, 0.444635)
    gaus2 = gaussian(0.405771, 0.535328)
    gaus3 = gaussian(2.22766, 0.154329)

    slater1 = slater(1.0, 1. ,0.)
    all = polynomial(gaus1, gaus2, gaus3)

    x_list = np.arange(0.,4,0.01)
    y_list = [all(i) for i in x_list]
    y2_list = [slater1(i) for i in x_list]

    another_gaus = gaussian(0.270950)
    plt.plot(x_list, y_list)
    plt.plot(x_list, y2_list)
    plt.show()
    print overlap_integral(slater1, another_gaus)
    print overlap_integral(slater1, slater1)
    # overlap_function = lambda x: (slater1(x) * all(x))
    # print integrate.quad(slater1, 0, np.inf)
