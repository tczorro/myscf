from __future__ import division
from math import pi, exp
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

class Gaussian(object):

    def __init__(self, alpha, nuclei_R=0.0, coeff=1.0):
        self.alpha = alpha
        self.nuclei_R = nuclei_R
        self.coeff = coeff

    @property
    def function(self):
        def fct(ele_r):
            function = (2 * self.alpha / pi) ** (3. / 4) * exp(-self.alpha * (abs(ele_r - self.nuclei_R)) ** 2)
            function *= self.coeff
            return function
        return fct

    def __call__(self, x):
        return self.function(x)

class Slater(object):

    def __init__(self, zeta, nuclei_R=0.0, coeff=1.0):
        self.zeta = zeta
        self.nuclei_R = nuclei_R
        self.coeff = coeff

    @property
    def function(self):
        def fct(ele_r):
            function = (self.zeta ** 3 / pi) ** (1. / 2) * exp(-self.zeta * abs(ele_r - self.nuclei_R))
            function *= self.coeff
            return function
        return fct

    def __call__(self, x):
        return self.function(x)

class Integral(object):

    def overlap_integral(self, function1, function2):
        if function1.nuclei_R == function2.nuclei_R:
            center_R = function1.nuclei_R
        func = lambda x: 4 * pi * (x - center_R) ** 2 * function1(x) * function2(x)
        overlap = integrate.quad(func, center_R, np.inf)
        return overlap

if __name__ == '__main__':
    gau1 = Gaussian(0.270950, 1.3)
    sla1 = Slater(1.0, 1.3)
    gau2 = Gaussian(0.270950, 0)
    sla2 = Slater(1.0, 0)
    func1 = lambda x: gau1(x) * sla1(x) * 4 * pi #* x**2
    func2 = lambda x: gau2(x) * sla2(x) * 4 * pi #* x**2
    print integrate.quad(gau1.function, 1.3, 5.3)
    print integrate.quad(gau2.function, 0., 4)
    print integrate.quad(sla1.function, 1.3, 5.3)
    print integrate.quad(sla2.function, 0., 4)
    print integrate.quad(func1, 1.3, 5.3)
    print integrate.quad(func2, 0., 4)
    integral = Integral()
    result = integral.overlap_integral(gau2, sla2)
    print result
    result2 = integral.overlap_integral(gau1, sla1)
    print result2