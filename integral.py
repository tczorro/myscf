from __future__ import division
from math import pi, exp
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

def gaussian(alpha, coeff=1.0, nuclei_R=0.0):
    def fct(ele_r):
        function = (2 * alpha / pi) ** (3 / 4) * exp(-alpha * (abs(ele_r - nuclei_R)) ** 2)
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

if __name__ == '__main__':
    gaus1 = gaussian(0.109818, 0.444635)
    gaus2 = gaussian(0.405771, 0.535328)
    gaus3 = gaussian(2.22766, 0.154329)


    all = polynomial(gaus1, gaus2, gaus3)

    x_list = np.arange(0,10,0.01)
    y_list = [all(i) for i in x_list]
    plt.plot(x_list, y_list,'ro')
    plt.show()