# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 20:49:21 2021

@author: jakes
"""

import sympy as smp
smp.init_printing()


def oneDReduction(function, variable, bounds, n, currentDim): # bounds is a len 2 list

    k1, k2, k3, k4, k5 = smp.symbols('k1 k2 k3 k4 k5')
    
    k = [k1, k2, k3, k4, k5]
    
    dummy = k[currentDim]

    a = bounds[0]
    b = bounds[1]
    

    approxPre = (b - a) / n
    approx1 = (function.subs(variable, a) + function.subs(variable, b)) / 2
    approx2 = function.subs(variable, (a + dummy) * ((b - a) / n))

    approx = approxPre * (approx1 + smp.Sum(approx2, (dummy, 1, n - 1)))

    return approx

def main():

    x, y, a, b = smp.symbols('x y a b')
        
    s, T, m, v, k = smp.symbols('sigma T m_a v k')
    
    exponent = -1 * (m * x**2) / (2 * k * T)
    norm = 4 * smp.pi * (m / (2 * smp.pi * k * T)) ** (3/2)
    maxwell =  norm * smp.exp(exponent) * x**2 
    
    numbersMaxwell = maxwell.subs([(T, 1E1), (k, 1.38E-23), (m, 1.67E-27)])
    
    function = numbersMaxwell * numbersMaxwell.subs(x, y)
    
    bounds = [0, 5000]
    n = 100000
    
    exact = smp.Integral(function, (x, 0, 5000), (y, 0, 5000))
    
    oneD = oneDReduction(function, x, bounds, n, 2)
    
    fullSum = oneDReduction(oneD, y, bounds, n, 1)
    
    print ("Exact:")
    smp.pprint (exact.doit().evalf())
    print ()
    print ("Approximate:")
    smp.pprint(fullSum.doit().evalf())


if __name__ == "__main__":
    
    main()