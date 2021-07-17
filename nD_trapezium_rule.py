# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 15:03:36 2021

@author: jakes
"""

import sympy as smp
smp.init_printing()


def oneDReduction(function, variableList, boundsList, nList, currentDim): # boundsList is a len n list of bound pairs, nList is a len n list

    ## variableList, boundsList, nList should be ordered with the outer-most values first (ie. the ones to do first go last)

    k1, k2, k3, k4, k5, k6, k7 = smp.symbols('k1, k2, k3, k4, k5, k6, k7')
    
    kList = [k1, k2, k3, k4, k5, k6, k7]
    
    #k = smp.symbols(f'k_{currentDim}')
    
    k = kList[currentDim - 1]
    
    bounds = boundsList[-1]
    n = nList[-1]
    variable = variableList[-1]
    
    a = bounds[0]
    b = bounds[1]
    
    ## Need to include handling for if the function does NOT explicitly depend on one of the variables

    approxPre = (b - a) / n
    approx1 = (function.subs(variable, a) + function.subs(variable, b)) / 2
    approx2 = function.subs(variable, (a + k) * ((b - a) / n))

    approx = approxPre * (approx1 + smp.Sum(approx2, (k, 1, n - 1)))
    '''
    print ()
    smp.pprint (approx)
    print ()
    '''
    
    currentDim -= 1
    
    if currentDim == 0:
        print (f"Reached dimension {currentDim}")
        
        print (approx.doit())
        
        return type(approx)
    
    else:
        print (f"Now working on dimension {currentDim}")
        
        variableList = variableList[0:currentDim]
        print (variableList)
        boundsList = boundsList[0:currentDim]
        print (boundsList)
        nList = nList[0:currentDim]
        print (nList)
        
        oneDReduction(approx, variableList, boundsList, nList, currentDim)


def main():

    x, y, z = smp.symbols('x y z')
        
    s, T, m, v, k = smp.symbols('sigma T m_a v k')
    
    exponent = -1 * (m * x**2) / (2 * k * T)
    norm = 4 * smp.pi * (m / (2 * smp.pi * k * T)) ** (3/2)
    maxwell =  norm * smp.exp(exponent) * x**2 
    
    numbersMaxwell = maxwell.subs([(T, 1E1), (k, 1.38E-23), (m, 1.67E-27)])
    
    # function = numbersMaxwell * numbersMaxwell.subs(x, y) * numbersMaxwell.subs(x, z)
    
    function = x * y # * z
    
    bounds = [[0, 10], [0, 10]]
    n = [50, 50]
    variables = [y, x]
    
    
    exact = smp.Integral(function, (x, 0, 10), (y, 0, 10))
    
    print ("Exact:")
    smp.pprint (exact.doit())
    print ()
    
    
    fullSum = oneDReduction(function, variables, bounds, n, 2)
    
    print ()
    print ("Approximate:")
    smp.pprint(fullSum)


if __name__ == "__main__":
    
    main()
