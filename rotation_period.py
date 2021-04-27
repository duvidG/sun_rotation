#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 17:34:45 2021

@author: dgiron
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit as cf
import math
from uncertainties import ufloat, unumpy
from uncertainties.umath import *
from uncertainties.unumpy import nominal_values as nv
from uncertainties.unumpy import std_devs as st

def ang_speed(x, y, lb):
    """
    Obtains the angular frequency (in deg/day) for a set of longitudes and times for a 
    given sunspot.
    
    ----------
    x : numpy.ndarray
        Day (first day is counted as 0)
    y : numpy.ndarray
        Longitude (degrees)
    lb : list
        Labels to identify the sunspot in the graph 

    Returns
    -------
    pol[0]: float
        Slope of the fit (angular speed)
    err[0]: float
        Uncertainty of the slope
    """
    pol, err = np.polyfit(x, nv(y), deg=1, cov=True, w=1/st(y))
    err = np.sqrt(np.diag(err))
    
    xx = np.linspace(x[0], x[-1])
    yy = np.polyval(pol, xx)
    
    plt.errorbar(x, nv(y), yerr=st(y), fmt='k.')
    plt.plot(xx, yy, label=lb)
    plt.xlabel(r'$t$/days')
    plt.ylabel(r'$\lambda$/deg')
    
    return pol[0], err[0]

def tabla_latex(tabla, ind, col, r):
    """
    Prints an array in latex format
    Args:
        tabla: array to print in latex format (formed as an array of the variables (arrays) that want to be displayed)
        ind: list with index names
        col: list with columns names
        r: number of decimals to be rounded
        
    Returns:
        ptabla: table with the data in a pandas.Dataframe format
    """
    tabla = tabla.T
    tabla = tabla.round(r)
    ptabla = pd.DataFrame(tabla, index=ind, columns=col)
    print("Tabla en latex:\n")
    print(ptabla.to_latex(index=False, escape=False))
    return ptabla


def main():
    datos = np.genfromtxt('datos2.txt', delimiter=',')
    datos = datos[1:]

    day = datos[:, 1]
    long = datos[:, 2]
    lat = datos[:, 3]
    letra = ['B', 'C', 'D', 'E']
    err = 0.7
    long_err = unumpy.uarray(long, err)
    
    day_rel = day - day[0] # Setting the first day as day 0
    
    plt.clf()
    
    # Period for the phi=30 sunspot
    
    pte_30, err_30 = ang_speed(day_rel[0:6], long_err[0:6], r'A: $\phi$ = '+str(np.round(np.mean(lat[0:6]), 3)))
    
    # With uncertainties package a variable is defined as var = (value, err), and error is propagated in every 
    # operation of this variable (it also works with arrays). nv() and st() fucntions return the value and the
    # uncertainty respectively.
    
    pte_30 = ufloat(pte_30, err_30)
    print('Non-equatorial sunspots:')
    print('Synodic period {:.3f}' .format(360 / pte_30))
    
    s_30 = 360/pte_30
    p_30 = (s_30 * 365.25) / (s_30 + 365.25)
    
    print('Sidereal period {:.3f}' .format(p_30))
    
    x = []
    err = []
    pte = []
    c = 6
    j = np.round(np.array([np.round(np.mean(lat[0:6]), 3), pte_30, s_30, p_30]), 2)
    print(j)
    
    # Period for every equatorial sunspot
    print('\nEquatorial sunspots:')
    for i in range(4):
        y = ang_speed(day_rel[c:c+6], long_err[c:c+6], letra[i]+r': $\phi$ = '+str(np.round(np.mean(lat[c:c+6]), 3)))
        pte.append(y[0])
        err.append(y[1])
        x.append(np.mean(lat[c:c+6]))
        c = c + 6
        
    plt.legend()
    plt.savefig('ajuste.png', dpi=720)
    plt.show()
    
    pte = np.array(pte)
    err = np.array(err)
    
    pte = unumpy.uarray(pte, err)
    s = 360 / pte
    x = np.round(x, 3)
    
    p =  (s  * 365.25) / (s + 365.25)
    
    tab2 = np.array([x, nv(pte), st(pte), nv(s), st(s), nv(p), st(p)])
    ptabla2 = tabla_latex(tab2, ['' for i in tab2[1]], ['$\phi$/deg','pte', 'err', '$\S$/day', '$\delta$ S','$\P$/day', '$\delta P$'], 3)
    print(ptabla2)
    
    # Mean synodic and sideral
    
    s_mn = np.mean(s)
    p_mn = np.mean(p)
    
    
    print('\nMean equatorial sunspots:')
    print('Synodic period {:.3f}' .format(s_mn))
    
    print('Sidereal period {:.3f}' .format(p_mn))
    
    # Prints a table with the data (also in latex format)
    
    tab = np.array([day_rel, lat, long])
    ptabla = tabla_latex(tab, ['' for i in tab[1]], ['$t$/day', '$\phi$/deg', '$\lambda$/deg'], 3)
    display(ptabla)
 
main()




