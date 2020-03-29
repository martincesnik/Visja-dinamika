# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 21:21:56 2020

@author: martin
"""

import numpy as np
import matplotlib.pyplot as plt


def get_F_exact (F_0, t, tau, show = False):
    F = np.ones_like(t_array)
    
    for (t, i) in zip(t_array, range(len(F))):
        if t % tau < tau/2:
            F[i] = F_0
        else:
            F[i] = -F_0
    
    if show:        
        plt.plot(t_array,F)
        plt.show()
        
    return F

def get_modal_params(m,d,k):
    
    w_0 = np.sqrt(k/m)
    delta = d/2/m/w_0
    
    return w_0, delta


def get_four_coeffs(F_0, n_max, compare_exact=True, t=None, F_exact=None, w=None):
    
    a_j=np.zeros(n_max)
    b_j=np.zeros(n_max)
    
    for j in range(n_max):
        
        if j==0: #a_0
            a_j[j] = 0
        
        else:
            a_j[j] = 0
            b_j[j] = 2*F_0/j/np.pi *(1-np.cos(j*np.pi))
            
    if compare_exact:
        
        F_approx = np.ones_like(F_exact) 
        
        for j in range(n_max):
            if j==0:
                F_approx *= a_j[j]/2 
            else:
                F_approx += a_j[j]*np.cos(j*w*t) + b_j[j]*np.sin(j*w*t)
        
        plt.plot(t,F_exact)
        plt.plot(t, F_approx)
        plt.show()
                
    return a_j, b_j





if __name__=='__main__':
    
    # input values
    m = 1           # [kg]
    k = 30*10**3    # [N/m]
    d = 1           # [Ns/m]
    F_0 = 10        # [N]
    tau = np.pi/10  # [s]
    
    w = 2*np.pi/tau 
    
    
    # time array 
    t_min = 0       # [s]
    t_max = 1       # [s]
    dt = 0.0001     # [s]
    t_array = np.arange(t_min,t_max,dt)
    
    # exact excitation force
    F_exact = get_F_exact(F_0, t_array, tau, show=True)
    
    # obtain modal parameters
    w_0, delta = get_modal_params(m,d,k)
    
    # get fourier series
    n_max = 10
    a_j, b_j= get_four_coeffs(F_0, n_max, True, t_array, F_exact, w)
    
    # the calculation of response using superposition principle is intentionally not coded yet
