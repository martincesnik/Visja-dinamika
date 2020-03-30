# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 22:09:44 2020

@author: marti
"""

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
        if (t % tau < tau/6 or t % tau > 5*tau/6):
            F[i] = F_0
        elif (t % tau > tau/6 and t % tau < tau/2):
            F[i] = F_0/2
        else:
            F[i] = 0
            
    if show:        
        plt.figure()
        plt.plot(t_array,F)
        plt.grid()
        plt.xlabel('Time [s]')
        plt.ylabel('Force [N]')
        plt.show()
        
    return F

def get_modal_params(m,d,k):
    
    w_0 = np.sqrt(k/m)
    delta = d/2/m/w_0
    
    return w_0, delta


def get_four_coeffs(F_0, n_max):
    
    a_j=np.zeros(n_max)
    b_j=np.zeros(n_max)
    
    for j in range(n_max):
        
        if j==0: #a_0
            a_j[j] = F_0
        
        else:
            a_j[j] = 3*F_0/(2*np.pi*j)*np.sin(j*np.pi/3)
            b_j[j] = F_0/(2*np.pi*j)*(np.cos(j*np.pi/3)-np.cos(j*np.pi))
            
    return a_j, b_j


    
def exact_vs_approx(t, F_exact, a_j, b_j, w):

    F_approx = np.ones_like(F_exact) 
    
    for j in range(len(a_j+1)):
        if j==0:
            F_approx *= a_j[j]/2 
        else:
            F_approx += a_j[j]*np.cos(j*w*t) + b_j[j]*np.sin(j*w*t)
    
    plt.plot(t,F_exact)
    plt.plot(t, F_approx)
    plt.xlabel('Time [s]')
    plt.ylabel('Force [N]')
    plt.grid()
    plt.show()






if __name__=='__main__':
    
    # input values
    m = 1           # [kg]
    k = 30*10**3    # [N/m]
    d = 2           # [Ns/m]
    F_0 = 10        # [N]
    tau = 2  # [s]
    
    w = 2*np.pi/tau 
    
    
    # time array 
    t_min = 0       # [s]
    t_max = 3*tau       # [s]
    dt = 0.0001     # [s]
    t_array = np.arange(t_min,t_max,dt)
    
    # exact excitation force
    F_exact = get_F_exact(F_0, t_array, tau)
    
    # obtain modal parameters
    w_0, delta = get_modal_params(m,d,k)
    
    # get fourier series
    n_max = 100
    
    a_j, b_j= get_four_coeffs(F_0, n_max)
    
    exact_vs_approx(t_array, F_exact, a_j, b_j, w)    # the calculation of response using superposition principle is intentionally not coded yet
    