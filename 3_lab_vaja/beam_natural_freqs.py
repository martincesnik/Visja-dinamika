
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


def BC_free_free(beta_L):

    return np.cos(beta_L)*np.cosh(beta_L)-1


def BC_free_mass(beta_L, L, E, I, m, c):
    
    beta = beta_L/L
    
    return (E*I*np.cos(beta_L)*np.cosh(beta_L) -
            E*I - 
            beta*c**2*m*np.sin(beta_L)*np.cosh(beta_L) + 
            beta*c**2*m*np.cos(beta_L)*np.sinh(beta_L))

def BC_free_mass_MMI(beta_L, L, E, I, m, J, c):
    
    beta = beta_L/L
    
    return (-E**2*I**2*np.cos(beta_L)*np.cosh(beta_L) + 
            E**2*I**2 + 
            E*I*J*beta**3*c**2*np.sin(beta_L)*np.cosh(beta_L) + 
            E*I*J*beta**3*c**2*np.cos(beta_L)*np.sinh(beta_L) + 
            E*I*beta*c**2*m*np.sin(beta_L)*np.cosh(beta_L) - 
            E*I*beta*c**2*m*np.cos(beta_L)*np.sinh(beta_L) + 
            J*beta**4*c**4*m*np.cos(beta_L)*np.cosh(beta_L) + 
            J*beta**4*c**4*m)

if __name__=='__main__':
    
    # define geometry and material properties
    
    L = 0             # beam length [m] 
    rho = 0          # beam density [kg/m^3]
    a = 0           # beam width [m]
    b = 0           # beam height [m]
    E = 0    # Young modulus [Pa]
    a_w = 0          # weight height [m]
    b_w = 0         # weight width [m]
    c_w = 0         # weight length [m]
    
    m = a_w * b_w * c_w * rho
    
    A_cs = a*b
    
    I = a*b**3/12
    J = m*(a_w**2+c_w**2)/12
    c = np.sqrt(E*I/rho/A_cs)
    
    beta_L = np.linspace(0.01, 40, 100000) 

    # =============================================================================
    #  plotting     
    # =============================================================================
    
    # plot characteristic determinant
    # plt.plot(beta_L, BC_free_free(beta_L), label = 'free-free')
    # plt.plot(beta_L, BC_free_mass(beta_L, L, E, I, m, c), label = 'free-mass')
    # plt.plot(beta_L, BC_free_mass_MMI(beta_L, L, E, I, m, J, c), label = 'free-mass-MMI')
    # plt.yscale('symlog')
    # plt.xlabel(r'$\beta L$')
    # plt.ylabel('characteristic determinant')
    # plt.tight_layout()
    # plt.legend()
    # plt.grid()
    # plt.show()  


    # =============================================================================
    #     free-free boundary conditions
    # =============================================================================
        
    ranges = [2.5, 6, 9, 12, 15, 18, 22,25, 28]
    ranges_MMI = [2.5, 6, 8, 11, 13, 16, 19,22, 25]
    
    roots = np.zeros(len(ranges)-1)
    
    for i in np.arange(len(roots)):
        roots[i] = brentq(BC_free_free, ranges[i], ranges[i+1])
        
    f_0i = (roots/L)**2*c/2/np.pi
    # print(roots)
    
    print('Natural frequencies [Hz] for free-free BC:')
    print(f_0i)

    # =============================================================================
    #     free-mass boundary conditions
    # =============================================================================
    
    roots = np.zeros(len(ranges)-1)
    for i in np.arange(len(roots)):
        roots[i] = brentq(BC_free_mass, ranges[i], ranges[i+1], args=(L, E, I, m, c))
        
    f_0i = (roots/L)**2*c/2/np.pi
    # print(roots)
    
    print('\nNatural frequencies [Hz] for free - mass BC:')
    print(f_0i)
    
    # =============================================================================
    #     free - mass+MMI boundary conditions
    # =============================================================================
    
    roots = np.zeros(len(ranges)-1)
    for i in np.arange(len(roots)):
        roots[i] = brentq(BC_free_mass_MMI, ranges_MMI[i], ranges_MMI[i+1], args=(L, E, I, m, J, c))
        
    f_0i = (roots/L)**2*c/2/np.pi
    # print(roots)
    
    print('\nNatural frequencies [Hz] for free - mass+MMI BC:')
    print(f_0i)
