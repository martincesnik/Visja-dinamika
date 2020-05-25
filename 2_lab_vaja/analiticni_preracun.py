# -*- coding: utf-8 -*-
"""
@author: Martin Česnik
"""

from sympy import symbols, diff, Function, sin, cos, Derivative, simplify, expand, Matrix
import numpy as np
import matplotlib.pyplot as plt


# define variable symbols for deduction ===============================================

t = symbols('t')    	                     # time variable
L, k, m_s = symbols('L, k, m_s')             # pole length, spring stiffness and spring mass
J_1, J_2, J_3 = symbols('J_1, J_2, J_3')     # MMI's of poles, exact values will be substituted later


# system's spatial properties ================================================

# values for variables in deduction
L_inp = 0.38            # pole length
k_inp = 25734           # spring stiffness
m_s_inp = 0.105         # spring mass

# other values
M = 1.513               # weight mass
m = 0.536               # pole mass
m_a = 0.052             # adapter mass
R = 0.0352              # weight outer radius
r = 0.0075              # weight inner radius
H = 0.04                # weight height

L_i = np.array([0.218, 0.088, 0.28])    # distance from weight mass centre to rotation axis

# define mass moments of inertia (MMI) ---------------------------------

J_p =1/3 * m*L_inp**2               # pole MMI
J_a = m_a*L_inp**2                  # adapter MMI
    
J_i_inp = np.zeros_like(L_i)

for i in range(len(L_i)):
    J_w = M*L_i[i]**2 + 1/12*M*( 3*(R**2+r**2) + H**2)
    
    J_i_inp[i] = J_p + J_w + J_a         # pole + weight + adapter MMI

    
# define generalized coordinates as functions =================================

phi_1 = Function('phi_1')
phi_2 = Function('phi_2')
phi_3 = Function('phi_3')

x_1 = Function('x_1')
x_2 = Function('x_2')
x_3 = Function('x_3')

# define Lagrangian ==========================================================

        
E_k = (J_1/2 * diff(phi_1(t),t)**2 +         # kinetic energy due to poles' rotation
        J_2/2* diff(phi_2(t),t)**2 +
        J_3/2* diff(phi_3(t),t)**2)

E_k +=  (diff(x_1(t),t)**2 + 
         (diff(x_1(t),t) - diff(x_2(t),t))**2 + 
         (diff(x_2(t),t)-diff(x_3(t),t))**2) *m_s/6     # kinetic energy due to spring velocity


    
E_p = (1/2*k*x_1(t)**2 +                    # potential energy of springs, gravitational potential energy is relatively small and is therefore neglected 
        1/2*k*(x_1(t)-x_2(t))**2 +  
        1/2*k*(x_2(t)-x_3(t))**2)
            
lagrangian = E_k - E_p                   # lagrangian


# include x_i = L*sin(phi_i) and do substitution in lagrangian    
lagrangian = lagrangian.subs([(x_1(t), L*sin(phi_1(t))), 
                              (x_2(t), L*sin(phi_2(t))), 
                              (x_3(t), L*sin(phi_3(t)))]).doit()

# linearize lagrangian under small rotation assumption
# if linearization is made after lagrangian derivation, the kinetic energy part of springs result in non-linear system of equations
assumptions = [(sin(phi_1(t)), phi_1(t)), (cos(phi_1(t)), 1),
                (sin(phi_2(t)), phi_2(t)), (cos(phi_2(t)), 1),
                (sin(phi_3(t)), phi_3(t)), (cos(phi_3(t)), 1)]

lagrangian = lagrangian.subs(assumptions)


# deduce equations of motion =================================================


eq_motion_1 = diff( diff( lagrangian, diff(phi_1(t),t) ) ,t) - diff(lagrangian, phi_1(t))
eq_motion_2 = diff( diff( lagrangian, diff(phi_2(t),t) ) ,t ) - diff(lagrangian, phi_2(t))
eq_motion_3 = diff( diff( lagrangian, diff(phi_3(t),t) ) ,t ) - diff(lagrangian, phi_3(t))

# construct mass and stiffness matrix =========================================

m_1, m_2, m_3 = [], [], []
k_1, k_2, k_3 = [], [], []

for mass_row, stiff_row, eq_motion in zip([m_1, m_2, m_3], [k_1, k_2, k_3], [eq_motion_1, eq_motion_2, eq_motion_3]):

    [mass_row.append(simplify(eq_motion).expand(eq_motion_1).coeff(diff(gen_coord,t,t)) )
            for gen_coord in [phi_1(t), phi_2(t), phi_3(t)]]
    
    [stiff_row.append(simplify(eq_motion).expand(eq_motion_1).coeff(gen_coord) )
            for gen_coord in [phi_1(t), phi_2(t), phi_3(t)]]
            

M = Matrix([m_1, m_2, m_3])
K = Matrix([k_1, k_2, k_3])


# substitute variables in M and K for values
M = M.subs([(L, L_inp), (m_s, m_s_inp), (J_1, J_i_inp[0]),
         (J_2, J_i_inp[1]), (J_3, J_i_inp[2])])
    
K = K.subs([(L, L_inp), (k, k_inp)])


# get natural frequencies and modeshapes ======================================
A = M**-1*K

A_n = np.array(A).astype(np.float64)

eigvals, eigvects = np.linalg.eig(A_n)

idx = eigvals.argsort()
w_0 = np.sqrt(eigvals[idx])/2/np.pi
Phi = eigvects[:,idx]/eigvects[:,idx][0,:]

print('System\'s natural frequencies')
print(w_0)
print('')
print(Phi)
