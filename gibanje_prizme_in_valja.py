'''
Created on Mar 23, 2015
author: Martin Cesnik
'''

# verified on interpreter python 3.7

# imports
import numpy as np
from numpy.linalg import inv
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def vectorfield(t, y, m_A, m_B, alpha, mu, g):

# Vectorfield function takes a state space vector as input and returns its derivative
# at given time instance. Arbitrary arguments can be passed for the derivative calculation
#
# inputs
# y:        state-space vector; in the current case of prism and cylinder it consists of
#           y = [y[0], y[1], y[2], y[3]] = [x, dx/dt, s, ds/dt] 
# t:        time instance; this is automatically provided by odeint function
# params:   a tuple of additional arguments
#
# outputs
# dy:       a derivative of a state-space vector; in the current case it consists of
#           dy = [dx/dt, d^2x/dt^2, ds/dt, d^2s/dt^2]

    
    #define matrix A according to devised equilibrium eqs
    A = np.array([[m_A + m_B, m_B*(np.cos(alpha)-mu*np.sin(alpha))],
                  [m_B*np.cos(alpha), 3./2*m_B]])
    
    # define matrix B according to devised equilibrium eqs
    # np.sign(dx/dt) is used for correct directionality of friction force vector 
    B = np.array([[-np.sign(y[1])*mu*g*(m_A+m_B)],
                 [m_B*g*np.sin(alpha)]])
    
    # A inversion
    A_inv = inv(A)

    # calculation of state-space vector derivative for given parameters
    dy = [y[1], 
        np.dot(A_inv,B)[0],
        y[3],
        np.dot(A_inv,B)[1]]
    
    return dy





#define data
g = 9.81
m_A = 1000
m_B = 1000
mu = 0. # if friction coefficient exceeds certain value, the equilibrium equations do not hold
alpha = 25*np.pi/180

#pack constants into tuple
params=(m_A, m_B, alpha, mu, g)

# set initial values 
x_0 = 0.0
dx_dt_0 = 0.0
s_0 = 0.0
ds_dt_0 = 0.0

init_vals = [x_0, dx_dt_0, s_0, ds_dt_0] 

# set time array
dt = 0.001
t_0 = 0.0
t_f = 1.5
t_points = np.arange(t_0, t_f, dt)


# solve a set of first order linear differential equations, as defined in vectorspace
result = solve_ivp(vectorfield,[t_0,t_f], init_vals, args=(params), t_eval=t_points)


# get rigid body displacements
t = result.t
x = result.y[0,:]
s = result.y[2,:]

plot = True
animate = True

if plot:
    
    # only a simple plot
    plt.plot(t,x)
    plt.plot(t,s)
    plt.show()

if animate:
        
    # basic animation
    
    #additional body dimensions
    r = 0.3
    l = 2
    h = l*np.tan(alpha)
    
    # displacement arrays of prism corners
    x_p_1 = x
    y_p_1 = np.zeros_like(x)
    
    x_p_2 = l + x
    y_p_2 = np.zeros_like(x)
    
    x_p_3 = x
    y_p_3 = np.zeros_like(x)+h 
    
    # displacement arrays of cylinder center and its rotation
    x_c = x+r*np.sin(alpha)+s*np.cos(alpha)
    y_c = np.zeros_like(x)+h+r*np.cos(alpha)-s*np.sin(alpha)
    phi = s/r
    
    # find a time instance, when a cylinder hits the ground
    end_frame = np.argmax(y_c<r)
    
    
    #  some prerequisites for animation plot, for detailed explanation search teh internets
    fig = plt.figure()
    ax = plt.axes(xlim=(-2,3), ylim=(-1, 3),aspect='equal')
    
    line, = ax.plot([], [], lw=2)
    prism = plt.Polygon([[0,0], [0,0], [0,0]])
    cylinder = plt.Circle((0,0), r, fc='y')

    ax.add_patch(prism) 
    ax.add_patch(cylinder)
    
    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    
    def init():
        line.set_data([], [])
        time_text.set_text('')
        prism.set_visible(False)
        cylinder.set_visible(False)
        ground, = ax.plot([-2,3], [-0.02,-0.02], lw=3, c='k')
        
        return prism, cylinder, line, ground, time_text
    
    def animate(i):
        
        if i == 1:
            prism.set_visible(True)
            cylinder.set_visible(True)
        
        x_line = [x_c[i], x_c[i]+r*np.sin(phi[i])]
        y_line = [y_c[i], y_c[i]+r*np.cos(phi[i])]
        
        line.set_data(x_line, y_line)
        
        time_text.set_text(time_template%(i*dt))
        
        prism.set_xy([[x_p_1[i], y_p_1[i]],
                      [x_p_2[i], y_p_2[i]],
                      [x_p_3[i], y_p_3[i]]])

        cylinder.center= ((x_c[i], y_c[i]))
        
        return prism, cylinder, line, time_text
    
    ani = animation.FuncAnimation(fig, animate, np.arange(1,end_frame),
        interval=dt, blit=True, init_func=init)
    
    # save maybe?
    #ani.save('D:\primer.mp4', fps=15)
    
    plt.show()
