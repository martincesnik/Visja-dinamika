import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# string properties

l = 10  # length
mass = 1.  # mass
P = 10000  # pretension force
c = np.sqrt(P*l/mass)


# time discr.--------------------------------------------

t_end = 1  # total time
dt = 0.0001  # time increment
t = np.arange(start=0, stop=t_end+dt, step=dt)  # time array
m = len(t)

# # space diskr. -----------------------------------------

n = 100
dx = l/n  # space discr.
x = np.arange(start=0, stop=l+dx, step=dx)  # space array


w = np.zeros((m, n+1))

# initial state
##
#l_1 = 0.25 * l
#l_2 = 0.5 * l
#l_3 = 0.25 * l
#w_0 = 0.1
#
#w_init_1 = np.arange(0, l_1, dx) * w_0/l_1
#w_init_2 = - np.arange(0, l_2, dx) * 2*w_0/l_2 + w_0
#w_init_3 = np.arange(0,l_3+dx, dx) * w_0/l_3-w_0
##print(y_init_3)
#
#w_init = np.concatenate((w_init_1, w_init_2, w_init_3))
#w[0,1:-1] = w_init[1:-1] # preserve zeros at edges

#####################################################
#l_1 = 0.5 * l
#l_2 = 0.5 * l
#w_0 = 0.1
#
#w_init_1 = np.arange(0, l_1, dx) * w_0/l_1
#w_init_2 = - np.arange(0, l_2+dx, dx) * w_0/l_2 + w_0
#
#w_init = np.concatenate((w_init_1, w_init_2))
#w[0,1:-1] = w_init[1:-1] # preserve zeros at edges

####################################
w[0,0]=1
####################################

plt.plot(w[0])
plt.show()
# first time increment
for i in np.arange(start=1, stop=n):

    w[1,i] = w[0,i]+1/2*(dt*c/dx)**2*(w[0,i+1]-2*w[0,i]+w[0,i-1])


# second and higher increment
for j in np.arange(start=1, stop=m-1):
    for i in np.arange(start=1, stop=n):
        w[j+1, i] = (dt*c/dx)**2*(w[j,i+1]-2*w[j,i]+w[j,i-1])+2*w[j,i]-w[j-1,i]


print('done calculating')


animate=True

if animate:


    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, l*1), ylim=(-0.1, 0.1))
    ax.grid()

    line, = ax.plot([], [], 'o-', lw=2)

    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    def init():
        line.set_data([], [])

        time_text.set_text('')
        return line, time_text

    def animate(i):
        thisx = x
        thisy = w[i,:]

        line.set_data(thisx, thisy)

        time_text.set_text(time_template%(t[i]))
        return line, time_text

    ani = animation.FuncAnimation(fig, animate, np.arange(1, len(t)),
        interval=1, blit=False, init_func=init)

    plt.show()





