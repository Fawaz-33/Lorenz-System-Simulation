import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from functools import partial
from itertools import count
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
from mpl_toolkits import mplot3d

# Lorenz 
sigma = 10
rho = 28
beta = 8/3

def rk4(function,h,t0,x0):
    f1 = function(t0,x0)
    f2 = function(t0+h/2, x0+(h/2)*f1)
    f3 = function(t0+h/2, x0+(h/2)*f2)
    f4 =  function(t0+h, x0+(h*f3))
    x_res = x0 + (h/6) * ((f1) + (2*f2) + (2*f3) + f4)
    return x_res

def my_lorenz(t,x):
    dx = [sigma*(x[1]-x[0]),x[0]*(rho-x[2])-x[1]
        ,(x[0]*x[1])-(beta*x[2])]
    dX = np.array(dx)
    return dX
# Initial Condition
x0 = np.array([-8, 8, 27])


h = 0.01
t_e = 30 # final time
t = np.arange(0,t_e,h)
# Store the values in X
X = np.zeros((3,len(t)))
X[:,0] = x0
x_in = x0
for i in range (len(t)-1):
    x_res = rk4(my_lorenz,h,t[i],x_in)
    X[:,i+1] = x_res
    x_in = x_res

# Plot

fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection="3d")
fig.set_facecolor('black')
ax.set_facecolor('black')
ax.grid(False)
ax.xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.set_title("Lorenz Simulator")

def update(i):
    ax.view_init(30, -10+ i/2)
    ax.plot3D(X[0,:i], X[1,:i], X[2,:i],'c')
    fig.set_facecolor('black')
    ax.set_facecolor('black')

ani =FuncAnimation(fig=fig,func=update,interval=5)
plt.show()




    
