import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

def f(x):
    return x**3 - x + 0.2

Roots_f = scipy.optimize.fsolve(f,[-1.5,0.5,1.5])
tol = 5e-3

xx = np.arange(-1.5,1.5,0.0001)

colors = list(np.zeros_like(xx))

for k in range(len(colors)):
    y = scipy.optimize.newton(func = f, x0 = xx[k], maxiter = 100 , tol = 5e-4)
    if np.abs(y-Roots_f[0]) < tol:
        colors[k] = "green"
    elif np.abs(y-Roots_f[1]) < tol:
        colors[k] = "blue"
    elif np.abs(y-Roots_f[2]) < tol:
        colors[k] = "red"
    else:
        colors[k] = "black"

fig = plt.figure()
ax = fig.add_subplot()
plt.plot([0,0] , [-0.5,1] , linestyle = "dashed" , color = "black")
plt.plot(xx,np.zeros_like(xx) , linestyle = "dashed" , color = "black")
plt.xlim(-1.5,1.5)
plt.ylim(-0.5,1)
ax.set_aspect("equal")
plt.scatter(xx,f(xx) , color = colors)
plt.grid()
plt.savefig("Bassin_Attraction_1D.pdf" , dpi=(1000) , bbox_inches="tight")
plt.show()