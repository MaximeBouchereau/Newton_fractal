import numpy as np
import matplotlib.pyplot as plt
import sys

Rmin = -1.5                 # Lower bound of the grid where function is studied
Rmax = 1.5                  # Upper bound of the grid where function is studied
N = 10000                   # Number of discretizations
tol = 1e-6                  # Tolerance for the Newton method
I = 100                     # Maximum number of iterations

print(100*"-")
print(" -> 1D Example of polynomial P = X^3 - X + 0.2")
print("  ")
print("  Parameters:")
print("  ")
print("    - Lower bound of the grid where function is studied:",Rmin)
print("    - Upper bound of the grid where function is studied:",Rmax)
print("    - Number of points:",(N+1))
print("    - Tolerance for the Newton method:",tol)
print("    - Maximum number of iterations:",I)
print(100*"-")

def f(x):
    """Function studied as example"""
    return x**3 - x + 0.2

def df(x):
    """Derivative of the function f"""
    return 3*x**2 - 1

def Roots():
    """Gives the roots of the function f"""
    return np.roots([1,0,-1,0.2])

def Newton_Points_1D(rmin = Rmin , rmax = Rmax , ax_equal = True ,  save = False):
    """Newton method over a grid corresponding to the interval [-rmin,rmax] with N+1 points, and determination of the root found
    from the initial point
    Input:
    - rmin: Float - Lower bound of the grid
    - rmax: Float - Upper bound of the grid
    - ax_equal: Boolean - Equal proportions for axes (default: True)
    - save: Boolean - Saved the figure or not (default: False)
    """

    Roots_f = Roots()

    xx = np.linspace(rmin,rmax,N+1) # Grid
    colors = list(np.zeros_like(xx))

    for k in range(len(colors)):
        y = xx[k]
        for n in range(I):
            y = y - f(y)/df(y)
        sys.stdout.write("\r%d   " % int(k) + "/" + str((N + 1)))
        sys.stdout.flush()
        if np.abs(y-Roots_f[0]) < tol:
            colors[k] = "green"
        elif np.abs(y-Roots_f[1]) < tol:
            colors[k] = "red"
        elif np.abs(y-Roots_f[2]) < tol:
            colors[k] = "blue"
        else:
            colors[k] = "black"

    fig = plt.figure()
    ax = fig.add_subplot()
    plt.plot([0,0] , [-0.5,1] , linestyle = "dashed" , color = "black")
    plt.plot(xx,np.zeros_like(xx) , linestyle = "dashed" , color = "black")
    plt.xlim(rmin,rmax)
    plt.ylim(-0.5,1)
    if ax_equal == True:
        ax.set_aspect("equal")
    plt.scatter(xx,f(xx) , color = colors)
    plt.grid()
    plt.savefig("Bassin_Attraction_1D"+str([rmin,rmax])+".pdf" , dpi=(1000) , bbox_inches="tight")
    plt.show()

    pass
