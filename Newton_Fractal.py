import numpy as np
import math as mt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import sys

# Newton fractal

# Parameters

R = 2                       # Size of the square where function is studied
N = 873                     # Number of discretizations
tol = 1e-6                  # Tolerance for the Newton method
I = 50                      # Maximum number of iterations
pol = [1,0,0,0,0,-1]        # Polynomial which is studied as function (List of coeffcients from highest to lowst degree)


# Class for polynomial print

class PrintPolynomial:
    """Prints polynomials"""

    def string_mon(self,a, i):
        """Returns a character string in order to print the monomial function az^i
        Inputs:
        - a: Float - Coefficient
        - i: Int - Degree of the monomial function"""
        string_coeff = ""
        string_deg = ""
        epsilon = 1e-12
        if np.abs(a - 1) < epsilon and i > 1:
            # print(1)
            string_coeff, string_deg = " + ", "z^" + str(i)
        elif np.abs(a + 1) < epsilon and i > 1:
            # print(2)
            string_coeff, string_deg = " - ", "z^" + str(i)
        elif np.abs(a - 1) < epsilon and i == 1:
            # print(3)
            string_coeff, string_deg = " + ", "z"
        elif np.abs(a + 1) < epsilon and i == 1:
            # print(4)
            string_coeff, string_deg = " - ", "z"
        elif np.abs(a - 1) < epsilon and i == 0:
            # print(5)
            string_coeff, string_deg = " + 1", ""
        elif np.abs(a + 1) < epsilon and i == 0:
            # print(6)
            string_coeff, string_deg = " - 1", ""
        elif np.abs(a) < epsilon:
            string_coeff, string_deg = "", ""
        else:
            # print(7)
            if i == 0:
                if a < 0:
                    string_coeff, string_deg = " - " + str(np.abs(a)), ""
                if a > 0:
                    string_coeff, string_deg = " + " + str(np.abs(a)), ""

            if i == 1:
                if a < 0:
                    string_coeff, string_deg = " - " + str(np.abs(a)), "z"
                if a > 0:
                    string_coeff, string_deg = " + " + str(np.abs(a)), "z"

            if i > 1:
                if a < 0:
                    string_coeff, string_deg = " - " + str(np.abs(a)), "z^" + str(i)
                if a > 0:
                    string_coeff, string_deg = " + " + str(np.abs(a)), "z^" + str(i)

        return string_coeff + string_deg

    def string(self,pol):
        """Returns a character string in order to print the polynomial whose coefficients are given by pol"""
        string = ""

        for i in range(len(pol)):
            string = string + self.string_mon(pol[i], len(pol) - i - 1)
        return string[2:]

# Print of parameters

print(100*"-")
print("  Parameters:")
print(" ")
print("    - Size of the square where function is studied:",R)
print("    - Number of points:",(N+1)**2)
print("    - Tolerance for the Newton method:",tol)
print("    - Maximum number of iterations:",I)
#print("    - Polynomial studied:","f(z) = "+''.join([" + ["+str(pol[i])+"] z^"+str(len(pol)-i-1) for i in range(len(pol))]))
print("    - Polynomial studied: f(z) =", PrintPolynomial().string(pol))
print(100*"-")

# Functions

def degree(pol):
    """Gives the degree of the polynomial (and the number of roots with multiplicity"""
    return len(pol)-1

def f(z):
    """Function whom roots are studied"""
    return sum([pol[i]*z**(len(pol)-i-1) for i in range(len(pol))])

def df(z):
    """Derivative of the studied function"""
    h = 1e-6
    return (f(z+h) - f(z-h))/(2*h)

def Roots_f():
    """Returns roots of the studied function"""
    return np.roots(pol[::-1])

def Points():
    """Computes points"""

    P =  - np.ones((N+1,N+1))

    Rf = Roots_f()
    LR = len(Rf)

    for i in range(N+1):
        for j in range(N+1):
            sys.stdout.write("\r%d   "% int(i*(N+1) + j + 1) + "/"+str((N+1)**2))
            sys.stdout.flush()
            z = - R + 2*R*i/(N) + (R - 2*R*j/(N))*(1j)
            It = 0

            while min([np.abs(z - Rf[k]) for k in range(LR)]) > tol:
                if np.abs(df(z))< 1e-12:
                    z = z - f(z)/(df(z) + 1e-11)
                else:
                    z = z - f(z)/df(z)

                It = It + 1

                if It > I:
                    break

            for k in range(LR):
                if np.abs(z-Rf[k]) < tol:
                    P[j,i] = k + It/I

    np.save("Points_"+str((N+1)**2)+"_Iter_"+str(I)+"_degree_"+str(len(pol)-1),P)
    pass

def Fractal(name_points,save=False):
    """Plot the Newton fractal"""
    POINTS = np.load(name_points)

    color_Roots = ["red", "green", "blue", "orange", "purple"]

    deg = degree(pol)
    List_Colors = [(0, "black") , (1/(deg+1),"black")]
    for k in range(deg):
        List_Colors = List_Colors + [((k+1)/(deg+1),color_Roots[k]),((k+2)/(deg+1),"white")]
        Cmap = mpl.colors.LinearSegmentedColormap.from_list("", List_Colors)
    plt.imshow(POINTS,cmap=Cmap,vmin=-1,vmax=len(Roots_f()),extent=[-R,R,-R,R],interpolation="bilinear")

    Real_Roots , Im_Roots = [z.real for z in Roots_f()] , [z.imag for z in Roots_f()]
    plt.scatter(Real_Roots , Im_Roots , color = color_Roots[0:len(Roots_f())] , edgecolors="black")
    plt.title("f(z) = "+PrintPolynomial().string(pol))
    if save == True:
        plt.savefig("Newton_Fractal_Points_"+str((N+1)**2)+"_Iter_"+str(I)+"_"+PrintPolynomial().string(pol)+".pdf",dpi=(1000),bbox_inches="tight")
    plt.show()
    pass






