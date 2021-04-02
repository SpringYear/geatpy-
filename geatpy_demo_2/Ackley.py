import numpy as np
import geatpy as ea

class Ackley(ea.Problem): # Inherited from Problem class.
    def __init__(self, D = 30):
        name = 'Ackley' # Problem's name.
        M = 1 # Set the number of objects.
        maxormins = [1] * M # All objects are need to be minimized.
        Dim = D # Set the dimension of decision variables.
        varTypes = [0] * Dim # Set the types of decision variables. 0 means continuous while 1 means discrete.
        lb = [-32.768] * Dim # The lower bound of each decision variable.
        ub = [32.768] * Dim # The upper bound of each decision variable.
        lbin = [1] * Dim # Whether the lower boundary is included.
        ubin = [1] * Dim # Whether the upper boundary is included.
        # Call the superclass's constructor to complete the instantiation
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    def aimFunc(self, pop): # Write the aim function here, pop is an object of Population class.
        x = pop.Phen # Get the decision variables
        n = self.Dim
        f = np.array([-20 * np.exp(-0.2*np.sqrt(1/n*np.sum(x**2, 1))) - np.exp(1/n * np.sum(np.cos(2 * np.pi * x), 1)) + np.e + 20]).T
        return f, CV
    def calReferObjV(self): # Calculate the global optimal solution here.
        realBestObjV = np.array([[0]])
        return realBestObjV