"""MyProblem.py"""
import numpy as np
import geatpy as ea

class MyProblem(ea.Problem): # Inherited from Problem class.
    def __init__(self, M): # M is the number of objects.
        name = 'DTLZ1' # Problem's name.
        maxormins = [1] * M # All objects are need to be minimized.
        Dim = M + 4 # Set the dimension of decision variables.
        varTypes = [0] * Dim # Set the types of decision variables. 0 means continuous while 1 means discrete.
        lb = [0] * Dim # The lower bound of each decision variable.
        ub = [1] * Dim # The upper bound of each decision variable.
        lbin = [1] * Dim # Whether the lower boundary is included.
        ubin = [1] * Dim # Whether the upper boundary is included.
        # Call the superclass's constructor to complete the instantiation
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    def aimFunc(self, pop): # Write the aim function here, pop is an object of Population class.
        Vars = pop.Phen # Get the decision variables
        XM = Vars[:,(self.M-1):]
        g = np.array([100 * (self.Dim - self.M + 1 + np.sum(((XM - 0.5)**2 - np.cos(20 * np.pi * (XM - 0.5))), 1))]).T
        ones_metrix = np.ones((Vars.shape[0], 1))
        pop.ObjV = 0.5 * np.fliplr(np.cumprod(np.hstack([ones_metrix, Vars[:,:self.M-1]]), 1)) * np.hstack([ones_metrix, 1 - Vars[:, range(self.M - 2, -1, -1)]]) * np.tile(1 + g, (1, self.M))
    def calReferObjV(self): # Calculate the theoretic global optimal solution here.
        uniformPoint, ans = ea.crtup(self.M, 10000) # create 10000 uniform points.
        realBestObjV = uniformPoint / 2
        return realBestObjV