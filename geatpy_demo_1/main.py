"""main.py"""
import geatpy as ea # Import geatpy
from MyProblem import MyProblem # Import MyProblem class

if __name__ == '__main__':
    """=========================Instantiate your problem=========================="""
    M = 3                      # Set the number of objects.
    problem = MyProblem(M)     # Instantiate MyProblem class
    """===============================Population set=============================="""
    Encoding = 'RI'            # Encoding type.
    NIND = 100                 # Set the number of individuals.
    Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # Create the field descriptor.
    population = ea.Population(Encoding, Field, NIND) # Instantiate Population class(Just instantiate, not initialize the population yet.)
    """================================Algorithm set==============================="""
    myAlgorithm = ea.moea_NSGA3_templet(problem, population) # Instantiate a algorithm class.
    myAlgorithm.MAXGEN = 500 # Set the max times of iteration.
    """===============================Start evolution=============================="""
    NDSet = myAlgorithm.run() # Run the algorithm templet.
    """=============================Analyze the result============================="""
    PF = problem.getReferObjV() # Get the global pareto front.
    GD = ea.indicator.GD(NDSet.ObjV, PF) # Calculate GD
    IGD = ea.indicator.IGD(NDSet.ObjV, PF) # Calculate IGD
    HV = ea.indicator.HV(NDSet.ObjV, PF) # Calculate HV
    Space = ea.indicator.spacing(NDSet.ObjV) # Calculate Space
    print('The number of non-dominated result: %s'%(NDSet.sizes))
    print('GD: ',GD)
    print('IGD: ',IGD)
    print('HV: ', HV)
    print('Space: ', Space)