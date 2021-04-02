import geatpy as ea # import geatpy
import numpy as np
from MyProblem import Ackley

if __name__ == '__main__':
    """=========================Instantiate your problem=========================="""
    problem = Ackley(30) # Instantiate MyProblem class.
    """===============================Population set=============================="""
    Encoding = 'RI'                # Encoding type.
    NIND = 20                      # Set the number of individuals.
    Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # Create the field descriptor.
    population = ea.Population(Encoding, Field, NIND) # Instantiate Population class(Just instantiate, not initialize the population yet.)
    """================================Algorithm set==============================="""
    myAlgorithm = ea.soea_DE_rand_1_bin_templet(problem, population) # Instantiate a algorithm class.
    myAlgorithm.MAXGEN = 1000      # Set the max times of iteration.
    myAlgorithm.mutOper.F = 0.5    # Set the F of DE
    myAlgorithm.recOper.XOVR = 0.2 # Set the Cr of DE (Here it is marked as XOVR)
    myAlgorithm.drawing = 1 # 1 means draw the figure of the result
    """===============================Start evolution=============================="""
    [population, obj_trace, var_trace] = myAlgorithm.run() # Run the algorithm templet.
    """=============================Analyze the result============================="""
    best_gen = np.argmin(obj_trace[:, 1]) # Get the best generation.
    best_ObjV = np.min(obj_trace[:, 1])
    print('The objective value of the best solution is: %s'%(best_ObjV))
    print('Effective iteration times: %s'%(obj_trace.shape[0]))
    print('The best generation is: %s'%(best_gen + 1))
    print('The number of evolution is: %s'%(myAlgorithm.evalsNum))