import random
import math

#0<= X <= 1 , 1<= len(parameter) <=30
def ZDT1(chrom):
    function1=chrom[0]
    quantal=0
    length=len(chrom)
    for i in range(1,length):
        quantal = quantal+chrom[i]
    g=(9*quantal/(length-1))+1
    h=1-((chrom[0]/g)**0.5)
    function2=g*h
    return function1,function2
    
#0<= X <= 1 , 1<= len(parameter) <=30
def ZDT2(chrom):
    function1=chrom[0]
    quantal=0
    length=len(chrom)
    for i in range(1,length):
        quantal = quantal+chrom[i]
    g=(9*quantal/(length-1))+1
    h=1-((chrom[0]/g)**2)
    function2=g*h
    return function1,function2
 
#0<= X <= 1 , 1<= len(parameter) <=30 
def ZDT3(chrom):
    function1=chrom[0]
    quantal=0
    length=len(chrom)
    for i in range(1,length):
        quantal = quantal+chrom[i]
    g=(9*quantal/(length-1))+1
    h=1-((chrom[0]/g)**0.5)-((chrom[0]/g)*(math.sin(10*math.pi*chrom[0])))
    function2=g*h
    return function1,function2
#0<= X1 <= 1 , -5<= X <= 5 , 2<= len(chrom) <=10
def ZDT4(chrom):
    function1=chrom[0]
    quantal=0
    length=len(chrom)
    for i in range(1,length):
        Summary=chrom[i]**2-10*math.cos*(4*math.pi*chrom[i])
        quantal = quantal+Summary
    g=(91+quantal/(length-1))+1
    h=1-((chrom[0]/g)**0.5)
    function2=g*h
    return function1,function2

#0<= X <= 1 , 1<= len(chrom) <=10
def ZDT6(chrom):
    function1=1-math.exp((-4)*chrom[0])*(math.sin(6*math.pi*chrom[0])**6)
    quantal=0
    length=len(chrom)
    for i in range(1,length):
        quantal = quantal+chrom[i]
    g=(quantal/(length-1))**0.25*(length-1)+1#注意length是变化还是固定的9
    h=1-((function1/g)**2)
    function2=g*h
    return function1,function2
    
function1,function2=ZDT1([0.1,0.2,0.3])
print(function1,function2)
