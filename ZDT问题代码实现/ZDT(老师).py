import random
import math

#0<= X <= 1 , 1<= len(chrom) <=30
def ZDT1(chrom):
    temp_total = 0
    chrom_len = len(chrom)
    fitness1 = chrom[0]
    for i in range(1,chrom_len):
        temp_total = chrom[i] + temp_total
    g = 1 + ( 9 * temp_total ) / ( chrom_len - 1 )
    fitness2 = g * (1 - (fitness1 / g ) ** 0.5)

    return fitness1 , fitness2

#0<= X <= 1 , 1<= len(chrom) <=30
def ZDT2(chrom):
    temp_total = 0
    chrom_len = len(chrom)
    fitness1 = chrom[0]
    for i in range(1, chrom_len):
        temp_total = chrom[i] + temp_total
    g = 1 + (9 * temp_total) / (chrom_len - 1)
    fitness2 = g * (1 - (fitness1 / g) ** 2)

    return fitness1 , fitness2

#0<= X <= 1 , 1<= len(chrom) <=30
def ZDT3(chrom):
    temp_total = 0
    chrom_len = len(chrom)
    fitness1 = chrom[0]
    for i in range(1, chrom_len):
        temp_total = chrom[i] + temp_total
    g = 1 + (9 * temp_total) / (chrom_len - 1)
    print((fitness1 / g) * math.sin(10*math.pi*chrom[0]))
    fitness2 = g * (1 - (fitness1 / g) - (fitness1 / g) * math.sin(10*math.pi*chrom[0]))

    return fitness1 , fitness2

#0<= X1 <= 1 , -5<= X <= 5 , 2<= len(chrom) <=10
def ZDT4(chrom):
    temp_total = 0
    chrom_len = len(chrom)
    fitness1 = chrom[0]
    for i in range(1, chrom_len):
        X2 = math.pow(chrom[i],2) - 10 * math.cos(4 * math.pi * chrom[i])
        print(X2)
        temp_total = X2 + temp_total
    g = 91 + temp_total
    fitness2 = g * (1 - (fitness1 / g) ** 0.5)
    return fitness1 , fitness2

#0<= X <= 1 , 1<= len(chrom) <=10
def ZDT6(chrom):
    x1 = chrom[0]
    temp_total = 0
    chrom_len = len(chrom)
    fitness1 = 1.0 - math.exp((-4.0) * x1) * math.pow(math.sin(6.0 * math.pi * x1), 6.0)
    for i in range(1, chrom_len):
        temp_total = chrom[i] + temp_total
    g = 1 + 9 * math.pow(temp_total/(chrom_len -1),0.25)
    fitness2 = g * (1.0 - math.pow((fitness1/g),2.0))
    return fitness1 , fitness2
# fitness1,fitness2 = ZDT1([0.1,0.2,0.3])
# print(fitness1,fitness2)
