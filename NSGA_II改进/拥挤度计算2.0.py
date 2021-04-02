#-*- coding: utf-8 -*-
import random
import math

#生成随机数组
def rando():
    a = random.randint(0,10)/10
    b = random.randint(0,10)/10
    c = random.randint(0,10)/10
    return [a,b,c]

#0<= X <= 1 , 1<= len(chrom) <=30


#ZDT问题，得出结果F1，F2
def ZDT1(chrom):
    temp_total = 0
    chrom_len = len(chrom)
    fitness1 = chrom[0]
    for i in range(1,chrom_len):
        temp_total = chrom[i] + temp_total
    g = 1 + ( 9 * temp_total ) / ( chrom_len - 1 )
    fitness2 = g * (1 - (fitness1 / g ) ** 0.5)
    return fitness1 , fitness2


F1=[]#在循环之前先创建好列表
for i in range(0,100):
    c = rando()
    fitness1,fitness2 = ZDT1(c)
    #F1.append(fitness1)     #列表增加元素
    F1.append([fitness1,fitness2])

Pdistance1 = 0

#计算距离临界值q
q = ((F1[99][1]-F1[0][1])**2+(F1[99][0]-F1[0][0])**2)**0.5/(2*100)

#淘汰个体
for i in range (1,98):
    h = ((F1[i+1][1]-F1[i][1])**2+(F1[i+1][0]-F1[i][0])**2)**0.5  #计算两个体的欧氏距离
    if h <= q:
        #比较 a b 距离中心点的距离
        #中心点e的坐标abs((F1[i-1][0]+F1[i+2][0])/2)  ,  abs((F1[i-1][1]+F1[i+2][1])/2)
        x=abs((F1[i-1][0]+F1[i+2][0])/2)
        y=abs((F1[i-1][1]+F1[i+2][1])/2)
        a1=(F1[i][0]-x)**2+(F1[i][1]-y)**2
        b1=(F1[i+1][0]-x)**2+(F1[i+1][1]-y)**2
        if a1<=b1:
            del(F1[i+1])
        else:
            del(F1[i]) 
print(F1)    #测试步骤
n=len(F1)


for i in range (1,n-1):
    Pdistance = abs(F1[i+1][0]-F1[i-1][0])+abs(F1[i+1][1]-F1[i-1][1])
    Pdistance1 = Pdistance1 + Pdistance
print(Pdistance1)