#-*- coding: utf-8 -*-
import random
import math
import numpy as np    # NumPy:数学基础库
import geatpy as ea

class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'ZDT1' # 初始化name（函数名称，可以随意设置）
        M = 2 # 初始化M（目标维数）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 30 # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [0] * Dim # 决策变量下界
        ub = [1] * Dim # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界
        ubin = [1] * Dim # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)  #考虑一下是否是递归函数

        #类中完成改进ZDT问题

    #生成随机数组
    def rando(self):
        a = random.randint(0,10)/10
        b = random.randint(0,10)/10
        c = random.randint(0,10)/10
        return [a,b,c]

    #ZDT问题，得出结果F1，F2
    def ZDT1(self,chrom):
        temp_total = 0
        chrom_len = len(chrom)
        fitness1 = chrom[0]
        for i in range(1,chrom_len):
            temp_total = chrom[i] + temp_total
        g = 1 + ( 9 * temp_total ) / ( chrom_len - 1 )
        fitness2 = g * (1 - (fitness1 / g ) ** 0.5)
        return fitness1 , fitness2

    def mainZDT(self):
        F1=[]#在循环之前先创建好列表
        for i in range(0,100):
            c = self.rando()
            fitness1,fitness2 = self.ZDT1(c)
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
        return F1
    


    def aimFunc(self, pop): # 目标函数,ZDT1目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        F1 = self.mainZDT()
        x = np.array(F1)
        ObjV1 = x[:, 0]
        ObjV2 = x[:, 1]
        #下面这句话的含义不懂
        pop.ObjV = np.array([ObjV1, ObjV2]).T # 把结果赋值给ObjV

    def calReferObjV(self): # 计算全局最优解
        N = 10000 # 生成10000个参考点
        ObjV1 = np.linspace(0, 1, N)
        ObjV2 = 1 - np.sqrt(ObjV1)
        globalBestObjV = np.array([ObjV1, ObjV2]).T
        return globalBestObjV