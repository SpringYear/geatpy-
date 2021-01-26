# -*- coding: utf-8 -*-
"""MyProblem.py"""
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
    '''
    在python中, 如果用一个列表list1乘一个数字n 会得到一个新的列表list2, 这个列表的元素是list1的元素重复n次, 例如
    list1 = [0]
    list2 = list1 * 5		# list2 = [0, 0, 0, 0, 0]
    '''
    
    def aimFunc(self, pop): # 目标函数,ZDT1目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        ObjV1 = Vars[:, 0]
        gx = 1 + 9 * np.sum(Vars[:, 1:30], 1)
        hx = 1 - np.sqrt(ObjV1 / gx)
        ObjV2 = gx * hx

        #下面这句话的含义不懂
        pop.ObjV = np.array([ObjV1, ObjV2]).T # 把结果赋值给ObjV
    '''
    X[:,0]是numpy中数组的一种写法，表示对一个二维数组，取该二维数组第一维中的所有数据，第二维中取第0个数据，
    直观来说，X[:,0]就是取所有行的第0个数据, X[:,1] 就是取所有行的第1个数据。

    举例说明：
    import numpy as np
    X = np.array([[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15],[16,17],[18,19]])
    print X[:,0]

    结果：
    [0,2,4,6,8,10,12,14,16,18]

    X[n,:]是取第1维中下标为n的元素的所有值。

    X[1,:]即取第一维中下标为1的元素的所有值

    X[:,  m:n]，即取所有数据的第m到n-1列数据，含左不含右


    [::-1]  取反
    a = "abcdef"
    print (a[::-1])
    结果：fedcba

    [:-1] 删除最后一位
    a = "abcdef"
    print (a[:-1])
    结果：abcde
    '''
    
    def calReferObjV(self): # 计算全局最优解
        N = 10000 # 生成10000个参考点
        ObjV1 = np.linspace(0, 1, N)
        ObjV2 = 1 - np.sqrt(ObjV1)
        globalBestObjV = np.array([ObjV1, ObjV2]).T
        return globalBestObjV

    '''
    np.linspace(start, stop, num=50, endpoint=True, retstep=False, dtype=None)
    在规定的时间内，返回固定间隔的数据。他将返回“num”个等间距的样本，在区间[start, stop]中。
    其中，区间的结束端点可以被排除在外。

    属性	说明
    start	队列的开始值
    stop	队列的结束值
    num	要生成的样本数，非负数，默认是50
    endpoint	若为True，“stop”是最后的样本；否则“stop”将不会被包含。默认为True
    retstep	若为False，返回等差数列；否则返回array([samples, step])。默认为False
    '''