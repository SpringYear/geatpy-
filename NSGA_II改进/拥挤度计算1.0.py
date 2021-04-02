#-*- coding: utf-8 -*-
import random
import math



'''
随机生成数
给出0到10之间的随机整数：
随机数包括边界，及包括0和10
import random
a = random.randint(0,10)

从9、19、29、39、……、99之间，随机选取一个实数：
a = random.randrange(9, 100, 10)

从列表[5,6,7,8,9]里面，随机选取一个数：
a = random.choice([5,6,7,8,9])

从一个字符串里面，随机选取一个字符：
a = random.choice("从一个字符串里面，随机选取一个字符！")

随机打乱列表里面的字符顺序：
a = ["p","q","r","s","t","p","q","r","s","t","p","q","r","s","t",]
random.shuffle(a)

从列表里面随机选取9个数字：
a = range(3,100,2)
b = random.sample(a, 9)

'''

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

'''
def ini():
    list = [[0,0,0],]
    for i in range(0,100):
        a = random.randint(0,10)/10
        b = random.randint(0, 10) / 10
        c = random.randint(0, 10) / 10
    return list


def sel():
    list = [[0,0],]
    for i in range(0,100):
        fitness1, fitness2 = ZDT1(list[i])
    return list
'''


F1=[]#在循环之前先创建好列表
F2=[]
F3=[]
for i in range(0,100):
    c = rando()
    fitness1,fitness2 = ZDT1(c)
    #F1.append(fitness1)     #列表增加元素
    #F2.append(fitness2)
    F3.append([fitness1,fitness2])
#Pdistance=abs(F1[0]-F1[3])+abs(F2[0]-F2[3])
#print(Pdistance)
Pdistance1 = 0
#开始排序  冒泡排序
for b in range(0,100):  # b从0到99
    for c in range(0,99):
        if  F3[c][0]<F3[c+1][0]:
            d = F3[c]
            F3[c] = F3[c+1]
            F3[c+1] = d
#计算距离临界值q
q = ((F3[99][1]-F3[0][1])**2+(F3[99][0]-F3[0][0])**2)**0.5/(2*100)
#淘汰个体
for i in range (1,98):
    print(i)
    h = ((F3[i+1][1]-F3[i][1])**2+(F3[i+1][0]-F3[i][0])**2)**0.5
    if h <= q:
        #比较 a b 距离中心点的距离
        #中心点e的坐标abs((F3[i-1][0]+F3[i+2][0])/2)  ,  abs((F3[i-1][1]+F3[i+2][1])/2)
        x=abs((F3[i-1][0]+F3[i+2][0])/2)
        y=abs((F3[i-1][1]+F3[i+2][1])/2)
        a1=(F3[i][0]-x)**2+(F3[i][1]-y)**2
        b1=(F3[i+1][0]-x)**2+(F3[i+1][1]-y)**2
        if a1<=b1:
            for j in range(i+1,98):
                F3[i+1][0]=F3[i+2][0]
                F3[i+1][1]=F3[i+2][1]
            del(F3[-1])
        else:
            for i in range(i,99):
                F3[i][0]=F3[i+1][0]
                F3[i][1]=F3[i+1][1]
            del(F3[-1]) 
n=len(F3)



for i in range (1,n-1):
    Pdistance = abs(F3[i+1][0]-F3[i-1][0])+abs(F3[i+1][1]-F3[i-1][1])
    Pdistance1 = Pdistance1 + Pdistance
#fitness3,fitness4 = ZDT1([0.3,0.4,0.5])
#Pdistance=abs(fitness3-fitness1)+abs(fitness4-fitness2)

print(Pdistance1)


"""
    def crowding_distance_assignment( I )
            nLen = len( I )        #I中的个体数量
        for i in I:
            i.distance = 0    #初始化所有个体的拥挤距离
        for objFun in M:        #M为所有目标函数的列表
                    I = sort( I, objFun )    #按照目标函数objFun进行升序排序
                    I[0] = I[ len[I]-1 ] = ∞    #对第一个和最后一个个体的距离设为无穷大
                    for i in xrange( 1, len(I) - 2 ):
                            I[i].distance = I[i].distance + ( objFun( I[i+1] ) - objFun( I[i-1] ) )/(Max(objFun()) - Min(objFun()) )
        
        一.sort函数

        1.sort函数包含在头文件为#include<algorithm>的c++标准库中，调用标准库里的排序方法可以实现对数据的排序，但是sort函数是如何实现的，我们不用考虑！

        2.sort函数的模板有三个参数：

        void sort (RandomAccessIterator first, RandomAccessIterator last, Compare comp);
        （1）第一个参数first：是要排序的数组的起始地址。

        （2）第二个参数last：是结束的地址（最后一个数据的后一个数据的地址）

        （3）第三个参数comp是排序的方法：可以是从升序也可是降序。如果第三个参数不写，则默认的排序方法是从小到大排序。




        Python中求数字的平方根和平方的几种方法

        方法一: 使用内置模块
        >>> import math
        >>> math.pow(12, 2)     # 求平方
        144.0
        
        >>> math.sqrt(144)      # 求平方根
        12.0
        
        方法二: 使用表达式
        >>> 12 ** 2             # 求平方
        144
        
        >>> 144 ** 0.5          # 求平方根
        12.0
        
        方法三: 使用内置函数
        >>> pow(12, 2)          # 求平方
        144
        
        >>> pow(144, .5)        # 求平方根
        12.0


        python中列表去掉最后一个元素
        在Python3中列表数据类型的内置方法里有三种方法可以删除列表的最后一个元素

        pop方法

        list = [1,2,3,4]
        list.pop()
        print(list)
        
        >>>[1, 2, 3]

        del方法

        list = [1,2,3,4]
        del(list[-1])
        print(list)
        
        >>>[1, 2, 3]

        切片

        list = [1,2,3,4]
        list = list[0:-1]
        print(list)

        >>>[1,2,3]
        
        总结：以上三种方法未在内存处理上进行测试，
        唯一区别，pop方法和del方法如果对空列表进行操作，会报错中断执行，
        切片方法不会因此报错，继续保持空列表向下运行

        if-elif-else  结构

        经常需要检查超过两个的情形，可以使用if-elif-else 结构，
        它依次检查每个条件测试，直到遇到通过了的条件测试.如果前一个测试通过了，
        就不进行下一个测试
"""