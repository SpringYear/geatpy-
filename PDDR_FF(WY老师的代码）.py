import random
from matplotlib import pyplot as plt

def PDDR_FF(fitness1,fitness2):
    Q = []
    P = []
    q = 0
    p = 0
    pddr_ff = []
    for d in range(len(fitness1)):
        for i in range(len(fitness1)):
            if i == d:
                continue
            else:
                if fitness1[d] <= fitness1[i]:
                    if fitness2[d] <= fitness2[i]:
                        p += 1
                elif fitness1[d] >= fitness1[i]:
                    if fitness2[d] >= fitness2[i]:
                        q += 1
        P.append(p)         #append() 方法用于在列表末尾添加新的对象。
        Q.append(q)
        pddr_ff.append(q + 1 / (1 + p))
        p = 0
        q = 0
    return Q,P,pddr_ff



# a = [None] * 10
# b = [None] * 10
# for i in range(len(a)):
#     a[i] = random.random()
#     b[i] = random.random()
#
# a1 = [None] * 10
# b1 = [None] * 10
# for i in range(len(a1)):
#     a1[i] = random.random()
#     b1[i] = random.random()
#
#
# # plt.title("Matplotlib demo")
# # plt.xlabel("x axis caption")
# # plt.ylabel("y axis caption")
# # plt.scatter(a,b)
# # plt.show()
#
# # plt.subplot(221, facecolor='r')  #facecolor指定背景颜色 在之前的Python版本使用的是axisbg 现在已经改成了facecolor
#
# plt.figure(figsize=(8, 5), dpi=80)
# ax = plt.subplot(221)#大图套小图
# p1 = ax.scatter(a,b,marker = '*',color = 'r',label='1',s=10)
# p2 = ax.scatter(a1,b1,marker = '+',color = 'g',label='2',s=20)
#
# ax1 = plt.subplot(222)#大图套小图
# p3 = ax1.scatter(a,b,marker = '*',color = 'r',label='1',s=10)
# p4 = ax1.scatter(a1,b1,marker = '+',color = 'g',label='2',s=20)#s大小
#
# ax.legend((p1, p2), (u'不喜欢', u'魅力一般', u'极具魅力'), loc=2)
# ax1.legend((p3, p4), (u'不喜欢', u'魅力一般', u'极具魅力'), loc=2)
# # plt.show()
#
#
# qq,ww ,ee= PDDR_FF(a,b)
# print(qq)
# print(ww)
# print(ee)
# plt.show()