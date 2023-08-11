# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 16:43:39 2021

@author: hanzy
"""

import numpy as np # 导入packages
import time

# 读取制作好的．fq数据集文件
with open("header_1.fq", 'rt') as f: #根据header_1.fq所在位置对应修改路径，默认数据集和代码文档在同一路径下
    text = f.readlines()
    qual_list = []
    for i in range(len(text)//4):
        qual_list.append(text[4*i+3])
    base_list = []
    for j in range(len(text)//4):
        base_list.append(text[4*j+1][:150]) #对于华大的数据集，第151列是无效数据，所以此处只取前150
        #若是illumina的数据集，其列数为101/103交替， 只取前100列

# Quality Scores 编码为uint
ascii_list = list()

for i in range(len(qual_list)):
    trial = list(str.encode(qual_list[i]))
    trial_np = np.array(trial, dtype=np.uint8) # Quality Scores 编码为uint
    ascii_list.append(trial_np)
    
    
ascii_list_100 = []

#截取前150列
for i in range(len(ascii_list)):
    ascii_list_100.append(ascii_list[i][:150])  ##########
    
# 转化为矩阵
matrix = np.array(ascii_list_100)

#查看matrix的大小
matrix.shape


# 截取前60000
matrix = matrix[:60000,:]
base_list = base_list[:60000]


# 对数据进行分段化处理
def quan(q):
    if q < 60:
        q = 60
    elif 60 <= q < 62:
        q = 62
    elif 62 <= q < 64:
        q = 64
    elif 64 <= q < 66:
        q = 66
    elif 66 <= q < 68:
        q = 68
    elif 68 <= q < 70:
        q = 70
    else:
        q = q
    return q

# 概率表
row_mean = matrix.mean(axis = 1)
total_entropy = 0
freq_dic = {}
starttime = time.time()

for k in range(len(row_mean)):
    row_mean[k] = quan(row_mean[k])

for j in range(4,matrix.shape[1]):#150
    for i in range(2, matrix.shape[0]):#10000
        #context_list = (max(matrix[i,j-4],matrix[i,j-3]), max(matrix[i,j-2],matrix[i,j-1]), min(matrix[i-2,j],matrix[i-1,j]),row_mean[i])
        context_list = (#max(matrix[i,j-4],matrix[i,j-3]),
                        max(matrix[i,j-2],matrix[i,j-1]), 
                         row_mean[i])    #base_list[i][j], base_list[i][j-1],
        #context_list = (max(matrix[i,j-2],matrix[i,j-1]),row_mean[i])
        if context_list not in freq_dic.keys():
            freq_dic[context_list] = {matrix[i,j]:1}
            total_entropy += 8
        else:
            if matrix[i,j] not in freq_dic[context_list].keys():
                freq_dic[context_list].update({matrix[i,j]:1})
                total_entropy += 8
            else:
                total_entropy += -np.log2(freq_dic[context_list][matrix[i,j]]/np.sum(list(freq_dic[context_list].values())))
                freq_dic[context_list][matrix[i,j]] += 1

endtime = time.time()
print("total_time: ",(endtime - starttime))

print('total_entropy：',total_entropy)
