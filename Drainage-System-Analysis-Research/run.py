#!/usr/bin/env python
# coding: utf-8

# In[1]:


import datetime
print(datetime.datetime.now())
import laspy
import numpy as np
import pandas as pd
from numpy.linalg import inv
import math
import datetime
import sys

var1 = int(sys.argv[1])
var2 = int(sys.argv[2])
var3 = int(sys.argv[3])

# In[2]:


infile = laspy.file.File("./3_2.las", mode="r")
ground_x = infile.x
ground_y = infile.y
ground_z = infile.z


# In[3]:


ground_x_2=ground_x-ground_x.min()
ground_y_2=ground_y-ground_y.min()


# In[4]:


len(ground_z)


# In[5]:


x = ground_x_2
y = ground_y_2
z = ground_z


# In[6]:


threedarray = np.vstack((x,y,z)).T
dictionary=pd.Series(threedarray.tolist(), index=map(lambda a: round(a,2), x.tolist()))


# In[7]:


dictionary[1.11]


# In[8]:


list_a=map(lambda a: 10000*int(a), x.tolist())


# In[9]:


list_b=map(lambda a: int(a), y.tolist())


# In[10]:


index_list=[sum(x) for x in zip(list_a, list_b)]


# In[11]:


x=(index_list[0])


# In[12]:


##开出2000*2000的dictionary，存每个里面1*1的值
##如果在里面就在格子内部找，如果在边界就找旁边两个格子，如果在4个交接就找旁边4个
##如果慢的话就做一个新的index=10000*x+y


# In[13]:


grid_dictionary=pd.Series(threedarray.tolist(), index=index_list)


# In[14]:


def find_elevation_new(x, y):
    #bilinear interploation
    storage = []
    curr_index = 10000*int(x)+int(y)
    diff= 0
    while(len(storage)<4):
        ## 找到9宫格中的所有的值
        temp = []
        
        for i in range(int(x-diff), int(x+diff+1)):
            for j in range(int(y-diff), int(y+diff+1)):
                if (i==int(x-diff) or i==int(x+diff) or j==int(y-diff) or j==int(y+diff)):
                    try:
                        rt = grid_dictionary[10000*i+j]
                        if (type(rt)==list):
                            temp.append(rt)
                        else:
                            for it in range(0,len(rt)):
                                temp.append(rt[it])
                    except:
                        useless=1
        temp.sort(key=lambda e: (e[0]-x)*(e[0]-x)+(e[1]-y)*(e[1]-y))
        if (len(storage)+len(temp)<=4):
            ## 全部加进去
            for i in range(0,len(temp)):
                storage.append(temp[i])
        else:
            k=0
            while(len(storage)!=4):
                storage.append(temp[k])
                k=k+1
        diff = diff+1
    new_storage = storage[:4]
    a=[[1,new_storage[0][0],new_storage[0][1],new_storage[0][0]*new_storage[0][1]],
       [1,new_storage[1][0],new_storage[1][1],new_storage[1][0]*new_storage[1][1]],
       [1,new_storage[2][0],new_storage[2][1],new_storage[2][0]*new_storage[2][1]],
       [1,new_storage[3][0],new_storage[3][1],new_storage[3][0]*new_storage[3][1]]]
    b=[new_storage[0][2], new_storage[1][2], new_storage[2][2],new_storage[3][2]]
    try:
        coef_matrix = np.matmul(inv(a), b)
        rt = coef_matrix[0]+coef_matrix[1]*x+coef_matrix[2]*y+coef_matrix[3]*x*y
        return rt
    except:
        return 10000





# In[ ]:



#10 - 100 * 10 - 100 data area
print(datetime.datetime.now())
min_x = var1*100
min_y = var2*100
max_x = var1*100+100
max_y = var2*100+100
increment = 1
angle = var3
x_coord = min_x
y_coord = min_y
rt = []
while x_coord!=(max_x+1):
    while y_coord!=(max_y+1):
        curr_elevation = find_elevation_new(x_coord, y_coord)
        curr_x = x_coord
        curr_y = y_coord
        #print("Starting Point:")
        print((curr_x,curr_y))
        curr_array = []
        while ((curr_x>=min_x and curr_x<=max_x) and (curr_y>=min_y and curr_y<=max_y)):
            #print((curr_x, curr_y, curr_elevation))
            curr_array.append((curr_x, curr_y))
            rt_x = curr_x
            rt_y = curr_y
            rt_elevation = curr_elevation
            angel_diff = 0
            while(angel_diff<360):
                new_x = curr_x+increment*math.sin(math.pi/180*angel_diff)
                new_y = curr_y+increment*math.cos(math.pi/180*angel_diff)
                new_elevation = find_elevation_new(new_x, new_y)
                #print((angel_diff,new_x,new_y,new_elevation,rt_x,rt_y,rt_elevation))
                if (new_elevation<rt_elevation):
                    rt_x = new_x
                    rt_y = new_y
                    rt_elevation = new_elevation
                angel_diff = angel_diff+angle
            if (rt_elevation<curr_elevation):
                curr_x = rt_x
                curr_y = rt_y
                curr_elevation = rt_elevation
            else:
                break
        rt.append(curr_array)
        #print("Result For:")
        #print((x_coord,y_coord))
        #print(curr_array)
        y_coord=y_coord+1
    y_coord = min_y
    x_coord=x_coord+1
print(datetime.datetime.now())



print(rt)

data = []
for i in range(0,101):
    data.append([0]*101)


# In[21]:


for i in range(0,len(rt)):
    for j in range(0,len(rt[i])):
        data[int(math.floor(rt[i][j][0]))-var1*100][int(math.floor(rt[i][j][1]))-var2*100]=data[int(math.floor(rt[i][j][0]))-var1*100][int(math.floor(rt[i][j][1]))-var2*100]+1

print("#################################")

print(data)


import numpy as np
import pandas as pd


filename="outfile_"+str(var1)+"_"+str(var2)+"_"+str(var3)+".csv"

mat = a = np.array(data)
df = pd.DataFrame(data=mat.astype(float))
df.to_csv(filename, sep=' ', header=False, float_format='%.2f', index=False)


print("finished")



