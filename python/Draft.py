# import libraries
from Bio import AlignIO
from Bio.AlignIO import MauveIO
import pprint
import matplotlib.pyplot as plt
import math
import numpy as np




# for i in range(30):
#     n = np.random.randint(1,7)
#     print(n)
t = np.arange(0 , 1.1 , 0.001)
y1 = np.sin(2 * t * 4 * np.pi)
y2 = np.cos(2 * t * 4 * np.pi)
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(t , y1 , color = 'yellow' , lw = 2 , label = 'sin')
plt.annotate("test" , xy = (0.6,1) , xytext = (0.3,1), arrowprops = dict(facecolor = 'black' , shrink = 0.25) ,
                                                                horizontalalignment = 'center' , verticalalignment = 'bottom')
plt.legend(loc = 0)
plt.subplot(1,2,2)
plt.plot(t , y2 , color = 'green' , lw = 4 , label = 'cos')
plt.xlabel('time')
plt.ylabel('value')
plt.title('sin and cos')
plt.legend(loc = 0)
plt.savefig('/home/nehleh/Downloads/sin-cos.pdf' , dpi = 300)
plt.show()

print()




print("test run")

# align = AlignIO.parse("/home/nehleh/0_Research/PhD/Data/100k/GTR_100k", "mauve")
#
# alignments = list(align)
#
#
# for id in range(len(alignments)):
#     # print("MYID************************",id)
#     for idx,record in enumerate(alignments[id]):
#         if record.name == '1':
#             print(record)
#
#
# lcb=arr = [[0 for i in range(5)] for j in range(len(alignments))]
#
#
#
# for id in range(len(alignments)):
#     # print("MYID************************",id)
#     for idx,record in enumerate(alignments[id]):
#             lcb[id][0] = record.name
#             lcb[id][1] = record.annotations['start']
#             lcb[id][2] = record.annotations['end']
#             lcb[id][3] = record.seq
#             lcb[id][4] = record.annotations['strand']
#
#
# pprint.pprint(lcb)
#
# name = []
# for i in range(len(lcb)):
#     name.append(lcb[i][0])
#
# # print(name)
# uniqename= list(set(name))

# for i in uniqename:
#     for id,item in enumerate(lcb):
#         if item[0] == i:
#             print(item)

#
# for record in alignments:
#         print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
#         print(record)