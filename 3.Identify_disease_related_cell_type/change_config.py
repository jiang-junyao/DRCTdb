# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 15:34:07 2023

@author: jjy
"""
### arg1: cluster_use_path
### arg2: original config path
### arg3: modified config output path
import sys
sample_name = []
with open(sys.argv[1],'r')as sample:
    line1 = sample.readlines()
    for i in line1:
        i = i.replace('\n','')
        sample_name.append(i)
    
text2 = '["' + sample_name[0] + '", '
for j  in sample_name[1:len(sample_name)]:
    text2 = text2 + '"' + j + '", '
text2 = text2[0:len(text2)-2]
text2 = text2 + ']'
outlist = []
with open(sys.argv[2],'r')as f1:
    line1 = f1.readlines()
    for i in line1:
        
        if i.startswith('let annotation_names'):
            text1 = i.split('=')[0]
            text1 = text1 + ' = ' +text2
            print(text1)
            outlist.append(text1)
        else:
            outlist.append(i)
            
with open(sys.argv[3],'w+')as f2:
    for i in outlist:
        f2.write(i)