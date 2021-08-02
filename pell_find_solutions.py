# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:51:05 2021

@author: luisa
"""

import numpy as np 

def pell(F0,F1,n):
    l = [F0,F1]
    for i in range(2,n):
        new = 2*l[i-1] + l[i-2]
        #print(new)
        l.append(new)
    return l

flist = pell(0,1,217)

print(flist)

for n in range(1,217):
    for m in range(1,n+1):
        for l in range(1,n+1):
            for k in range(1,n+1):
                for a in range(n+1):
                    for b in range(n+1):
                        for c in range(n+1):
                            if flist[n] + flist[m] + flist[l] + flist[k] - (2**a)*(3**b)*(5**c) == 0:
                                print(n,m,l,k,a,b,c)
            

