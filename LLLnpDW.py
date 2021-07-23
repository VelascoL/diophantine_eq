# -*- coding: utf-8 -*-
"""
Solves the family of Diophantine equations of the form
u_{n_1} + ... + u_{n_k} = wp_1^{z_1} ... p_s^{z_s},
under some mild technical restrictions.
@authors: Brian Ha, Lily McBeath, Luisa Velasco
"""
import math
import numpy as np
from functools import reduce
from sage.all import *
#from sage.rings.number_field.S_unit_solver import minimal_vector
from constants1 import Constants
from fractions import Fraction
from decimal import Decimal
#import mpmath as mp
from fpylll import *
import operator


def generate_lattice(large_constant, eta_list):
    size = len(eta_list)
    lattice = np.identity(size-1)
    lattice = np.concatenate((lattice, np.zeros((size-1,1))),axis=1)
    row = []
    for i in range(len(eta_list)):
        row.append(math.floor(Fraction(large_constant) * Fraction(math.log(eta_list[i]))))
    print(row)
    row = np.reshape(row,(1,size))
    lattice = np.concatenate((lattice,row),axis=0)
    return lattice
    
def initial(basis):
    D = []
    D.append(1)
    casis = np.zeros((np.shape(basis)))
    casis = Matrix(ZZ,casis)
    Lambda = np.zeros((np.shape(basis)))
    Lambda = Matrix(ZZ,Lambda)
    for i in range(0,np.shape(basis)[0]):
        casis[:,i] = basis[:,i]
        for j in range(0,i):
            print(i,j,basis.column(i), casis.column(j))
            Lambda[i,j] = basis.column(i).dot_product(casis.column(j))
            casis[:,i] = (D[j+1]*casis.column(i) - Lambda[i][j]*casis.column(j))/D[j]
        D.append((casis.column(i).dot_product(casis.column(i)))/D[i])
    D = vector(ZZ,D)
    Lambda = Matrix(ZZ,Lambda)
    return Lambda,D

def prod_A(k,l,basis,Lambda,D):
    if 2*abs(Lambda[k-1,l-1]) > D[l]:
        r = (Lambda[k-1,l-1]/D[l]).round()
        print('r',r)
        basis[:,k-1] = basis[:,k-1] - r*basis[:,l-1]
        for j in range(0,l-1):
            Lambda[k-1,j] = Lambda[k-1,j] - r*Lambda[l-1,j]
        Lambda[k-1,l-1] = Lambda[k-1,l-1] - r*D[l]
    return k,l,basis,Lambda,D
        
def prod_B(k,basis,Lambda,D):
    basis[:,[k-2,k-1]] = basis[:,[k-1,k-2]]
    for j in range(0,k-2):
        Lambda[[k-2,k-1],j] = Lambda[[k-1,k-2],j]
    for i in range(k,np.shape(basis)[0]):
        t = Lambda[i,k-2]
        Lambda[i,k-2] = (Lambda[i,k-2]*Lambda[k-1,k-2] + Lambda[i,k-1]*D[k-2])/D[k-1]
        Lambda[i,k-1] = (t*D[k] - Lambda[i,k-1]*Lambda[k-1,k-2])/D[k-1]
    D[k-1] = (D[k-2]*D[k] + Lambda[k-1,k-2]**2)/D[k-1]
    return k,basis,Lambda,D
    
def LLLDW(basis):
    Lambda,D = initial(basis)
    print('basis',basis)
    print('L',Lambda)
    print('D', D)
    k=2
    while k <= np.shape(basis)[0]: 
        l = k-1
        print('doing prod A1')
        k,l,basis,Lambda,D = prod_A(k,l,basis,Lambda,D)
        print(Lambda)
        if 4*D[k-2]*D[k] < (3*D[k-1]**2 - 4*Lambda[k-1,k-2]**2):
            print('doing prod B',k)
            k,basis,Lambda,D = prod_B(k,basis,Lambda,D)
            print(Lambda,D)
            if k > 2:
                k = k-1
        else:
            print('doing the other thing',k)
            for l in range(k-2,2):
                print('doing prod A2')
                k,l,basis,Lambda,D = prod_A(k,l,basis,Lambda,D)
                print(Lambda)
            k = k + 1
    print(basis)
    print(Lambda)
    print(D)
    return basis

def calculate_norm_squared(v):
    return v.dot_product(v)

def calculate_max_frac_norm_squared(lll_reduced_basis, gs_basis):
    return max([calculate_norm_squared(lll_reduced_basis.column(0))/calculate_norm_squared(gs_basis.column(i)) for i in range(size)])

def c1(red):
    red = red.transpose()
    gs,mu = red.gram_schmidt()
    print(gs)
    c1 = calculate_max_frac_norm_squared(red,gs)
    return c1

def c4(c1,basis):
    b1 = basis.column(0)
    b1sq = calculate_norm_squared(b1)
    c4 = (c1**(-1))*b1sq
    return c4


    
def generate_approximation_matrix(large_constant):
        n = len(primes)
        primes_row = [round(large_constant * log(p)) for p in primes]
        approximation_matrix = []
        for i in range(n):
            zero_row = [0] * n
            zero_row[i] = 1
            approximation_matrix.append(zero_row)
        approximation_matrix[n - 1] = primes_row
        return Matrix(ZZ, approximation_matrix)
        
    
if __name__ == "__main__":
    large_constant = 10**(10)
    primes = [2,3,5]
    size = len(primes)
    basis = Matrix(ZZ,generate_approximation_matrix(large_constant))
    red = LLLDW(basis)
    print(red)
    
    


