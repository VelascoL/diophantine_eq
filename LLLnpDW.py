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


def generate_lattice_approx(large_constant,primes):
    """
    Generates the lattice approximation matrix not using fyplll
    """
    size = len(primes)
    
    approximation_matrix = np.identity(size-1)
    # Append an empty column to the identity matrix
    approximation_matrix = np.concatenate((approximation_matrix,np.zeros((size-1,1))),axis = 1)
    # Create the row that represents the approximation of the linear form
    row_approx = []
    for i in range(len(primes)):
        #if i == len(primes)-1:
         #   row_approx.append(-math.floor(Fraction(large_constant) * Fraction(math.log(primes[i]))))
        #else:
        row_approx.append(np.rint(Fraction(large_constant) * Fraction(math.log(primes[i]))))
    row_approx = np.reshape(row_approx,(1,size))
    print(np.shape(approximation_matrix))
    print(np.shape(row_approx), np.shape(approximation_matrix))
    approximation_matrix = np.concatenate((approximation_matrix, row_approx),axis = 0)
    return approximation_matrix

def generate_lattice(large_constant, eta_list):
    size = len(eta_list)
    lattice = np.identity(size-1)
    lattice = np.concatenate((lattice, np.zeros((size-1,1))),axis=1)
    row = []
    for i in range(len(primes)):
        row.append(math.floor(Fraction(large_constant) * Fraction(math.log(eta_list[i]))))
    print(row)
    row = np.reshape(row,(1,size))
    lattice = np.concatenate((lattice,row),axis=0)
    return lattice
    
                   

def projection_scale(u, v):
    '''Computes <u,v>/<u,u>, which is the scale used in projection.'''
    return np.dot(u, v) / np.dot(u, u)

def proj(u, v):
    '''Computes the projection of vector v onto vector u. Assumes u is not zero.'''
    print(projection_scale(u,v))
    return projection_scale(u,v)*u

def gram_schmidt(basis): 
    '''Computes Gram Schmidt orthoganalization (without normalization) of a basis.'''
    orthobasis = basis.copy()
    orthobasis[:,0] = basis[:,0]
    mu = np.empty((np.shape(basis)))
    for i in range(1, np.shape(basis)[1]):  # Loop through dimension of basis.
        orthobasis[:,i] = basis[:,i] 
        for j in range(0, i):
            orthobasis[:,i] -= proj(orthobasis[:,j], basis[:,i])
            mu[i][j] = projection_scale(orthobasis[:,j],basis[:,i])
    return orthobasis,mu

def initialize(basis):
    D = get_D(basis)
    GS,mu = gram_schmidt(basis)
    Lambda = basis.copy()
    for i in range(np.shape(basis)[0]):
        for j in range(np.shape(basis)[1]):
            Lambda[i][j] = D[j]*mu[i][j]
    return Lambda,D
    
    

def get_D(basis):
    D = np.empty((np.shape(basis)[0]+1,1),dtype=np.float64)
    ortho = basis.copy()
    GS = gram_schmidt(basis)[0]  
    #print(GS)
    for i in range(np.shape(GS)[0]+1):
        if i == 0:
            D[i] = 1
        else:
            D[i] = reduce(operator.mul,[np.linalg.norm(GS[:,j])**2 for j in range(i)])
    print(D)
    return D

#===================================================

import sys
import json
import numpy as np
from numpy import linalg as la

#k = 1 # Initialize the working index.
#DELTA = 0.75 



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
            Lambda[i,j] = basis.column(i).dot_product(casis.column(j))
            casis[:,i] = (D[j+1]*casis.column(i) - Lambda[i][j]*casis.column(j))/D[j]
            print(casis)
        D.append((casis.column(i).dot_product(casis.column(i)))/D[i])
    print('D',D)
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
    print('a',Lambda)
    return k,l,basis,Lambda,D
        
def prod_B(k,basis,Lambda,D):
    basis[:,[k-2,k-1]] = basis[:,[k-1,k-2]]
    for j in range(0,k-2):
        Lambda[[k-2,k-1],j] = Lambda[[k-1,k-2],j]
    for i in range(k,np.shape(basis)[0]):
        t = Lambda[i,k-2]
        s = (Lambda[i,k-2]*Lambda[k-1,k-2] + Lambda[i,k-1]*D[k-2])/D[k-1]
        print(s.is_integer())
        Lambda[i,k-2] = (Lambda[i,k-2]*Lambda[k-1,k-2] + Lambda[i,k-1]*D[k-2])/D[k-1]
        Lambda[i,k-1] = (t*D[k] - Lambda[i,k-1]*Lambda[k-1,k-2])/D[k-1]
    D[k-1] = (D[k-2]*D[k] + Lambda[k-1,k-2]**2)/D[k-1]
    return k,basis,Lambda,D
    
def main(basis):
    Lambda,D = initial(basis)
    print('L',Lambda)
    print('D', D)
    k=2
    while k <= np.shape(basis)[0]: 
        l = k-1
        print('doing prod A')
        k,l,basis,Lambda,D = prod_A(k,l,basis,Lambda,D)
        print(Lambda)
        if 4*D[k-2]*D[k] < (3*D[k-1]**2 - 4*Lambda[k-1,k-2]**2):
            print('doing prod B',k)
            k,basis,Lambda,D = prod_B(k,basis,Lambda,D)
            print('D',D)
            if k > 2:
                k = k-1
        else:
            print('doing the other thing',k)
            for l in range(k-2,2):
                print('doing prod A')
                k,l,basis,Lambda,D = prod_A(k,l,basis,Lambda,D)
                #print(Lambda)
            k = k + 1
    return basis
    
if __name__ == "__main__":
    large_constant = 10**(10)
    primes = [2,3,5]
    basis = Matrix(ZZ,generate_lattice(large_constant,primes))
    print(main(basis))



