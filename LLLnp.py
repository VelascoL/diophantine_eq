# -*- coding: utf-8 -*-
"""
Solves the family of Diophantine equations of the form
u_{n_1} + ... + u_{n_k} = wp_1^{z_1} ... p_s^{z_s},
under some mild technical restrictions.
@authors: Brian Ha, Lily McBeath, Luisa Velasco
"""
import math
import numpy as np
from sage.all import *
from sage.rings.number_field.S_unit_solver import minimal_vector
from constants1 import Constants
from fractions import Fraction
from decimal import Decimal
import mpmath as mp
from fpylll import *

def generate_lattice_approx(large_constant,primes):
    """
    Generates the lattice approximation matrix not using fyplll
    """
    size = len(primes)
    
    approximation_matrix = np.identity(size-1)
    # Append an empty column to the identity matrix
    approximation_matrix = np.concatenate((approximation_matrix,np.zeros((1,size-1))),axis = 0)
    # Create the row that represents the approximation of the linear form
    row_approx = []
    for i in range(len(primes)):
        #if i == len(primes)-1:
         #   row_approx.append(-math.floor(Fraction(large_constant) * Fraction(math.log(primes[i]))))
        #else:
        row_approx.append(math.floor(Fraction(large_constant) * Fraction(math.log(primes[i]))))
    row_approx = np.reshape(row_approx,(size,1))
    print(np.shape(approximation_matrix))
    print(np.shape(row_approx), np.shape(approximation_matrix))
    approximation_matrix = np.concatenate((approximation_matrix, row_approx),axis = 1)
    return approximation_matrix



"""LLL.py: Implements the LLL Algorithm originally developed by Lenstra, Lenstra, Lovasz in 1982. Python 2.7."""
__author__ = "Raj Kane"
__version__ = "Spring 2019"

import sys
import json
import numpy as np
from numpy import linalg as la

k = 1 # Initialize the working index.
DELTA = 0.75   

def projection_scale(u, v):
    '''Computes <u,v>/<u,u>, which is the scale used in projection.'''
    return np.dot(u, v) / np.dot(u, u)

def proj(u, v):
    '''Computes the projection of vector v onto vector u. Assumes u is not zero.'''
    return np.dot(projection_scale(u, v), u)

def gram_schmidt(orthobasis,basis): 
    '''Computes Gram Schmidt orthoganalization (without normalization) of a basis.'''
    orthobasis[0] = basis[0]
    for i in range(1, basis.shape[1]):  # Loop through dimension of basis.
        orthobasis[i] = basis[i]
        for j in range(0, i):
            orthobasis[i] -= proj(orthobasis[j], basis[i])
    return orthobasis

def reduction(orthobasis,basis):
    '''Performs length reduction on a basis.'''
    total_reduction = 0 # Track the total amount by which the working vector is reduced.
    for j in range(k-1, -1, -1):   # j loop. Loop down from k-1 to 0.
        m = round(projection_scale(orthobasis[j], basis[k]))
        total_reduction += np.dot(m, basis[j])[0]
        basis[k] -= np.dot(m, basis[j]) # Reduce the working vector by multiples of preceding vectors.
    if total_reduction > 0:
        gram_schmidt(orthobasis,basis) # Recompute Gram-Scmidt if the working vector has been reduced. 

def lovasz(orthobasis,basis):
    global k
    '''Checks the Lovasz condition for a basis. Either swaps adjacent basis vectors and recomputes Gram-Scmidt or increments the working index.'''
    c = DELTA - projection_scale(orthobasis[k-1], basis[k])**2
    if la.norm(orthobasis[k])**2 >= np.dot(c, la.norm(orthobasis[k-1]**2)): # Check the Lovasz condition.
        k += 1  # Increment k if the condition is met.
    else: 
        basis[[k, k-1]] = basis[[k-1, k]] # If the condition is not met, swap the working vector and the immediately preceding basis vector.
        gram_schmidt(orthobasis,basis) # Recompute Gram-Schmidt if swap
        k = max([k-1, 1])

def main(lattice):
    basis = np.array(json.loads(lattice)).astype(float) # Initialize the basis as the user input.
    orthobasis = basis.copy()   # Initialize the Gram-Schmidt basis
    while True:
        x = input("Would you like to see the steps? Press [Y/N] and Enter. ")
        if x in ['Y','N', 'y', 'n']: break
        else: input("Would you like to see the steps? Press [Y/N] and Enter. ")
    if x in ['Y', 'y']:
        gram_schmidt(orthobasis,basis)
        steps = 0
        while k <= basis.shape[1] - 1:
            reduction(orthobasis,basis)
            steps += 1
            print('Step ', steps,'. After the reduction step, the basis is\n', basis)
            input("")
            lovasz(orthobasis,basis)
            steps +=1
            print('Step ', steps,'. After checking the Lovasz condition, the basis is\n', basis)
            input("")
        print('LLL Reduced Basis:\n', basis)
    else:
        gram_schmidt(orthobasis, basis)
        while k<= basis.shape[1] - 1:
            reduction(orthobasis,basis)
            lovasz(orthobasis,basis)
        print('LLL Reduced Basis:\n', basis)

if __name__ == "__main__":
    large_constant = 10
    primes = [2,3,5,7]
    basis = generate_lattice_approx(large_constant,primes)
    print(basis)
    b = basis.tolist()
    
    main(str(b))



