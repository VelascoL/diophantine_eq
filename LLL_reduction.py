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
from LLLnpDW import LLLDW

def get_primes_list(min_num, max_num):
    return list(primes(max_num))

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
    
def calculate_norm_squared(v):
    return v.dot_product(v)

def calculate_max_frac_norm_squared(lll_reduced_basis, gs_basis):
    return max([calculate_norm_squared(lll_reduced_basis.column(0))/calculate_norm_squared(gs_basis.column(i)) for i in range(size)])

def calculate_distance_to_nearest_int(x):
    """
    Used in calculating sigma.
    """
    y = mp.fabs(x)
    frac_part = mp.fsub(y, mp.nint(y))
    if frac_part == 0:
        return ValueError("frac_part is 0.")
    return min(frac_part, mp.fsub(1, frac_part))

def find_last_nonzero_index(v,primes):
    """
    Used in calculating sigma.
    """ 
    last_index = math.inf
    for i in range(len(primes)):
        test = v[-(i + 1)]
        frac = abs(test.numer()) % (abs(test.denom()))
        t = N(frac / test.denom())
        if t != 0:
            return len(primes) - (i + 1)
    raise ValueError("Zero vector passed in.")

def calculate_sigma(basis, v):
    """
    Calculate sigma as defined in the paper.
    """
    # STEP 1: Calculate the vector z
    inverted_matrix = basis.inverse()
    z = inverted_matrix*v
    #print('z',z)

    # STEP 2: Find the largest index such that the entry is non-zero
    last_index = find_last_nonzero_index(z,primes)

    # STEP 3: Calculate the distance to the nearest integer.
    return z[last_index]#calculate_distance_to_nearest_int(z[last_index])

def calculate_S_and_T(x_list,c20):
    """
    Calculates S and T and returns it as a tuple (S, T)
    """
    # STEP 1: Calculate S
    S = sum([x ** 2 for x in x_list])
    
    # STEP 2: Calculate T
    T = (sum(x_list) + c20)/2
    
    print('S',S)
    
    return (S, T)


def get_bound(bound,x_list,primes):
    size = len(primes)
    L = generate_approximation_matrix(large_constant)
    B = LLLDW(L)
    #print('reduced',B)
    #do gram-schmidt
    GS = B.transpose().gram_schmidt()[0] #tranpose bc sage is weird
    GS = Matrix(QQ,GS)
    GS = GS.transpose()
    #find c_2
    c_2 = calculate_max_frac_norm_squared(B,GS)
    #print('c_2',n(c_2))
    vy = [0 for _ in range(size)]
    vy[-1] = -math.floor(Decimal(large_constant)*Decimal(math.log(math.sqrt(5))))
    y = vector(ZZ,vy)
    sig = calculate_sigma(B,y)
    frac = abs(sig.numer()) % (abs(sig.denom()))
    sig = N(frac/sig.denom())
    sig = min(sig, 1-sig)
    #calculate c_1
    c_1 = (1/n(c_2))*1*calculate_norm_squared(B.column(0))
    #print('c_1',c_1)
    S,T = calculate_S_and_T(x_list,bound)
    print("S",S)
    if c_1 < T**2 + S:
        #print(T**2 + S)
        raise ValueError("Need to choose larger C")
    c_3 = 2*(1 + k*abs(b)/abs(a))
    c_4 = math.log(min((abs(alpha)/abs(beta), abs(alpha))))
    #get bound on H
    new_bound = (1/c_4)*(math.log(large_constant*c_3) - math.log(math.sqrt(c_1 - S) - T))
    if new_bound < 0:
        raise ValueError("new bound is a negative number")
    return new_bound
    

if __name__ == "__main__":
    constants = Constants(
        a = 1,
        b = 1,
        A = 1,
        B = 1,
        delta = 5,
        alpha = (1 + sqrt(5))/2,
        beta = (1 - sqrt(5))/2,
        num_terms = 2,
        w = 1,
        primes = [2,3,5]
    )
    k = 3
    b = 1
    a = 1
    alpha = (1 + sqrt(5))/2
    beta = (1 - sqrt(5))/2
    primes = [2,3,5]
    size = len(primes)
    large_constant = 10**(100)
    c = constants.calculate_constants()
    nbound = c['n1_bound']
    Z_list = c['Z_bounds']
    print(Z_list)
    
    diff_bound = get_bound(nbound,Z_list,primes)
    print(diff_bound)
    
    
    print("Handling n1 = ... = nk case...")
    print("Done.")

    #print("Performing first reduction case...")
    #LLL_real_reduce(100, [1, 2, 3])
    