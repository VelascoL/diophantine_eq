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

def get_primes_list(min_num, max_num):
    return list(primes(max_num))


def LLL_real_reduce(bound,x_list,primes):
    size = len(primes)
    large_constant = 10
    approximation_matrix = Matrix(ZZ,generate_lattice_approx(large_constant,primes))
    print(approximation_matrix)
    #approximation_matrix = approximation_matrix.transpose()
    #find LLL reduced basis
    #print(approximation_matrix)
    B = approximation_matrix.LLL()
    print(B)
    #gram schmidt approximation
    GS,C = B.gram_schmidt()
    print(C)
    #print(GS)
    #find c_1
    vy = [0 for _ in range(size)]
    vy[-1] = -math.floor(Decimal(large_constant)*Decimal(math.log(math.sqrt(5))))
    y = vector(ZZ,vy)
    print(y)
    #v = y.transpose()
    sig = calculate_sigma(B,y) 
    frac = abs(sig.numer()) % (abs(sig.denom()))
    sig = N(frac / sig.denom())
    print(sig)
    c_2 = max([((B[0].norm())**2)/((GS[i].norm())**2) for i in range(size)])
    print('c2',c_2)
    print('b0',B[0].norm()**2)
    c_1 = (1/c_2)*sig*B[0].norm()**2
    print(c_1)
    test = minimal_vector(approximation_matrix,y)
    print(test)
    #calculate values needed for Lemma 6
    S,T = calculate_S_and_T(x_list,bound)
    print(T**2 + S)
    if c_1**2 < T**2 + S:
    #    raise ValueError("Need to choose larger C")
    #define constants 
    c_3 = 2*(1 + k*abs(b)/abs(a))
    c_4 = math.log(min(abs(alpha), abs(alpha)/abs(beta)))
    #get bound on H
    new_bound = (1/c_4)*(math.log(large_constant*c_3) - math.log(math.sqrt(c_1**2 - S) - T))
    if new_bound < 0:
        raise ValueError("new bound is a negative number.")
    return new_bound

def LLL_real_reduce1(bound,x_list,primes):
    size = len(primes)
    large_constant = 10**(1000)
    approximation_matrix = generate_lattice_approximation_matrix(large_constant,primes)
    approximation_matrix = approximation_matrix.transpose()
    #find LLL reduced basis
    B = LLL.reduction(IntegerMatrix.from_matrix(approximation_matrix))
    LLL.is_reduced(B)
    #B_matrix = [[0 for _ in range(B.nrows)] for _ in range(B.ncols)]
    #B.to_matrix(B_matrix)
    
    #B_cols = get_column_vectors(B_matrix)
    
    #GS = GSO.Mat(B,flags=GSO.INT_GRAM); _ = GS.update_gso()
    #GSM = GS.G
    #GS_matrix = [[0 for _ in range(GSM.nrows)] for _ in range(GSM.ncols)]
    #GSM.to_matrix(GS_matrix)
    
    #GS_cols =get_column_vectors(GS_matrix)
    #print(GS_cols)
   
    #vy = [0 for _ in range(size)]
    #vy[-1] = -math.floor(Decimal(large_constant)*Decimal(math.log(math.sqrt(5))))
    #c2 = calculate_max_frac_norm_squared(B_cols,GS_cols)
    #sigma = calculate_sigma1(B_matrix,vy)
    #print('sig',sigma)
    #print('norm',calculate_norm_squared(B_cols[0]))
    #print('c2',c2)
    #c1 = sigma*calculate_norm_squared(B_cols[0])/c2
    #print('c1',c1)
    #calculate values needed for Lemma 6
    #S,T = calculate_S_and_T(x_list,bound)
    #print(T**2 + S)
    #if c1**2 < T**2 + S:
    #    raise ValueError("Need to choose larger C")
    #define constants 
    #c_3 = 2*(1 + k*abs(b)/abs(a))
    #c_4 = math.log(min(abs(alpha), abs(alpha)/abs(beta)))
    #get bound on H
    #new_bound = (1/c_4)*(math.log(large_constant*c_3) - math.log(math.sqrt(c_1**2 - S) - T))
    #if new_bound < 0:
    #    raise ValueError("new bound is a negative number.")
    #return new_bound



def LLL_padic_reduce(test):
    return

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
    
def calculate_sigma1(matrix, v):
    """
    Calculate sigma as defined in the paper.
    """
    # STEP 1: Calculate the vector z
    #print(matrix)
    np_inverted_matrix = np.linalg.inv(np.array(matrix,dtype=np.float64))
    print(np_inverted_matrix)
    z = np.dot(np_inverted_matrix, v)

    # STEP 2: Find the largest index such that the entry is non-zero
    last_index = find_last_nonzero_index(z)

    # STEP 3: Calculate the distance to the nearest integer
    # TODO: Fix floating point errors.
    return calculate_distance_to_nearest_int(z[last_index])

def calculate_sigma(matrix, v):
    """
    Calculate sigma as defined in the paper.
    """
    # STEP 1: Calculate the vector z
    inverted_matrix = matrix.inverse()
    z = inverted_matrix*v
    print('z',z)

    # STEP 2: Find the largest index such that the entry is non-zero
    last_index = find_last_nonzero_index(z,primes)

    # STEP 3: Calculate the distance to the nearest integer
    # TODO: Fix floating point errors.
    return z[last_index]#calculate_distance_to_nearest_int(z[last_index])

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

def calculate_S_and_T(x_list,c20):
    """
    Calculates S and T and returns it as a tuple (S, T)
    """
    # STEP 1: Calculate S
    S = sum([x ** 2 for x in x_list])
    
    # STEP 2: Calculate T
    T = ((sum(x_list) + c20)/2)
    
    return (S, T)

def calculate_norm_squared(v):
    return sum([x ** 2 for x in v])

def calculate_max_frac_norm_squared(lll_reduced_basis, gs_basis):
    return max([calculate_norm_squared(lll_reduced_basis[0])/calculate_norm_squared(v) for v in gs_basis])

def calculate_lattice_min_distance_squared(lll_reduced_basis, vy):
    """
    Calculates l(L, y)^2, where y is a vector.
    """
    distances = []
    for v in lll_reduced_basis:
        distances.append([x - y for x, y in zip(v, vy)])
    lattice_min_distance = min(map(calculate_norm_squared, distances))
    
    print("lattice min distance: %d" % lattice_min_distance)
    # TODO: Check if vy is in the lattice.
    return lattice_min_distance

def calculate_euclidean_distance_squared(vx, vy):
    """
    Calculates the Euclidean distance
    """
    return sum([(x - y) ** 2 for x, y in zip(vx, vy)])

def get_column_vectors(matrix):
    """
    Extract the column vectors from a matrix (in list representation)
    """
    column_vectors = []
    for i in range(len(matrix)):
        col = []
        for j in range(len(matrix)):
            col.append(matrix[j][i])
        column_vectors.append(col)
    return column_vectors

def calculate_new_bound():
    return

def calculate_d1():
    return

def calculate_constants(alpha, beta, a, b):
    """
    Calculates constants that are defined within our paper.
    """
    constants = {}
    constants["c1"] = 2
    constants["c2"] = 2
    constants["c3"] = 2
    constants["c4"] = 2
    constants["c5"] = 2

    constants["d1"] = 2

    return constants

def generate_lattice_approximation_matrix(large_constant, primes):
    """
    Generates the lattice approximation matrix.
    """
    size = len(primes)
    approximation_matrix = [[0 for _ in range(size - 1)] for _ in range(size - 1)]
    IntegerMatrix.identity(size - 1).to_matrix(approximation_matrix)

    # Append an empty column to the identity matrix
    for row in approximation_matrix:
        row.append(0)

    # Create the row that represents the approximation of the linear form
    row_approx = []
    for i in range(len(primes)):
        row_approx.append(math.floor(large_constant * primes[i]))
    approximation_matrix.append(row_approx)
    return approximation_matrix

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
        primes = [2,3,5,7]
    )

    primes = [2,3,5,7]
    c = constants.calculate_constants()
    nbound = c['n1_bound']
    Z_list = c['Z_bounds']
    new_bound = LLL_real_reduce(nbound,Z_list,primes)
    print(new_bound)
    print("Handling n1 = ... = nk case...")
    print("Done.")

    #print("Performing first reduction case...")
    #LLL_real_reduce(100, [1, 2, 3])
    