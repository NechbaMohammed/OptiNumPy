#=================================================\\Documentation//==================================================

"""
This module contains some types of numerical analysis methods for solving systems of equations & matrix decompositions, 
such as Gauss-Jordan, LU decomposition and Cholesky decomposition.

"""

#================================================\\Libraries needed//================================================

import numpy as np
import numdifftools as nd
import math as m

#=========================================\\The Elimination Of Gauss-Jordan//========================================

# 1.The Inverse of a Matrix:
    
def inverse( A ):
    """
    Calculer l'inverse d'une matrice carrée.

    Parameters
    ----------
    A : (N, N) array_like
        Matrice a invirse.

    Returns
    -------
    A_1 : (N, N) array_like
        L'inverse de A.

    """
    
    N = len( A )
    I = np.identity( N )
    A = np.c_[ A , I ]
    for n in range(0,N):
        p = A[ n ][ n ]
        for i in range(n+1,N):
            x = A[ i ][ n ]
            for j in range(n,2*N):
                A[ i ][ j ] = A[ i ][ j ] - x * A[ n ][ j ] / p
    for n in range(N-1,0,-1):
        p = A[ n ][ n ]
        for i in range(n-1,-1,-1):
            x = A[ i ][ n ]
            for j in range(n,2*N):
                A[ i ][ j ] = A[ i ][ j ] - x * A[ n ][ j ]/ p
    for i in range(0,N):
        p = A[ i ][ i ]
        for j in range(i,2*N):
            A[ i ][ j ] = A[ i ][ j ] / p
    A_1 = I   
    for i in range(0,N):
        for j in range(0,N):
            A_1[ i ][ j ] = A[ i ][ j+N ]
    return A_1


# 2.Gauss_Jordan method to solve system of linear equations:


def Gauss_Jordan_to_Solvea_system( A , b ):
    """
    Résoudre un système d'équations, A x = b.

    Parameters
    ----------
    A : (N, N) array_like
        Matrice de système.
    b : ( N ) array
        Côté droit.

    Returns
    -------
    x : ( N ) array
        Solution au système.

    """
    A = np.c_[ A , b ]
    N = len( A )
    for n in range(0,N):
        p = A[ n ][ n ]
        for i in range(n+1,N):
            x = A[ i ][ n ]
            for j in range(n,N+1):
                A[ i ][ j ] = A[ i ][ j ] - x * A[ n ][ j ] / p
    for n in range(N-1,0,-1):
        p = A[ n ][ n ]
        for i in range(n-1,-1,-1):
            x = A[ i ][ n ]
            for j in range(n,N+1):
                A[ i ][ j ] = A[ i ][ j ] - x * A[ n ][ j ] / p
    for i in range(0,N):
        p = A[ i ][ i ]
        for j in range(i,N+1):
            A[ i ][ j ] = A[ i ][ j ] / p
    x=np.zeros(N)
    for i in range(0,N):
        x[ i ] = A[ i ][ N ]
    return x

# =============================================================================
# La décomposition PA=LU d'un matrice 
# =============================================================================

def permutation( A , n , N ) :
    T = np.identity( N )
    m = A[ n ][ n ]
    j = n
    for i in range(n+1,N):
        if m < abs( A[ i ][ n ] ):
            m = A[ i ][ n ]
            j = i
    T[ n ][ n ] = 0
    T[ j ][ n ] = 1
    T[ j ][ j ] = 0
    T[n][j]=1
    A0 = T @ A
    return T,A0

def Decomposition_LU_PA( A ) :
    """
    Calculer la décomposition LU d'une matrice PA.

    La décomposition est : :
       P A = L U
    où P est une matrice de permutation, L triangulaire inférieure avec des éléments diagonaux unitaires et U triangulaire supérieure.

    Parameters
    ----------
    A : (N, N) array_like
        Matrice à décomposer.

    Returns
    -------
    P : (N, N) ndarray
        Matrice de permutation.
    L : (N, N) ndarray
        Matrice triangulaire inférieure.
    U : (N, N) ndarray
        Matrice triangulaire supérieure.

    """
    N = len( A )
    A0 = A
    P = L = np.identity( N )
    n = 0
    while n < N-1:
        T,A0 = permutation( A0 , n , N )
        P = T @ P
        L = T @ L @ T
        L0 = np.identity( N )
        for i in range( n+1 , N ):
            L0[ i ][ n ] = -A0[ i ][ n ] / A0[ n ][ n ]
        A0 = L0 @ A0
        L = L @ np.linalg.inv( L0 )
        n = n+1
    U = A0
    return P,L,U  

# =============================================================================
# Résolution d'un système d'équations linéaires en utilisant la décomposition QLU 
# =============================================================================

def Dec_PA_LU_to_Solvea_system( A , b ) :
    """
    Résoudre un système d'équations, A x = b, en utilisant la factorisation LU de P A.

    Parameters
    ----------
    A : (N, N) array_like
        Matrice de système.
    b : ( N ) array
        Côté droit.

    Returns
    -------
    X : ( N ) array
        Solution au système.

    """
    P,L,U = Decomposition_LU_PA( A ) 
    N = len( A )
    b = P @ b
    Y = X = np.zeros( N )
    for i in range(0,N):
        Y[ i ] = b[ i ]
        for j in range(0,i):
            Y[ i ] = Y[ i ] - Y[ j ] * L[ i ][ j ]
        Y[ i ] = Y[ i ] / L[ i ][ i ]
    for i in range(N-1,-1,-1):
        X[ i ] = Y[ i ]
        for j in range(i+1,N):
            X[ i ] = X[ i ] - X[ j ] * U[ i ][ j ]
        X[ i ] = X[ i ] / U[ i ][ i ]
    return X

# =============================================================================
# La décomposition de Choleski d'un matrice
# =============================================================================

def determinant_sous_matrice( A , N ) :
    SA = np.zeros(( N , N ))
    for i in range(0,N):
            for j in range(0,N):
                SA[ i ][ j ] = A[ i ][ j ]
    return np.linalg.det( SA ) 

def define_positif( A ) :
    N = len( A )
    for n in range(1,N+1):
        if determinant_sous_matrice( A , n ) <= 0 :
            return False
    return True

def symetrique( A ) :
    N = len( A )
    for i in range(0,N):
        for j in range(i+1,N):
            if( A[ i ][ j ] != A[ j ][ i ] ) :
              return False
    return True

def Dec_Choleski( A ) :
    """
    Calculer la décomposition ``A =L L.T`` où L triangulaire inférieure et L.T la transposée de L.

    Parameters
    ----------
    A : (N, N) array_like
        Matrice à décomposer.

    Returns
    -------
    L : (N, N) array_like
        Matrice triangulaire inférieure de la décomposition.

    """
    if symetrique( A ) :
        if define_positif( A ) :
            N = len( A )
            L = np.zeros(( N , N ))
            for i in range(0,N):
                for k in range(0,i):
                    L[ i ][ i ] = L[ i ][ i ] + L[ i ][ k ] * L[ i ][ k ]
                L[ i ][ i ] = A[ i ][ i ] - L[ i ][ i ]
                L[ i ][ i ] = np.sqrt( L[ i ][ i ] )
                for j in range(i+1,N):
                    L[ j ][ i ] = A[ i ][ j ]
                    for k in range(0,i):
                       L[ j ][ i ] = L[ j ][ i ] - L[ i ][ k ] * L[ j ][ k ]
                    L[ j ][ i ] = L[ j ][ i ] / L[ i ][ i ]

            return L
        else:
            print("Donnes un matrice define positif\n")
    else:
          print("Donnes un matrice symetrique\n")

# =============================================================================
# Résolution d'un système d'équations linéaires en utilisant la décomposition de Choleski
# =============================================================================

def Choleski_to_Solvea_system( A , b ) :
    """
    Résoudre un système d'équations, A x = b, en utilisant la décomposition de Choleski de A.

    Parameters
    ----------
    A : (N, N) array_like
        Matrice de système.
    b : ( N ) array
        Côté droit.

    Returns
    -------
    X : ( N ) array
        Solution au système.

    """
    L = Dec_Choleski( A )
    LT = np.transpose( L )
    N = len( A )
    Y = X = np.zeros( N )
    for i in range(0,N):
        Y[ i ] = b[ i ]
        for j in range(0,i):
            Y[ i ] = Y[ i ] - Y[ j ] * L[ i ][ j ]
        Y[ i ] = Y[ i ] / L[ i ][ i ]
    for i in range(N-1,-1,-1):
        X[ i ] = Y[ i ]
        for j in range(i+1,N):
            X[ i ] = X[ i ] - X[ j ] * LT[ i ][ j ]
        X[ i ] = X[ i ] / LT[ i ][ i ]
    return X

