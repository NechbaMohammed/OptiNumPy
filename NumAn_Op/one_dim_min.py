#=================================================\\Documentation//==================================================

"""
This module contains different types of optimization methods that search for the min of a one dimensional unimodal function, 
using both the search with elimination and search with interpolation methods.

"""

#================================================\\Libraries needed//================================================

import numpy as np
import numdifftools as nd
import math as m


#=======================================\\Searching with elimination methods//=======================================

# 1.UNRESTRICTED SEARCH:

# Search with fixed step size:


def fixed_step_size( x0 , p ,f ) :
    """ Minimisation de la fonction scalaire d'une variable réel.
    
    Parameters
    ----------
    x0 : float
       Point initiale.
    
    p : float
      Le pas.
    f : callable
      La fonction objective à minimiser.
      
         ``f( x ) -> float``
         
      où ``x`` est un nombre réel 
      
    Returns
    -------
    res : float
       Le résultat de l'optimisation.
    List_sol_iteration : list
                        liste des solution 

    """    
   
    f0 = f( x0 )
    if f0 > f( x0 + p ) :
        i = 1
    if f0 > f( x0 - p ) :
        i = -1
 
    x1 = x0 + i * p
    f1 = f( x1 )
    List_sol_iteration = [ x0 ]
    while f1 < f0 :
         x0 = x1
         x1 = x0 + i * p
         f0 = f1
         f1 = f( x1 )
         List_sol_iteration.append( x0 )
    res = x0
    return res,List_sol_iteration

# Search with accelerated step size:

def accelerated_step_size( x0 , p , f ) :
         """ Minimisation de la fonction scalaire d'une variable réel.
         
         Parameters
         ----------
         x0 : float
            Point initiale.
         
         p : float
           Le pas.
         f : callable
           La fonction objective à minimiser.
           
              ``f( x ) -> float``
              
           où ``x`` est un nombre réel 
           
         Returns
         -------
         res : float
            Le résultat de l'optimisation.
         List_sol_iteration : list
                             liste des solution 

         """
         
         f0 = f( x0 )
         x1 = x0
         if f0 > f( x0 + p ) :
            i = 1
         if f0 > f( x0 - p ) :
            i = -1
         List_sol_iteration = [ x0 ]
         while  round( abs( x1 - x0 ) , 2 ) != p :
            k = i * p
            x1 = x0 + k
            f1 = f( x1 )
            while f1 < f0 :
              k = 2 * k
              x0 = x1
              x1 = x0 + k
              f0 = f1
              f1 = f( x1 )
              List_sol_iteration.append( x0 )
            res = x0
         return res,List_sol_iteration
     
# 2.EXHAUSTIVE SEARCH:
    
def Exhaustive_search( a , b , n ,f ) :
            """ Minimisation de la fonction scalaire d'une variable réel.
            

            Parameters
            ----------
            a : float
                Bornes inférieure de domine de recherche.
            b : float
                Bornes supérieure de domine de recherche.
            n : int
                Nombre de point evalue.
            f : callable
                La fonction objective à minimiser.
              
                 ``f( x ) -> float``
                 
              où ``x`` est un nombre réel 

            Returns
            -------
            res : float
               Le résultat de l'optimisation.
            List_sol_iteration : list
                                 liste des solution 


            """
           
            x1 = a
            h = ( b - a ) / n
            x2 = x1 + h
            x3 = x2 + h
            f1 = f( x1 )
            f2 = f( x2 )
            f3 = f( x3 )
            List_sol_iteration = [ ]
            while 1:
                List_sol_iteration.append( ( x1+ x3)/2)

                if f2 <= f1 and f2 <= f3 :
                     res = ( x1 +x3 ) / 2
                     return res , List_sol_iteration
                else:
                    x1 = x2
                    f1 = f2
                    x2 = x3
                    f2 = f3
                    x3 = x2 + h
                    f3 = f( x3 )
                    

# 3.DICHOTOMOUS SEARCH:
    
def Dichotomous_search( a , b , f , tol = 1e-3 ) :
            """ Minimisation de la fonction scalaire d'une variable réel.
            

            Parameters
            ----------
            a : float
                Bornes inférieure de domine de recherche.
            b : float
                Bornes supérieure de domine de recherche.
            f : callable
                La fonction objective à minimiser.
              
                 ``f( x ) -> float``
                 
              où ``x`` est un nombre réel 
              
            tol : float, optional
                Tolérance pour la terminaison. The default is 1e-3.

            Returns
            -------
            res : float
               Le résultat de l'optimisation.
            List_sol_iteration : list
                                  liste des solution
               
            """

            a1 = a
            b1 = b
            List_sol_iteration = [ ]
            while round( b1 - a1 , 2 ) > tol :
                m = ( a1 + b1 ) / 2
                List_sol_iteration.append( m )
                x1 = m - tol / 2
                x2 = m + tol / 2
                if f( x1 ) < f( x2 ) :
                    b1 = x2
                else:
                    a1 = x1

            res=( a1 + b1 ) / 2
            return res , List_sol_iteration


# 4.INTERVAL HALVING METHOD:

    
def Interval_Halving_Method( a , b , f , tol = 1e-3 ) :
            """ Minimisation de la fonction scalaire d'une variable réel.
            

            Parameters
            ----------
            a : float
                Bornes inférieure de domine de recherche.
            b : float
                Bornes supérieure de domine de recherche.
            f : callable
                La fonction objective à minimiser.
              
                 ``f( x ) -> float``
                 
              où ``x`` est un nombre réel 
            tol : float, optional
                Tolérance pour la terminaison. The default is 1e-3.

            Returns
            -------
            res : float
               Le résultat de l'optimisation.
            List_sol_iteration : list
                                   liste des solution
                
               
            """
            
            m = ( a + b ) / 2
            List_sol_iteration = [ ]
            L = b - a
            f0 = f( m )
            while 1 :
                List_sol_iteration.append( m )
                x1 = a + L / 4
                x2 = b - L / 4
                f1 = f( x1 )
                f2 = f( x2 )
                if f1 < f0 :
                    b = m
                    m = x1
                    f0 = f1
                elif f2 < f0 :
                    a = m
                    m = x2
                    f0 = f2
                else:
                    a = x1
                    b = x2
                L = b - a
                if abs( L ) < tol :
                     res = ( a + b ) / 2
                     return res , List_sol_iteration
                    
# 5.FIBONACCI METHOD:

def suite_Fibonacci( n ) :
            """
             

             Parameters
             ----------
             n : int
                 Le nombre de termes.

             Returns
             -------
             res : int
                   n ième terme d'une suite.

            """
            f0 = f1 = 1
            for i in range( 2 , n + 1 ):
               f = f0 + f1
               f0 = f1
               f1 = f
            res = f1
            return res



def Fibonacci_method( n , a , b , f ) :
            """ Minimisation de la fonction scalaire d'une variable réel.
            

            Parameters
            ----------
            n : int
                Le nombre d'itérations.
            a : float
                Bornes inférieure de domine de recherche.
            b : float
                Bornes supérieure de domine de recherche.
            f : callable
                La fonction objective à minimiser.
              
                 ``f( x ) -> float``
                 
              où ``x`` est un nombre réel

            Returns
            -------
            res : float
               Le résultat de l'optimisation.
            List_sol_iteration : list
                                  liste des solution    

            """
           
            l = b - a
            k = 1
            List_sol_iteration = [ ]
            while k != n :
                List_sol_iteration.append( (a+b)/2)
                k = k + 1
                L= ( suite_Fibonacci( n - k + 1 ) / suite_Fibonacci( n + 1 ) ) * l
                x1 = a + L
                x2 = b - L
                f1 = f( x1 )
                f2 = f( x2 )
                if f1 < f2 :
                    b = x2
                if f2 < f1 :
                    a = x1
                if f1 == f2 :
                    a = x1
                    b = x2
                    k = k + 1
            res = ( a + b ) / 2
            return res , List_sol_iteration

# 6.GOLDEN SECTION METHOD:
    
    
def Golden_section( a , b , f , tol = 1e-3 ) :
            """ Minimisation de la fonction scalaire d'une variable réel.
            

            Parameters
            ----------
            a : float
                Bornes inférieure de domine de recherche.
            b : float
                Bornes supérieure de domine de recherche.
            f : callable
                La fonction objective à minimiser.
              
                 ``f( x ) -> float``
                 
                   où ``x`` est un nombre réel 
            tol : float, optional
                Tolérance pour la terminaison. The default is 1e-3.

            Returns
            -------
            res : float
               Le résultat de l'optimisation.
            List_sol_iteration : list
                                   liste des solution    
            """
            inv_nbr_or = 0.618 # L'inverse du nombre d'or
            a1 = a
            b1 = b
            x1 = b1 - inv_nbr_or * ( b1 - a1 )
            x2 = a1 + inv_nbr_or * ( b1 - a1 )
            f1 = f( x1 )
            f2 = f( x2 )
            List_sol_iteration = [ ]
            while 1 :
                List_sol_iteration.append((a1+b1)/2)
                if f2 > f1:
                    b1 = x2
                    x2 = x1
                    f2 = f1
                    x1 = b1 - inv_nbr_or * ( b1 - a1 )
                    f1 = f( x1 )
                else:
                    a1 = x1
                    x1 = x2
                    f1 = f2
                    x2 = a1 + inv_nbr_or * ( b1 - a1 )
                    f2 = f( x2 )
                if round( b1 - a1 , 3 ) <= tol :
                        res = ( a1 + b1 ) / 2
                        return  res , List_sol_iteration



#======================================\\Searching with interpolation methods//======================================

# 1.NEWTON-RAPHSON METHOD:


def Newton_method( x0 , fp , fpp , tol = 1e-3 ) :
    """ Minimisation de la fonction scalaire d'une variable réel.
    
    Parameters
    ----------
    x0 : float
        Point initiale.
    fp : callable
        La derivee premiere de la fonction objective.
        
        ``fp( x ) -> float``
        
          où ``x`` est un nombre réel 
    fpp : callable
        La derivee seconde de la fonction objective.
        
        ``fpp( x ) -> float``
        
          où ``x`` est un nombre réel 
    tol : float, optional
        Tolérance pour la terminaison. The default is 1e-3.

    Returns
    -------
    res : float
       Le résultat de l'optimisation.
    List_sol_iteration : list
                           liste des solution

    """
    List_sol_iteration = [ x0 ]
    f0 = fp( x0 )
    while abs( f0 ) > tol :
        x = x0 - f0 / fpp ( x0 )
        x0 = x
        f0 = fp( x0 )
        List_sol_iteration.append(x0 )
    res = x0
    return res , List_sol_iteration


# 2.QUASI-NEWTON METHOD:

def Quasi_Newton( x0 , f , tol = 1e-3 ) :
    """ Minimisation de la fonction scalaire d'une variable réel.
    
    Parameters
    ----------
    x0 : float
        Point initiale.
     f : callable
         La fonction objective à minimiser.
       
          ``f( x ) -> float``
          
            où ``x`` est un nombre réel 
     tol : float, optional
         Tolérance pour la terminaison. The default is 1e-3.

     Returns
     -------
     x : float
        Le résultat de l'optimisation.
     List_sol_iteration : list
                            liste des solution


    """
  
    h=1e-2
    List_sol_iteration = [ ]
    while True:
        List_sol_iteration.append( x0 )
        f0 = f( x0 )
        f1 = f( x0 - h )
        f2 = f( x0 + h )
        x = x0 - h * ( f2 - f1 ) / ( 2 * ( f2 - 2 * f0 + f1 ) )
        x0 = x
        if abs ( ( f2 - f1 ) / ( 2 * h ) ) <= tol :
            return x , List_sol_iteration
        
        
# 3.SECANT METHOD:
    
    
def Secant_method( x0  , fp , tol = 1e-3 ) :
    """ Minimisation de la fonction scalaire d'une variable réel.
    
    Parameters
    ----------
    x0 : float
        Point initiale.
     fp : callable
        La derivee premiere de la fonction objective.
       
          ``fp( x ) -> float``
          
            où ``x`` est un nombre réel 
     tol : float, optional
         Tolérance pour la terminaison. The default is 1e-3.

     Returns
     -------
     x : float
        Le résultat de l'optimisation.
    List_sol_iteration : list
                           liste des solution

    """
    A = x0
    h = 1e-2
    while fp( h ) < 0 :
      A = h
      h = 2 * h
    B = h
    fA = fp( A )
    fB = fp( B )
    List_sol_iteration = []
    while True :
        x = A - fA * ( B - A ) / ( fB - fA )
        List_sol_iteration.append( x )
        fx = fp( x )
        if abs( fx ) <= tol :
          return x , List_sol_iteration
        if fx < 0 :
          A = x
          fA = fp( A )
        else:
          B = x
          fB = fp( B )
