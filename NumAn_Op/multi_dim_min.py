#=================================================\\Documentation//==================================================

"""
This module contains different types of optimization methods that search for the min of a multi-dimensional function, 
such as Newton method, Gradient descent, Conjugate Gradient...

"""

#================================================\\Libraries needed//================================================

import numpy as np
import numdifftools as nd
import math as m

#=================================================\\Gradient methods//===============================================

# 1.Gradient Descent:

    
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
            return x 
        

def methode_gradient( x0 , f , tol =1e-3 ):
    """ Minimisation de la fonction scalaire d'une ou plusieurs variable.
    

    Parameters
    ----------
    x0 : ndarray, shape (n,)
        Point initiale.
    f : callable
         La fonction objective à minimiser.
       
          ``f( x ) -> float``
          
            x`` is an 1-D array with shape (n,)
    tol : float, optional
         Tolérance pour la terminaison. The default is 1e-3.

    Returns
    -------
     x : ndarray, shape (n,)
        Le résultat de l'optimisation.
     List_sol_iteration : list
                           liste des solution

    """
    df = nd.Gradient( f )
    grad = df( x0 )
    print(grad)
    g = lambda x , y : f( x - y * df( x ))
    List_sol_iteration = [ ]
    while True:
        List_sol_iteration.append( x0)
        g1 = lambda y : g( x0 , y )
        y = Quasi_Newton( 0 , g1 ) 
        print(y)

        x = x0 - y * grad
        grad = df( x )

        if np.linalg.norm( grad ) <= tol or  np.linalg.norm( x-x0 ) <= tol :
           return x , List_sol_iteration
        x0 = x
        

# 2.Conjugate Gradient:

def gradient_conjugue( x0 , f ):
    """ Minimisation de la fonction scalaire d'une ou plusieurs variable.
    

    Parameters
    ----------
    x0 : ndarray, shape (n,)
        Point initiale.
    f : callable
         La fonction objective à minimiser.
       
          ``f( x ) -> float``
          
            x`` is an 1-D array with shape (n,)
    

    Returns
    -------
     res : ndarray, shape (n,)
        Le résultat de l'optimisation.

    """
    List = [ ]
    df=nd.Gradient( f )
    d2f=nd.Hessian( f )
    d0=-df(x0)
    n=len(x0)
    hess=d2f(x0)
    for k in range(n):
        List += [ x0 ]
        Td0 = d0.T
        a = ( Td0 @ d0 )/( Td0 @ hess @ d0 )
        x0 = x0 + a * d0
        grd = df(x0)
        hess = d2f(x0)
        b =( grd.T @ hess @ d0 ) / ( Td0 @ hess @ d0 )
        d0 = -grd + b * d0
    res = x0
   
    return  res , List

#==================================================\\Newton methods//================================================

# 1.Newton Method:

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

def methode_Newton( x0 , f , tol = 1e-3 ):
    """ Minimisation de la fonction scalaire d'une ou plusieurs variable.
    

    Parameters
    ----------
    x0 : ndarray, shape (n,)
        Point initiale.
    f : callable
         La fonction objective à minimiser.
       
          ``f( x ) -> float``
          
            x`` is an 1-D array with shape (n,)
    tol : float, optional
         Tolérance pour la terminaison. The default is 1e-3.

    Returns
    -------
     res : ndarray, shape (n,)
        Le résultat de l'optimisation.

    """
    
    df = nd.Gradient( f )
    d2f = nd.Hessian( f )
    n = len( x0 )
    d0 = -np.linalg.inv( d2f( x0 ) ) @ df( x0 )
    g = lambda x , y :  f( x + y * np.linalg.inv( d2f( x0 ) ) @ df( x ) )

    while np.linalg.norm( d0 ) > tol :

        g1 = lambda y : g( x0 , y )
        y = Quasi_Newton( 0 , g1 ) 
        x = x0 - y * d0;
        hess = d2f( x )
        if define_positif( hess )==True:
          d0 = -np.linalg.inv( hess ) @ df( x )
        else:
            va,ve = np.linalg.eig( d2f(x0) )
            eps = min ( va )
            d0 = -np.linalg.inv( eps * np.eye( n ) @ hess ) @ df( x )
        x0 = x 
    res = x0
    return res 



# 2.Armijo Method:

def armijo(phi,alpha0,ita,epsilon,d,grad):
    
    
    approximation = lambda alpha : phi(0)+ epsilon* (d.T @ grad) * alpha

    alpha=alpha0 
    if  phi( alpha )<approximation(alpha) :
      while phi( alpha )<approximation(alpha) :
            ancienne_v =alpha
            alpha=ita*alpha
      return ancienne_v 
    
    
    else :
       while phi( alpha ) >= approximation(alpha) :
            ancienne_v =alpha
            alpha= alpha/ita           
       return ancienne_v

# 3.Davidon Fletcher Powell

def Davidon_Fletcher_Powell( H , d , y , alpha ):
   A= (alpha* (d @ d.T) )/ (d.T @ y )
   B= - ( H @ y ) @ (( H @ y ).T) / ( y.T @ H @ y )
   return H + A + B

# 4.Quasi-Newton with DFP and armijo:


def Quasi_Newton_and_armijo( x0 , f ,tol = 1e-3 ):
  """ Minimisation de la fonction scalaire d'une ou plusieurs variable.
    

    Parameters
    ----------
    x0 : ndarray, shape (n,)
        Point initiale.
    f : callable
         La fonction objective à minimiser.
       
          ``f( x ) -> float``
          
            x`` is an 1-D array with shape (n,)
    tol : float, optional
         Tolérance pour la terminaison. The default is 1e-3.

    Returns
    -------
     res : ndarray, shape (n,)
        Le résultat de l'optimisation.
     List_sol_iteration : list
                          liste des solution

  """
  ita=2
  alpha0 = tol 
  n = len( x0 ) 
  df = nd.Gradient( f )
  grad = df( x0 )
  H = np.eye( n )
  List_sol_iteration= [ x0 ]
  while ( np.linalg.norm( grad ) > tol ):
     d = - H @ grad
     phi= lambda alpha : f(x0+alpha*d)  
     alpha_min  = armijo( phi , alpha0 , ita , tol , d , grad )
     x0 = x0 + alpha_min * d
     List_sol_iteration.append( x0 )
     y = - grad
     grad = df( x0 )  
     y = y + grad
     H = Davidon_Fletcher_Powell( H , d , y , alpha_min )
  res = x0
  return res , List_sol_iteration

#====================================================================================================================
