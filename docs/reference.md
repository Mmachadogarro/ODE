The following module shows the the usage of the numerical solutions. 


import numpy as np 
import matplotlib.pyplot as plt



def func(x,t): 
    """ 
    Diferencial equation to be solved
    
    Parameters
    ----------
    x : numerical values such as floats or doubles
        First  argument 
    t : numerical values such as floats or doubles 
        Second argument 

    Returns
    -------
    output : numerical values such as floats or doubles
        Returns the value of the diferential equation evaluated at (x,t)
    
    Example 
    --------
    --- func(1,np.pi/2) 
    0
    """
    return -x**3 + np.sin(t)


def euler(f,x_0,t): 

    """
    Euler method used to solve diferential equations 

    Parameters
    -----------

    f : numerical values such as floats or doubles
        First argument

    x_0 : numerical values such as floats or doubles
        Second argument 

    t : lists
        Third argument 
    
    Returns 
    --------
    output : lists 
        Returns the values of the solved diferential equation at each point
    
    Example
    -------
   
    t = [ 0.          0.52631579  1.05263158  1.57894737  2.10526316  2.63157895
  3.15789474  3.68421053  4.21052632  4.73684211  5.26315789  5.78947368
  6.31578947  6.84210526  7.36842105  7.89473684  8.42105263  8.94736842
  9.47368421 10.        ]

   --- [ 0.          0.          0.26439534  0.71189381  1.04830651  0.89488942
  0.77464602  0.52141013  0.17502349 -0.28921312 -0.80263945 -0.97897518
 -0.73458315 -0.50879969 -0.16038526  0.30726691  0.81787725  0.97386719
  0.72957635  0.49945695]

        """

    h = t[1]-t[0]
    x = np.zeros(t.size)
    x[0] = x_0 
    for i in range(1,t.size):
        x[i] = x[i-1] + h*f(x[i-1],t[i-1])
    return x



def rk2(f,x_0,t):

 """
    Runge-Kutta2  method used to solve diferential equations 

    Parameters
    -----------

    f : numerical values such as floats or doubles
        First argument

    x_0 : numerical values such as floats or doubles
        Second argument 

    t : lists
        Third argument 
    
    Returns 
    --------
    output : lists 
        Returns the values of the solved diferential equation at each point
    
    Example
    -------
t = [ 0.          0.52631579  1.05263158  1.57894737  2.10526316  2.63157895
  3.15789474  3.68421053  4.21052632  4.73684211  5.26315789  5.78947368
  6.31578947  6.84210526  7.36842105  7.89473684  8.42105263  8.94736842
  9.47368421 10.        ]

---rk2(func,0,t)

[ 0.          0.13691106  0.50040599  0.83221858  0.89696773  0.83638552
  0.684369    0.42791827  0.0377284  -0.46988042 -0.7896376  -0.78706438
 -0.65422559 -0.40234264 -0.0089812   0.49847514  0.79691697  0.78634182
  0.64914835  0.39298319]
    

"""
    h = t[1]-t[0]
    x = np.zeros(t.size)
    x[0] = x_0
    for i in range(1,t.size):
        k_1 = h*f(x[i-1],t[i-1])
        k_2 = h*f(x[i-1] + 0.5*k_1,t[i-1] + h*0.5)
        x[i] = x[i-1] + k_2
    return x 

def rk4(f,x_0,t):

    """
    Runge-Kutta4 method used to solve diferential equations

    Parameters
    -----------
f : numerical values such as floats or doubles
        First argument

    x_0 : numerical values such as floats or doubles
        Second argument

    t : lists
        Third argument

    Returns
    --------
    output : lists
        Returns the values of the solved diferential equation at each point

    Example
    -------
t = [ 0.          0.52631579  1.05263158  1.57894737  2.10526316  2.63157895
  3.15789474  3.68421053  4.21052632  4.73684211  5.26315789  5.78947368
  6.31578947  6.84210526  7.36842105  7.89473684  8.42105263  8.94736842
  9.47368421 10.        ]

---rk4(func,0,t)

[ 0.          0.13505937  0.48487586  0.81979626  0.93102462  0.88144889
  0.72370945  0.46026479  0.06848691 -0.42773679 -0.78084333 -0.83096351
 -0.6993985  -0.43991777 -0.04506765  0.45089026  0.79057875  0.83047637
  0.69376603  0.43014634]
"""

    h = t[1]-t[0]
    x = np.zeros(t.size)
    x[0] = x_0
    for i in range(1,t.size):
        k_1 = h*f(x[i-1],t[i-1])
        k_2 = h*f(x[i-1] + 0.5*k_1,t[i-1] + h*0.5)
        k_3 = h*f(x[i-1] + 0.5*k_2,t[i-1] + h*0.5)
        k_4 = h*f(x[i-1] + k_3,t[i-1] + h)
        x[i] = x[i-1] + (k_1 + 2*k_2 + 2*k_3 + k_4)/6
    return x

